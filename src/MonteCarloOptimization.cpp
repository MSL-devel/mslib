/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#include "MonteCarloOptimization.h"

using namespace MSL;
using namespace std;


MonteCarloOptimization::MonteCarloOptimization(){
	// Set defaults
	selfEnergy              = NULL;
	pairEnergy              = NULL;

	totalNumPositions       = 0;
	initType                = LOWESTSELF;
	numStoredConfigurations = 10;

	responsibleForEnergyTableMemory = false;
	pRng = new RandomNumberGenerator;
	deleteRng = true;
	verbose = false;
}

MonteCarloOptimization::~MonteCarloOptimization(){
	deletePointers();
}

void MonteCarloOptimization::deletePointers() {
	deleteEnergyTables();
	if(deleteRng && pRng) {
		delete pRng;
		pRng = NULL;
	}

}

void MonteCarloOptimization::deleteEnergyTables() {
	if(responsibleForEnergyTableMemory ) {
		if(selfEnergy) {
			delete selfEnergy;
			selfEnergy = NULL;
		}
		if(pairEnergy) {
			delete pairEnergy;
			pairEnergy = NULL;
		}
	}
}



	

void MonteCarloOptimization::readEnergyTable(string _filename){

	/*
	  ENERGY TABLE FILE FORMAT:
	  
	  # COMMENT LINES
	  # SELF TERMS
	  #  POSITION_INDEX ROTAMER_INDEX SELF_ENERGY
	  0 0 -0.5
	  0 1 -0.3
	  1 0  0.1
	  ..

	  # PAIR TERMS
	  #  POSITION_INDEX ROTAMER_INDEX POSITION_INDEX ROTAMER_INDEX PAIR_ENERGY
	  0 0 1 0 -0.2
	  0 0 1 1  0.5
	  ...
	  
	*/

	// This object is now responsible for the energy table memory.

	deleteEnergyTables();

	responsibleForEnergyTableMemory = true;

	selfEnergy = new vector<vector<double> >();
	pairEnergy = new vector<vector<vector<vector<double> > > >();

	ifstream fin;
	string line;
	try {
		fin.open(_filename.c_str());
		if (fin.fail()){
			cerr << "ERROR 8904 in MonteCarloOptimization::readEnergyTable opening file "<<_filename<<endl;
			exit(8904);
		}
		
		bool firstPair = true;
		while (!fin.fail()){
			getline(fin,line);
			if (line[0] == '#') continue;
			if (line.length() <= 1) continue;



			vector<string> toks = MslTools::tokenize(line);

			
			// self Energy Line has 3 numbers on it
			if (toks.size() == 3){
				int posIndex = MslTools::toInt(toks[0]);
				if (selfEnergy->size() < posIndex+1){
					vector<double> tmp;
					selfEnergy->push_back(tmp);
				}

				(*selfEnergy)[posIndex].push_back(MslTools::toDouble(toks[2]));
				
			}



			// pair Energy Line has 5 numbers on it
			/*
			  TODO:
			     Make sure pairEnergy is not a full sized table.
			      it should only be a triangle - lower triangular.
			      resize(i+1) for pairEnergy[i]?
			 */
			if (toks.size() == 5){

				if (firstPair){
					pairEnergy->resize(selfEnergy->size());
					// For each position resize to
					// num rotamers
					for (uint i = 0; i < selfEnergy->size();i++){
						(*pairEnergy)[i].resize((*selfEnergy)[i].size());

						// For each rotamer,
						// resize to num positions in lower triangle
						for (uint j = 0; j < (*selfEnergy)[i].size();j++){
							(*pairEnergy)[i][j].resize(i);

							// For each position resize to num rotamers
							for (uint k = 0; k < (*pairEnergy)[i][j].size();k++){
								(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

							}
						}

					}
				}
				firstPair = false;
				
				// Add to pairEnergy object. Indices pairEnergy[POS1][ROT1][POS2][ROT2] = ENERGY;
				int p1 = MslTools::toInt(toks[0]);
				int r1 = MslTools::toInt(toks[1]);
				int p2 = MslTools::toInt(toks[2]);
				int r2 = MslTools::toInt(toks[3]);
				if(p1 > p2) {
					(*pairEnergy)[p1][r1][p2][r2] = MslTools::toDouble(toks[4]);
				} else if(p2 > p1) {
					(*pairEnergy)[p2][r2][p1][r1] = MslTools::toDouble(toks[4]);
				} else {
					cerr << "WARNING 12343: IGNORING LINE: " << line << endl;
				}

				// Add symmetric entries into the pairEnergy table ****NO NEED****
				// pairEnergy[toInt(toks[2])][toInt(toks[3])][toInt(toks[0])][toInt(toks[1])] = toDouble(toks[4]);
				
			}
		}
	} catch (...){
		cerr << "ERROR 7809 in MonteCarloOptimization::readEnergyTable in file: "<<_filename<<endl;
		exit(7809);
	}


	totalNumPositions = (*selfEnergy).size();
	masks.clear();
	masks.resize(totalNumPositions);
	currentState.clear();
	currentState.resize(totalNumPositions);
	for (uint i = 0; i < masks.size();i++){
		masks[i].resize((*selfEnergy)[i].size());
		currentState[i] = 0;
		for (uint j = 0; j < masks[i].size();j++){
			masks[i][j] = true;
		}
	}

}

void MonteCarloOptimization::addEnergyTable(vector<vector<double> > &_selfEnergy, vector<vector<vector<vector<double> > > > &_pairEnergy){
	/*
	  Pair Energy Table Layout:
	  
	  Example 3 positions [ 2 1 2 ] = # of rotamers at each position


	  ====== LowerTriangular - Used here =====
	  Indices:

	  0 

	  1 -> 0

	       |

	       0 -> 0   1


	 2 -> 0             1
	    
	      |             |
	      
	      0 -> 0 1      0 -> 0 1
	       
	      1 -> 0        1 -> 0
	  

	 so:
	  pairEnergy[0].size() = 0;
	  pairEnergy[1].size() = 1;
	  pairEnergy[1][0].size() = 1;
	  pairEnergy[1][0][0].size() = selfEnergy[0].size();
	  
	  pairEnergy[2].size() = 2;
	  pairEnergy[2][0].size() = 2;
	  pairEnergy[2][0][0].size() = selfEnergy[0].size();
	  pairEnergy[2][0][1].size() = selfEnergy[1].size();
	  pairEnergy[2][1][0].size() = selfEnergy[0].size();
	  pairEnergy[2][1][1].size() = selfEnergy[1].size();



	  ====== UpperTriangular - Not Used =====
	  Indices:

	  0 -> 0          1            
	            
               |          |            
		    	  	     	       
	       0  	  0            
	       		  	     	       
               1 -> 0	  1 -> 0       
			                       
               2 -> 0 1	  2 -> 0 1     


	  1 -> 0

	       |

	       0
	       
	       1

	       2 -> 0 1


	 2 -> 0             1
	    


	  

	 so:
	  pairEnergy[0].size() = selfEnergy[0].size();
	  pairEnergy[0][0].size() = selfEnergy.size();
	  pairEnergy[0][0][0].size() = 0;
	  pairEnergy[0][0][1].size() = selfEnergy[1].size();
	  pairEnergy[0][0][2].size() = selfEnergy[2].size();
	  pairEnergy[0][1][0].size() = 0;
	  pairEnergy[0][1][1].size() = selfEnergy[1].size();
	  pairEnergy[0][1][2].size() = selfEnergy[2].size();
	  

	  pairEnergy[1].size() = selfEnergy[1].size();
	  pairEnergy[1][0].size() = selfEnergy.size();
	  pairEnergy[1][0][0].size() = 0;
	  pairEnergy[1][0][1].size() = 0;
	  pairEnergy[1][0][2].size() = selfEnergy[2].size();
	  
	  
	  pairEnergy[2].size() = 2;
	  pairEnergy[2][0].size() = 0;
	  pairEnergy[2][1].size() = 0;


	  
	 */

	deleteEnergyTables();

	selfEnergy = &_selfEnergy;
	pairEnergy = &_pairEnergy; 
	responsibleForEnergyTableMemory = false;
	
	totalNumPositions = _selfEnergy.size();
	masks.clear();
	masks.resize(totalNumPositions);
	currentState.clear();
	currentState.resize(totalNumPositions);
	for (uint i = 0; i < masks.size();i++){
		masks[i].resize((*selfEnergy)[i].size());
		currentState[i] = 0;
		for (uint j = 0; j < masks[i].size();j++){
			masks[i][j] = true;
		}
	}

	
}


vector<unsigned int> MonteCarloOptimization::runMC(double _startingTemperature, double _endingTemperature, int _scheduleCycles, int _scheduleShape, int _maxRejectionsNumber, int _convergedSteps, double _convergedE){
	
	if (verbose) {
		cout << "===================================" << endl;
		cout << "Run Unbiased Monte Carlo Optimization" << endl;
		cout << endl;
	}


	time_t startMCOtime, endMCOtime;
	double MCOTime;

	time (&startMCOtime);
	
	// Create a starting configuration.. sets initState appropriately and sets currentState to initState
	initialize();

	MonteCarloManager MCMngr(_startingTemperature, _endingTemperature,_scheduleCycles, _scheduleShape, _maxRejectionsNumber, _convergedSteps, _convergedE);
	MCMngr.setRandomNumberGenerator(pRng);
	
	unsigned int cycleCounter = 0;
	unsigned int moveCounter = 0;

	vector<unsigned int> bestState = initState;
	double bestEnergy = getStateEnergy(bestState);

	vector<unsigned int> prevStateVec = bestState;
	vector<unsigned int> stateVec;

	
	MCMngr.setEner(bestEnergy);
	string state = getRotString();
	configurationMap[state] = bestEnergy;
	sampledConfigurations.push(pair<double,string>(configurationMap[state],state));				

	while (!MCMngr.getComplete()) {
		// make atleast 1 move and update the current state
		stateVec = moveRandomState(pRng->getRandomInt(totalNumPositions-1) + 1);  // random number between 1 and totalNumPositions
		//stateVec = moveRandomState();  
		if (verbose) {
			for (int i = 0; i < stateVec.size(); i++){
				cout << stateVec[i] << ",";
			}
			cout << endl;
		}
		double oligomerEnergy = getStateEnergy(stateVec);

		if (verbose) {
			cout << "MCO [" << cycleCounter << "]: ";
		}

		if(oligomerEnergy < bestEnergy) {
			bestEnergy = oligomerEnergy;
			bestState = stateVec;
		}
		if (!MCMngr.accept(oligomerEnergy)) {
			setCurrentState(prevStateVec);
			if (verbose) {
				cout << "MCO: State not accepted, E=" << oligomerEnergy << "\n";
			}
		} else {
			prevStateVec = stateVec;
			if (verbose) {
				cout << "MCO: State accepted, E=" << oligomerEnergy << "\n";
			}
			// Update configurationMap with energy
			map<string,double>::iterator it;
			string rotStr = getRotString();
			it = configurationMap.find(rotStr);
			if (it == configurationMap.end()){
				configurationMap[rotStr] = oligomerEnergy;

				// Store configuration ...
				if (sampledConfigurations.size() >= numStoredConfigurations){
					if (oligomerEnergy < sampledConfigurations.top().first){
						// Remove highest energy , then add
						sampledConfigurations.pop();
						sampledConfigurations.push(pair<double,string>(oligomerEnergy,rotStr));				
					}
				} else {
					sampledConfigurations.push(pair<double,string>(oligomerEnergy,rotStr));				
				}
			}

		}
		cycleCounter++;
		moveCounter++;
	}

	time (&endMCOtime);
	MCOTime = difftime (endMCOtime, startMCOtime);
	if (verbose) {
		cout << endl;
		cout << "Best MCO Energy: " << bestEnergy << endl;
		cout << "MCO Time: " << MCOTime << " seconds" << endl;
		cout << "===================================" << endl;
	}
	return bestState;
}

void MonteCarloOptimization::printSampledConfigurations(){
	while (!sampledConfigurations.empty()){
		cout << sampledConfigurations.top().first <<" "<<sampledConfigurations.top().second<<endl;
		sampledConfigurations.pop();
	}
}


string MonteCarloOptimization::getRotString(){
	string result = "";

	for (uint i = 0; i < totalNumPositions;i++){
		result += MslTools::intToString(currentState[i])+":";
	}

	return result;
}
string MonteCarloOptimization::getRotString(int _pos, int _rot){

	string result = "";

	for (uint i = 0; i < totalNumPositions;i++){

		if (i == _pos){
			result += MslTools::intToString(_rot)+":";
		} else {
			result += MslTools::intToString(currentState[i])+":";
		}
	}

	//cout << "rotString: "<<result<<endl;

	return result;
}
void MonteCarloOptimization::initialize(){
	
	// Randomly select a starting state..
	if (initType == RANDOM) {
		if(verbose) {
			cout << "Initializing MCO to random state" << endl;
		}
		getRandomState();	
	} else if (initType == LOWESTSELF){
		if(verbose) {
			cout << "Initializing MCO to lowest self energy state" << endl;
		}
		// Pick lowest self energy state...
		for (uint i = 0; i < totalNumPositions;i++){
			double energy = MslTools::doubleMax;
			int rot       = pRng->getRandomInt((*selfEnergy)[i].size()-1);
			for (uint j = 0; j < (*selfEnergy)[i].size();j++){
				if (masks.size() != 0 && !masks[i][j]){
					continue;
				}


				if ((*selfEnergy)[i][j] < energy){
					energy = (*selfEnergy)[i][j];
					rot = j;
				}
			}
			selectRotamer(i,rot);
		}
	} else if (initType == QUICKSCAN) {    
	// Quickly scan:
	//    Random start. 
	//    Go through each position (highest energy first) and best energy lowering rotamer

		if(verbose) {
			cout << "Initializing MCO using quickscan" << endl;
		}
		vector<pair<int,double> > energies;
		for (uint pos = 0; pos < totalNumPositions;pos++){


			int rot = pRng->getRandomInt((*selfEnergy)[pos].size()-1);
			while (masks.size() != 0 && !masks[pos][rot]){
				rot = pRng->getRandomInt((*selfEnergy)[pos].size()-1);
			}

			selectRotamer(pos,rot);
			
			double energy = getEnergy(pos,rot);

			energies.push_back(pair<int,double>(pos,energy));
		}

		std::sort(energies.begin(),energies.end(),MslTools::sortPairIntDoubleDecending);

		// Now get lowest energy rotamer for each position
		for (uint i = 0; i < energies.size();i++){

			double minEnergy = MslTools::doubleMax;
			int    minRotamer = MslTools::intMax;
			for (uint j = 0; j < (*selfEnergy)[energies[i].first].size();j++){
				if (masks.size() != 0 && !masks[energies[i].first][j]){
					continue;
				}

				double energy = getEnergy(energies[i].first,j);
				if (energy < minEnergy){
					minEnergy = energy;
					minRotamer = j;
				}
			}

			selectRotamer(energies[i].first,minRotamer);

		}
	} else 	if (initType == USERDEF){
		// initState should have been already updated
		if(verbose) {
			cout << "Initializing MCO using userdefined state" << endl;
		}
		currentState = initState;	
		return;
	} else {
		cerr << "ERROR 2454 Unknown initType in MonteCarloOptimization::initialize() - initializing randomly" << endl;	
		getRandomState();	
	}

	initState = currentState;	
	return;
	
}

void MonteCarloOptimization::selectRotamer(int _pos, int _rot){
	//cout << "POS,ROT: "<<_pos<<","<<_rot<<endl;
	if(_pos < currentState.size()) {
	// What if the state is masked? - DONT CARE?
		currentState[_pos] = _rot;
	} else {
		cerr << "ERROR 2456 _pos out of range[0," << currentState.size() << ")" << endl;
	}
}
double MonteCarloOptimization::getEnergy(int _pos, int _rot){

        double energy = 0.0;
        energy += (*selfEnergy)[_pos][_rot];

        for (uint i = 0;i <currentState.size();i++){

                int pos2 = i;
                int rot2 = currentState[i];

                if (_pos == pos2) continue;

                // IF _pos is LINKED THEN WHAT?
                // Then if pos2 is linked to _pos, we need to use _rot instead of rot2?


                if (_pos > pos2){
                        energy += (*pairEnergy)[_pos][_rot][pos2][rot2];
                } else {
                        energy += (*pairEnergy)[pos2][rot2][_pos][_rot];
                }
        }

        return energy;
}

double MonteCarloOptimization::getStateEnergy(vector<unsigned int> _states) {
	double energy = 0.0;
	// Check that the states are indeed valid
	if(_states.size() != selfEnergy->size()) {
		cerr << "ERROR: The number of positions doesnot agree with selfEnergy Table in double MonteCarloOptimization::getStateEnergy(vector<unsigned int> _states) " << endl;
		return energy;
	}
	for (uint i = 0 ;i < _states.size();i++){
		
		//cout << "Adding self: "<<(*selfEnergy)[i][currentState[i]]<<endl;
		if( _states[i] < (*selfEnergy)[i].size()) {
			energy += (*selfEnergy)[i][_states[i]];
		} else {
			cerr << "ERROR: SelfEnergyTable for Position " << i << " doesnot contain rotamer " << _states[i] << " in double MonteCarloOptimization::getStateEnergy(vector<unsigned int> _states) " << endl;
			return energy;
		}
		for (uint j = 0; j < i;j++){
			//cout << "Adding Position "<<i<<" to "<<j<<" which is "<<(*pairEnergy)[i][currentState[i]][j][currentState[jxo]]<<endl;
			
			//energy += (*pairEnergy)[j][currentState[j]][i][currentState[i]];
			energy += (*pairEnergy)[i][_states[i]][j][_states[j]];
		}
	}

	return energy;
}

double MonteCarloOptimization::getStateEnergy(){
	
	double energy = 0.0;
	for (uint i = 0 ;i < totalNumPositions;i++){
		
		//cout << "Adding self: "<<(*selfEnergy)[i][currentState[i]]<<endl;
		energy += (*selfEnergy)[i][currentState[i]];
		for (uint j = 0; j < i;j++){
			//cout << "Adding Position "<<i<<" to "<<j<<" which is "<<(*pairEnergy)[i][currentState[i]][j][currentState[jxo]]<<endl;
			
			//energy += (*pairEnergy)[j][currentState[j]][i][currentState[i]];
			energy += (*pairEnergy)[i][currentState[i]][j][currentState[j]];

		}
	}

	return energy;
}

void MonteCarloOptimization::printMe(bool _selfOnly){
	fprintf(stdout,"Self terms:\n");
	for (uint i = 0; i < (*selfEnergy).size();i++){
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){

			fprintf(stdout, "    %4d %4d %8.3f", i, j, (*selfEnergy)[i][j]);
			// alive Rotamers	
			if (masks[i][j]) {
				fprintf(stdout, " **** ");
			}
			fprintf(stdout,"\n");


		}

	}

	if (_selfOnly){
		return;
	}

	fprintf(stdout,"Pair terms:\n");
	for (uint i = 0; i < (*pairEnergy).size();i++){
		for (uint j = 0; j < (*pairEnergy)[i].size();j++){
			//for (uint k = i+1 ; k < (*pairEnergy)[i][j].size();k++){	
			for (uint k = 0 ; k < (*pairEnergy)[i][j].size();k++){	
				if (i == k){
					continue;
				}
				for (uint l = 0 ; l < (*pairEnergy)[i][j][k].size();l++){	
					fprintf(stdout, "    %4d %4d %4d %4d %8.3f", i, j, k, l, (*pairEnergy)[i][j][k][l]);

					// alive Rotamers	
					if (masks[i][j] && masks[k][l]) {
						fprintf(stdout, " **** ");
					}
					fprintf(stdout,"\n");
				}
			}
		}

	}
}

void MonteCarloOptimization::setInitializationState(vector<unsigned int>& _state) {
	initType = USERDEF;
	initState = _state;
}
void MonteCarloOptimization::setInitializationState(INITTYPE _type, std::string _userDef){
	initType = _type;
	if(_userDef != "") {
		std::vector<string> toks = MslTools::tokenize(_userDef,":");
		if(toks.size() != selfEnergy->size()) {
			cerr << "ERROR 2456 no: of states in _userDef does not match with EnergyTables" << endl;
			return;
		}
		initState.resize(toks.size());
		for(int i = 0; i < toks.size(); i++) {
			initState[i] = MslTools::toInt(toks[i]);
		}
	}
}
vector<unsigned int> MonteCarloOptimization::getRandomState() {
	for (int i=0; i<totalNumPositions; i++) {
		currentState[i] = selectRandomStateAtPosition(i);
	}
	return currentState;
}
vector<unsigned int> MonteCarloOptimization::moveRandomState()  {
	return moveRandomState(1);
}

vector<unsigned int> MonteCarloOptimization::moveRandomState(unsigned int _numOfMoves)  {
	/****************************************************
	 *  this function selects one or more moves 
	 ****************************************************/

	if (_numOfMoves > currentState.size()) {
		cerr << "WARNING 3205: the number of step is larger than the number of positions in vector<unsigned int> MonteCarloOptimization::moveRandomState(unsigned int _numOfMoves)" << endl;
		_numOfMoves = currentState.size();
	}
	
	/*********************************************
	 *  Pick a position
	 *********************************************/
	vector<bool> alreadySelected(currentState.size(), false);
	for (unsigned int i=0; i<_numOfMoves; i++) {
		vector<double> residualP;
		double sumP = 0.0;
		for (unsigned int j=0; j<currentState.size(); j++) {
			if (alreadySelected[j]) {
				// if the position was already choosen, no prob to it
				residualP.push_back(0);
			} else {
				// all other positions have the same probability 
				residualP.push_back(1);
			}
			sumP += residualP.back();
		}


		if (sumP == 0.0) {
			// there is no available move, stay there
			continue;
		}
	
		pRng->setDiscreteProb(residualP);
		unsigned int randomPos = pRng->getRandomDiscreteIndex();
		currentState[randomPos] = selectRandomStateAtPosition(randomPos);
		alreadySelected[randomPos] = true;
	}

	return currentState;
}

int MonteCarloOptimization::selectRandomStateAtPosition(int _position) const {
	if (_position >= currentState.size()) {
		cerr << "ERROR 7210: position " << _position << " out of range in int MonteCarloOptimization::selectRandomStateAtPosition(int _position) const" << endl;
		return 0;
	}

	if ((*selfEnergy)[_position].size() == 1) {
		return 0;
	}

	vector<double> residualP;
	double sumP = 0.0;
	for (int i=0; i<(*selfEnergy)[_position].size(); i++) {
		if (i == currentState[_position] || !masks[_position][i]) {
			residualP.push_back(0.0);
		} else {
			residualP.push_back(1.0);
		}
		sumP = residualP.back();
	}

	if (sumP == 0.0) {
		// no move avalable, stay there
		return currentState[_position];
	}

	pRng->setDiscreteProb(residualP);
	return pRng->getRandomDiscreteIndex();
}

