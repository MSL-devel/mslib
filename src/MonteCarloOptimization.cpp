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

	totalNumRotamers        = 0;
	totalNumPositions       = 0;
	cycles                  = 1000;

	annealType              = LIN_TEMP_ANNEAL;
	startTemp               = 1000;
	endTemp                 = 1;
	temp                    = 1000;
	initType                = LOWESTSELF;
	numStoredConfigurations = 10;
	randomSeed              = 0;

	responsibleForEnergyTableMemory = false;
	verbose = false;

	masks.clear();
}

MonteCarloOptimization::~MonteCarloOptimization(){

	
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
						// resize to num positions
						for (uint j = 0; j < (*selfEnergy)[i].size();j++){
							(*pairEnergy)[i][j].resize((*selfEnergy).size());

							// For each position resize to num rotamers
							for (uint k = 0; k < (*selfEnergy).size();k++){
								(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

							}
						}

					}
				}
				firstPair = false;
				
				// Add to pairEnergy object. Indices pairEnergy[POS1][ROT1][POS2][ROT2] = ENERGY;
				(*pairEnergy)[MslTools::toInt(toks[0])][MslTools::toInt(toks[1])][MslTools::toInt(toks[2])][MslTools::toInt(toks[3])] = MslTools::toDouble(toks[4]);

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
	rotamerSelection.clear();
	rotamerSelection.resize(totalNumPositions);
	for (uint i = 0; i < masks.size();i++){
		masks[i].resize((*selfEnergy)[i].size());
		rotamerSelection[i] = 0;
		totalNumRotamers += masks[i].size();
		for (uint j = 0; j < masks[i].size();j++){
			masks[i][j] = false;
		}
	}

}

void MonteCarloOptimization::addEnergyTable(vector<vector<double> > &_selfEnergy, vector<vector<vector<vector<double> > > > &_pairEnergy){
	/*
	  Pair Energy Table Layout:
	  
	  Example 3 positions [ 2 1 2 ] = # of rotamers at each position


	  ====== LowerTriangular =====
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



	  ====== UpperTriangular =====
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

	selfEnergy = &_selfEnergy;


	pairEnergy = new vector<vector<vector<vector<double> > > >();
	pairEnergy->resize(selfEnergy->size());

	totalNumPositions = _selfEnergy.size();

	// For each position resize to
	// num rotamers
	for (uint i = 0; i < selfEnergy->size();i++){
		(*pairEnergy)[i].resize((*selfEnergy)[i].size());
		totalNumRotamers += (*selfEnergy)[i].size();

		// For each rotamer,
		// resize to num positions
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){
			(*pairEnergy)[i][j].resize((*selfEnergy).size());

			// For each position resize to num rotamers
			for (uint k = 0; k < (*selfEnergy).size();k++){
				(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

			}
		}

	}


	// Convert to upper triangular
	for (uint pos1 = 0 ; pos1 < _pairEnergy.size();pos1++){
		for (uint pos2 = pos1+1 ; pos2 < _pairEnergy.size();pos2++){
			for (uint rot2 = 0 ; rot2 < _pairEnergy[pos2].size();rot2++){
				for (uint rot1 = 0 ; rot1 < _pairEnergy[pos2][rot2][pos1].size();rot1++){
					(*pairEnergy)[pos1][rot1][pos2][rot2] =  _pairEnergy[pos2][rot2][pos1][rot1];
				}
			}
		}
	}


	
	masks.clear();
	masks.resize(totalNumPositions);
	rotamerSelection.clear();
	rotamerSelection.resize(totalNumPositions);
	for (uint i = 0; i < masks.size();i++){
		masks[i].resize((*selfEnergy)[i].size());
		rotamerSelection[i] = 0;
		for (uint j = 0; j < masks[i].size();j++){
			masks[i][j] = false;
		}
	}

	
}


void MonteCarloOptimization::runMC(){
	
	// Create a starting configuration..
	initialize();

	// Tabulate Energies
	double totalEnergy = 0.0;
	double newEnergy   = 0.0;
	totalEnergy = getTotalEnergy();
	newEnergy   = totalEnergy;

	configurationMap[getRotString()] = totalEnergy;
	sampledConfigurations.push(pair<double,string>(totalEnergy,getRotString()));				

	//cout << "TotalEnergy: "<<totalEnergy<<endl;

	
	
	// Gas-constant
	double R =  1.9872e-3;

	// Loop
	for (currentStep = 0 ; currentStep < cycles ;currentStep++){

		
		// Select Random Change
		vector<int> nextRotamer = getRandomRotamer();

		// get "new" energy : total - current + next;
		newEnergy = totalEnergy - getEnergy(nextRotamer[0],rotamerSelection[nextRotamer[0]]) + getEnergy(nextRotamer[0],nextRotamer[1]);

		// Linked Positions need to considered here both in the new rotamer and the previous!
		
		// Metropolis..
		bool accept = true;
		if (newEnergy > totalEnergy) {

			// Anneal the temperature according to schedule..
			annealTemperature(startTemp, endTemp, currentStep, cycles);

			// Reject based on Boltzmann
			double p    = exp( -(newEnergy-totalEnergy) / (R*temp) );
			double rand = rng.getRandomDouble();
			if (rand > p){
				accept = false;
			}
			
		}

		if (accept){

			// Set rotamer active..
			selectRotamer(nextRotamer[0],nextRotamer[1]);

			// Linked positions set rotamers active.

			// Keep track of totalEnergy
			totalEnergy = newEnergy;

			// Update configurationMap with energy
			map<string,double>::iterator it;
			string rotStr = getRotString();
			it = configurationMap.find(rotStr);
			if (it == configurationMap.end()){
				configurationMap[rotStr] = newEnergy;

				// Store configuration ...
				if (sampledConfigurations.size() >= numStoredConfigurations){
					if (totalEnergy < sampledConfigurations.top().first){
						// Remove highest energy , then add
						sampledConfigurations.pop();
						sampledConfigurations.push(pair<double,string>(totalEnergy,getRotString()));				
					}
				} else {
					sampledConfigurations.push(pair<double,string>(totalEnergy,getRotString()));				
				}
			}


			



		}


		// Update configurationMap with energy
		//configurationMap[getRotString(nextRotamer[0],nextRotamer[1])] = newEnergy;
		
	}

	//cout << "OptimizedEnergy: "<<totalEnergy<<" "<<getTotalEnergy()<<endl;


	
}

void MonteCarloOptimization::printSampledConfigurations(){
	while (!sampledConfigurations.empty()){
		cout << sampledConfigurations.top().first <<" "<<sampledConfigurations.top().second<<endl;
		sampledConfigurations.pop();
	}
}

vector<int> MonteCarloOptimization::getRandomRotamer(){

	vector<int> rot;

	int nextPosition = 0;
	int nextRotamer  = 0;

	// Try 10,000 times to get a random position and rotamer that has not been tried before..
	int i = 0;
	while (i < 10000){

		// increment counter
		i++;

		nextPosition = rng.getRandomInt(totalNumPositions-1);
		nextRotamer  = rng.getRandomInt((*selfEnergy)[nextPosition].size()-1);

		// Check for already sampled configuration..
		map<string,double>::iterator it;


		string rotString = getRotString(nextPosition,nextRotamer);

		// Check our configurationMap
		/*
		it = configurationMap.find(rotString);
		if (it == configurationMap.end()){
			configurationMap[rotString] = MslTools::doubleMax;
			break;
		} 
		*/

		
		if (inputMasks.size() != 0 && !inputMasks[nextPosition][nextRotamer]){
			continue;
		}

		if (nextRotamer != rotamerSelection[nextPosition]){
			break;
		}


	}

	// If we can't find a unique configuration
	//    This can mean 2 things:
	//    1. We have completely sampled configuational space
	//    2. We are stuck and need to change 2 rotamers to get to a new configuration.
	/*
	if (i == 10000){
		cout << "Possible completion of configurational sampling..."<<endl;
		map<string,double>::iterator it;
		for (it = configurationMap.begin();it != configurationMap.end();it++){
			cout << "RotString: "<<it->first<<" "<<it->second<<endl;
		}
	}
	*/
	if (inputMasks.size() != 0 && !inputMasks[nextPosition][nextRotamer]){
		cerr << "ERROR MonteCarloOptimization::getRandomRotamer() TRIED TO PICK A ROTAMER 10000 times and only one I found was a EXCLUDED rotamer!!!!!!!!!!!!"<<endl;
		exit(3455);

	}
	// Store nextPostion,nextRotamer..
	rot.push_back(nextPosition);
	rot.push_back(nextRotamer);


	return rot;

	
}

string MonteCarloOptimization::getRotString(){
	string result = "";

	for (uint i = 0; i < totalNumPositions;i++){
			result += MslTools::intToString(rotamerSelection[i])+":";
	}

	return result;
}
string MonteCarloOptimization::getRotString(int _pos, int _rot){

	string result = "";

	for (uint i = 0; i < totalNumPositions;i++){

		if (i == _pos){
			result += MslTools::intToString(_rot)+":";
		} else {
			result += MslTools::intToString(rotamerSelection[i])+":";
		}
	}

	//cout << "rotString: "<<result<<endl;

	return result;
}
void MonteCarloOptimization::initialize(){
	
	// Set seed.
	//rng.setRNGType("knuthran2002");

	if (randomSeed == 0){
		rng.setTimeBasedSeed();
		randomSeed = rng.getSeed();
	} else {
		rng.setSeed(randomSeed);
	}
	
	// Randomly select a starting state..
	if (initType == RANDOM) {
		

		for (uint pos = 0; pos < totalNumPositions;pos++){

			int rot = rng.getRandomInt((*selfEnergy)[pos].size()-1);

			while (inputMasks.size() != 0 && !inputMasks[pos][rot]){
				rot = rng.getRandomInt((*selfEnergy)[pos].size()-1);
			}
			selectRotamer(pos,rot);
		}
	}

	// Pick lowest self energy state...
	if (initType == LOWESTSELF){
		for (uint i = 0; i < totalNumPositions;i++){
			double energy = MslTools::doubleMax;
			int rot       = rng.getRandomInt((*selfEnergy)[i].size()-1);
			for (uint j = 0; j < (*selfEnergy)[i].size();j++){
				if (inputMasks.size() != 0 && !inputMasks[i][j]){
					continue;
				}


				if ((*selfEnergy)[i][j] < energy){
					energy = (*selfEnergy)[i][j];
					rot = j;
				}
			}

			selectRotamer(i,rot);
		}
	}

	// Quickly scan:
	//    Random start. 
	//    Go through each position (highest energy first) and best energy lowering rotamer
	if (initType == QUICKSCAN){    

		vector<pair<int,double> > energies;
		for (uint pos = 0; pos < totalNumPositions;pos++){


			int rot = rng.getRandomInt((*selfEnergy)[pos].size()-1);
			while (inputMasks.size() != 0 && !inputMasks[pos][rot]){
				rot = rng.getRandomInt((*selfEnergy)[pos].size()-1);
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
				if (inputMasks.size() != 0 && !inputMasks[energies[i].first][j]){
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
	}

	if (initType == USERDEF){
		vector<string> toks = MslTools::tokenize(initConf,":");
		if (toks.size() != totalNumPositions){
			cerr << "ERROR 2453 Number of tokens: "<<toks.size()<<" Number of positions: "<<totalNumPositions<<endl;
			exit(2453);
		}

		for (uint i =0; i < totalNumPositions;i++){

			// Should I check here for inputMasks ?
			//  What if a user asked for something that has been 
			//   excluding via the inputMasks?
			// Now we allow user to do what he/she wants..
			selectRotamer(i,MslTools::toInt(toks[i]));
		}
	}
	
}

void MonteCarloOptimization::selectRotamer(int _pos, int _rot){

	//cout << "POS,ROT: "<<_pos<<","<<_rot<<endl;
	rotamerSelection[_pos];

	masks[_pos][rotamerSelection[_pos]] = false;
	rotamerSelection[_pos]              = _rot;
	masks[_pos][_rot]                   = true;
	
}

double MonteCarloOptimization::getEnergy(int _pos, int _rot){

	double energy = 0.0;
	energy += (*selfEnergy)[_pos][_rot];

	for (uint i = 0;i <rotamerSelection.size();i++){

		int pos2 = i;
		int rot2 = rotamerSelection[i];
		
		if (_pos == pos2) continue;
		
		// IF _pos is LINKED THEN WHAT?
		// Then if pos2 is linked to _pos, we need to use _rot instead of rot2?


		if (_pos < pos2){
			energy += (*pairEnergy)[_pos][_rot][pos2][rot2];
		} else {
			energy += (*pairEnergy)[pos2][rot2][_pos][_rot];
		}

		/*
		if (_pos > pos2){
			energy += (*pairEnergy)[_pos][_rot][pos2][rot2];
		} else {
			energy += (*pairEnergy)[pos2][rot2][_pos][_rot];
		}
		*/
	}

	return energy;
}
double MonteCarloOptimization::getTotalEnergy(){
	
	double energy = 0.0;
	for (uint i = 0 ;i < totalNumPositions;i++){
		
		//cout << "Adding self: "<<(*selfEnergy)[i][rotamerSelection[i]]<<endl;
		energy += (*selfEnergy)[i][rotamerSelection[i]];
		for (uint j = i+1; j < totalNumPositions;j++){
			//cout << "Adding Position "<<i<<" to "<<j<<" which is "<<(*pairEnergy)[i][rotamerSelection[i]][j][rotamerSelection[jxo]]<<endl;
			
			//energy += (*pairEnergy)[j][rotamerSelection[j]][i][rotamerSelection[i]];
			energy += (*pairEnergy)[i][rotamerSelection[i]][j][rotamerSelection[j]];

		}
	}

	return energy;
}


void MonteCarloOptimization::printMe(bool _selfOnly){
	fprintf(stdout,"Self terms:\n");
	for (uint i = 0; i < (*selfEnergy).size();i++){
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){

			fprintf(stdout, "    %4d %4d %8.3f", i, j, (*selfEnergy)[i][j]);
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

					if (masks[i][j] && masks[k][l]) {
						fprintf(stdout, " **** ");
					}
					fprintf(stdout,"\n");
				}
			}
		}

	}
}


void MonteCarloOptimization::annealTemperature(double initialTemp, double finalTemp, int step, int totalsteps){


	double FREQ;
	double scaleStep;
	double cycleNumber;
	double expFactor;
	double expVal;
	cycleNumber = expVal = expFactor =0.0;
	int numStepsInCycle;
	int lastStart      ;
	double zeroingFactor  ;

	switch (annealType) {
	case LIN_TEMP_ANNEAL:
		setCurrentTemp((initialTemp - (step * ((initialTemp - finalTemp) / (double)totalsteps))));  
		break;
	case EXP_TEMP_ANNEAL:
		setCurrentTemp((initialTemp * pow((finalTemp/initialTemp), ((double)step/(double)totalsteps))));
		break;
	case SAWTOOTH_TEMP_ANNEAL:
		FREQ = 5;// a parameter to set from outside the object..
		scaleStep   = step / FREQ;
		setCurrentTemp((initialTemp/2) * ((-2*(scaleStep - floor(scaleStep +.5)))+1));
		break;

	case EXPCYCLE_TEMP_ANNEAL:

		/*
		  FREQ = 5; // a parameter to set from outside the object..
		  cycleNumber = floor(double(step) / double(FREQ)) ; 
		  expFactor   = double(FREQ) / 5;
		  expVal  =  (step - (FREQ * cycleNumber)) / expFactor; 
		  setCurrentTemp(startTemp* exp(-expVal)); 

		*/
		// 
		numStepsInCycle = int(cycles / numAnnealCycles);
		lastStart       = (int(step / numStepsInCycle) * numStepsInCycle)+1;
		if (step % numStepsInCycle == 0){
			lastStart = step;
		}
		zeroingFactor   = double(numStepsInCycle)/10 + double(numStepsInCycle)/20; // seems to work ok for deciding when we are close to zero

		setCurrentTemp(initialTemp * exp( - (step-lastStart) / zeroingFactor));
		//cout << "FORMULA: "<<initialTemp<<" * exp( - ("<<step<<" - "<< lastStart<<") / "<<zeroingFactor<<")"<<endl;
		//cout << "TEMP: "<<step<<" "<<temp<<endl;

		break;
	case NO_ANNEAL:
		break;
	default:
		cerr << "No annealing will be done: " << annealType << endl;
	}	

}
