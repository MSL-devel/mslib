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

#include "PairwiseEnergyCalculator.h"
#include <omp.h>

PairwiseEnergyCalculator::PairwiseEnergyCalculator(string _charmmParameterFile){
	pairwiseEnergy = new AtomicPairwiseEnergy(_charmmParameterFile);
	pairwiseEnergy->setEnergyByType(true);
	
}

PairwiseEnergyCalculator::~PairwiseEnergyCalculator(){ 
	delete(pairwiseEnergy);
}

double PairwiseEnergyCalculator::calculateTotalEnergy(System &_sys){


	double energy = 0.0;
	vector<double> selfEnergy(_sys.positionSize(),0.0);
	vector<double> templateEnergy(_sys.positionSize(),0.0);
        vector<vector<double> > pairEnergy(_sys.positionSize(),vector<double>(_sys.positionSize(),0.0));
	map<int,int>::iterator it;
	
	uint i;
	int tid;
	tid =0;
	//#pragma omp parallel private(i,tid) 
	//{

	//#pragma omp for schedule(dynamic)
	for (i = 0; i < _sys.positionSize();i++){

		Position &p1 = _sys.getPosition(i);

		
		for (uint j = i; j < _sys.positionSize();j++){
			Position &p2 = _sys.getPosition(j);


			//cout << "POSITIONS: "<<i<<" "<<j<<endl;
			// Compute self and template terms..
			if (i == j){
				selfEnergy[i] += pairwiseEnergy->calculateSelfEnergy(_sys, i, 0);
				templateEnergy[i] += pairwiseEnergy->calculateTemplateEnergy(_sys, i, 0);
			} else {

				// Only compute pair energies between variable positions..
				if (p1.getTotalNumberOfRotamers() > 1 && p2.getTotalNumberOfRotamers() > 1){
					pairEnergy[i][j] += pairwiseEnergy->calculatePairEnergy(_sys, i, 0, j, 0); 
				}
			}

		}

	}
	//}

	for (i = 0; i < selfEnergy.size();i++){
		energy += templateEnergy[i] + selfEnergy[i];  
		
		for (uint j = 0; j < pairEnergy.size();j++){
			energy += pairEnergy[i][j];
		}
	}



	return energy;
	
}



double PairwiseEnergyCalculator::calculateStateEnergy(System &_sys,vector<unsigned int> &_stateVector){


	double energy = 0.0;
	vector<double> selfEnergy(_sys.positionSize(),0.0);
	vector<double> templateEnergy(_sys.positionSize(),0.0);
        vector<vector<double> > pairEnergy(_sys.positionSize(),vector<double>(_sys.positionSize(),0.0));
	map<int,int>::iterator it;
	
	uint i;
	int tid;
	tid =0;

	map<int,int> rotamersOfVariablePositions;
	int rotamerIndex = 0;
	for (i = 0; i < _sys.positionSize();i++){
		if (_sys.getPosition(i).getTotalNumberOfRotamers() > 1){
			//cout << "Position : "<<_sys.getPosition(i).getChainId()<< " "<<_sys.getPosition(i).getResidueName()<<" "<<_sys.getPosition(i).getResidueNumber()<<" is rotamer "<<_stateVector[rotamerIndex]<<endl;
			rotamersOfVariablePositions[i] = _stateVector[rotamerIndex];
			rotamerIndex++;
		}
	}
	
	//#pragma omp parallel private(i,tid) 
	//{

	//#pragma omp for schedule(dynamic)
	for (i = 0; i < _sys.positionSize();i++){

		Position &p1 = _sys.getPosition(i);

		int rotamer1 = 0;
		it = rotamersOfVariablePositions.find(i);
		if (it != rotamersOfVariablePositions.end()){
			rotamer1 = it->second;
		}
		for (uint j = i; j < _sys.positionSize();j++){
			Position &p2 = _sys.getPosition(j);

			int rotamer2 = 0;
			it = rotamersOfVariablePositions.find(j);
			if (it != rotamersOfVariablePositions.end()){
				rotamer2 = it->second;
			}

			//cout << "POSITIONS: "<<i<<" "<<j<<endl;
			// Compute self and template terms..
			if (i == j){
				double s = pairwiseEnergy->calculateSelfEnergy(_sys, i, rotamer1);
				double t = pairwiseEnergy->calculateTemplateEnergy(_sys, i, rotamer1);
				selfEnergy[i]     += s;
				templateEnergy[i] += t;

				//fprintf(stdout, "Self[%4d] = %8.3f\n",i,s);
				//fprintf(stdout, "Template[%4d] = %8.3f\n\n",i,t);
			} else {

				// Only compute pair energies between variable positions..
				if (p1.getTotalNumberOfRotamers() > 1 && p2.getTotalNumberOfRotamers() > 1){
					double p = pairwiseEnergy->calculatePairEnergy(_sys, i, rotamer1, j, rotamer2); 
					pairEnergy[i][j] += p;

					//fprintf(stdout, "Pair[%4d][%4d] = %8.3f\n\n",i,j,p);
				}
			}
			

		}

	}
	//}

	for (i = 0; i < selfEnergy.size();i++){
		energy += templateEnergy[i] + selfEnergy[i];  
		
		for (uint j = 0; j < pairEnergy.size();j++){
			energy += pairEnergy[i][j];
		}
	}



	return energy;
	
}


double PairwiseEnergyCalculator::calculateEnergyTable(System &_sys){


	double energy = 0.0;


	// Size of tables will be total positions - slave positions (that leaves unlinked and master positions)
	int sysSize = _sys.positionSize()- _sys.slavePositionSize();
	selfEnergy.clear();
	selfEnergy.resize(sysSize);

	templateEnergy.clear();
	templateEnergy.resize(sysSize);

	pairEnergy.clear();
	pairEnergy.resize(sysSize);
	
	uint i;
	int tid;
	tid =0;

	//#pragma omp parallel private(i,tid) 
	//{

	//#pragma omp for schedule(dynamic)
	int posIndexI = -1;
	for (i = 0; i < _sys.positionSize();i++){

		// Get position
		Position &p1 = _sys.getPosition(i);

		// Skip linked positions that are not the "master"
		if (p1.getLinkedPositionType() == Position::SLAVE) { continue; }

		// Index of non-linked and master positions
		posIndexI++;
		
		cout << p1.getChainId()<<" "<<p1.getResidueNumber()<< " "<<p1.getResidueName()<< " is index "<<posIndexI<<" ("<<p1.getTotalNumberOfRotamers()<<")"<<endl;

		//Resize vectors at position posIndexI
		pairEnergy[posIndexI].resize(p1.getTotalNumberOfRotamers());
		selfEnergy[posIndexI].resize(p1.getTotalNumberOfRotamers());
		templateEnergy[posIndexI].resize(p1.getTotalNumberOfRotamers());

		// For each rotamer at position p1
		for (uint r1 = 0; r1 < p1.getTotalNumberOfRotamers();r1++){


			// Resize vectors at rotamer r1
			pairEnergy[posIndexI][r1].resize(sysSize);


			// Compute self and template energy
			double s = pairwiseEnergy->calculateSelfEnergy(_sys, i, r1);
			double t = pairwiseEnergy->calculateTemplateEnergy(_sys, i, r1);
			selfEnergy[posIndexI][r1]     = s;
			templateEnergy[posIndexI][r1] = t;

			
			// Add self and template linked position energy
			vector<Position *> &linkedPositionsP1 = p1.getLinkedPositions();
			for (uint linked = 0;linked < linkedPositionsP1.size();linked++){
				selfEnergy[posIndexI][r1]     += pairwiseEnergy->calculateSelfEnergy(_sys,linkedPositionsP1[linked]->getIndex(),r1);
				templateEnergy[posIndexI][r1] += pairwiseEnergy->calculateTemplateEnergy(_sys,linkedPositionsP1[linked]->getIndex(),r1);

				// Add pairwise energy between linked positions (i.e. Master to Slave "self" energy), part of self energy for super-rotamers
				selfEnergy[posIndexI][r1]     += pairwiseEnergy->calculatePairEnergy(_sys, i, r1, linkedPositionsP1[linked]->getIndex(), r1); 
				
			}


			// No pair interactions for non-variable positions
			if (p1.getTotalNumberOfRotamers() <= 1) continue;

			int posIndexJ = posIndexI;

			// For each position j
			for (uint j = i+1; j < _sys.positionSize();j++){
				Position &p2 = _sys.getPosition(j);


				// Skip linked positions that are not "master"
				if (p2.getLinkedPositionType() == Position::SLAVE) { continue; }

				// Resize vector with new posIndexJ and number of rotamers.
				posIndexJ++;

				// Only compute pair energies between variable positions..
				if (p2.getTotalNumberOfRotamers() <= 1) continue;
				
				

				pairEnergy[posIndexI][r1][posIndexJ].resize(p2.getTotalNumberOfRotamers());

				// Get list of positions linked to p2.
				vector<Position *> &linkedPositionsP2 = p2.getLinkedPositions();


				// For each rotamer r2
				for (uint r2 = 0; r2 < p2.getTotalNumberOfRotamers();r2++){
			
					// Compute pair energy (when linked it is considered MasterP1-to-MasterP2)
					double e = pairwiseEnergy->calculatePairEnergy(_sys, i, r1, j, r2); 
					pairEnergy[posIndexI][r1][posIndexJ][r2] = e;


					// For linked positions...
					for (uint linked1 = 0;linked1 < linkedPositionsP1.size();linked1++){

						// Add SlavesP1-to-MasterP2 pairEnergy
						e = pairwiseEnergy->calculatePairEnergy(_sys, linkedPositionsP1[linked1]->getIndex(), r1, j, r2); 
						pairEnergy[posIndexI][r1][posIndexJ][r2] += e;
						

						for (uint linked2 = 0;linked2 < linkedPositionsP2.size();linked2++){

							// Add MasterP1-to-SlaveP2 pairEnergy
							if (linked1 == 0){
								e = pairwiseEnergy->calculatePairEnergy(_sys, i, r1, linkedPositionsP2[linked2]->getIndex(), r2); 
								pairEnergy[posIndexI][r1][posIndexJ][r2] += e;
							}			

							// Add SlaveP1-to-SlaveP2 pairEnergy
							e = pairwiseEnergy->calculatePairEnergy(_sys, linkedPositionsP1[linked1]->getIndex(), r1, linkedPositionsP2[linked2]->getIndex(), r2); 
							pairEnergy[posIndexI][r1][posIndexJ][r2] += e;
								
						} // FOR each linked2 
					} // FOR each linked1						
				} // FOR each R2
			} // FOR each p2
		} // FOR each r1
	} // FOR EACH p1
	//}

	for (i = 0; i < selfEnergy.size();i++){
		for (uint j = 0; j < selfEnergy[i].size();j++){
			energy += templateEnergy[i][j] + selfEnergy[i][j];  

			for (uint ii = 0; ii < pairEnergy[i][j].size();ii++){
				for (uint jj = 0; jj < pairEnergy[i][j][ii].size();jj++){
					energy += pairEnergy[i][j][ii][jj];
				}
			}
		}
		

	}



	return energy;
	
}
double PairwiseEnergyCalculator::getStateEnergy(System &_sys, vector<unsigned int> &_stateVector){


	double energy = 0.0;
	map<int,int> rotamersOfVariablePositions;
	int rotamerIndex = 0;
	for (uint i = 0; i < _sys.positionSize();i++){


		if (_sys.getPosition(i).getTotalNumberOfRotamers() > 1){
			cout << "Position : "<<_sys.getPosition(i).getChainId()<< " "<<_sys.getPosition(i).getResidueName()<<" "<<_sys.getPosition(i).getResidueNumber()<<" is position "<< i<< " and rotamer "<<_stateVector[rotamerIndex]<<endl;
			rotamersOfVariablePositions[i] = _stateVector[rotamerIndex];


			energy += selfEnergy[i][_stateVector[rotamerIndex]];
			energy += templateEnergy[i][_stateVector[rotamerIndex]];

			//fprintf(stdout, "\tSelf[%4d] = %8.3f\n",i,selfEnergy[i][_stateVector[rotamerIndex]]);
			//fprintf(stdout, "\tTemplate[%4d] = %8.3f\n",i,templateEnergy[i][_stateVector[rotamerIndex]]);

			// Increment the rotamer index 
			rotamerIndex++;
		} else {

			energy += selfEnergy[i][0];
			energy += templateEnergy[i][0];

			//fprintf(stdout, "\tSelf[%4d] = %8.3f\n",i,selfEnergy[i][0]);
			//fprintf(stdout, "\tTemplate[%4d] = %8.3f\n",i,templateEnergy[i][0]);
		}

	}

	map<int,int>::iterator itI;
	map<int,int>::iterator itJ;

	for (itI = rotamersOfVariablePositions.begin();itI != rotamersOfVariablePositions.end();itI++){
		itJ = itI;
		itJ++;
		for (;itJ != rotamersOfVariablePositions.end();itJ++){
			if (itI->first <= itJ->first){
				energy += pairEnergy[itI->first][itI->second][itJ->first][itJ->second];
				//fprintf(stdout, "\tPair[%4d][%4d] = %8.3f\n",itI->first,itJ->first,pairEnergy[itI->first][itI->second][itJ->first][itJ->second]);
			} else {
				energy += pairEnergy[itJ->first][itJ->second][itI->first][itI->second];
				//fprintf(stdout, "\tPair[%4d][%4d] = %8.3f\n",itJ->first,itI->first,pairEnergy[itJ->first][itJ->second][itI->first][itI->second]);
			}
		}

	}


	return energy;


}
void PairwiseEnergyCalculator::printSummary(){

	
	map<string,double> energies = pairwiseEnergy->getAllComputedEnergiesByType();
	map<string,double>::iterator it;

	fprintf(stdout, "================  ======================  ===============\n");
	fprintf(stdout, "%-20s%-22s\n", "Interaction Type", "Energy");
	fprintf(stdout, "================  ======================  ===============\n");
	double energyTotal  = 0.0;
	double storedEnergyTotal = 0.0;
	for (it = energies.begin();it != energies.end();it++){
		if (it->first == "TOTAL"){
			storedEnergyTotal = it->second;
		}else {
			fprintf(stdout,"%-20s%20f\n",it->first.c_str(),it->second);
			energyTotal += it->second;
		}
	}

	fprintf(stdout, "================  ======================  ===============\n");
	fprintf(stdout, "%-20s%20f%20f\n","TOTAL",energyTotal,storedEnergyTotal);
	fprintf(stdout, "================  ======================  ===============\n\n");
}

void PairwiseEnergyCalculator::printPairwiseTable(){
	
	for (uint i = 0; i < pairEnergy.size();i++){
		for (uint r1 = 0; r1 < pairEnergy[i].size();r1++){
			for (uint ii = 0; ii < pairEnergy[i][r1].size();ii++){
				for (uint r2 = 0; r2 < pairEnergy[i][r1][ii].size();r2++){

					fprintf(stdout,"%4d %4d %4d %4d %8.3f\n",i,r1,ii,r2,pairEnergy[i][r1][ii][r2]);
				}
			}
		}
		
	}
}
