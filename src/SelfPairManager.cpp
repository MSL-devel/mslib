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

#include "SelfPairManager.h"

using namespace MSL;
using namespace std;


SelfPairManager::SelfPairManager() {
	setup();
}

SelfPairManager::SelfPairManager(System * _pSystem) {
	setup();
	setSystem(_pSystem);
}

SelfPairManager::~SelfPairManager() {
}

void SelfPairManager::setup() {
	pSys = NULL;
	pESet = NULL;
}

void SelfPairManager::copy(const SelfPairManager & _sysBuild) {
}

void SelfPairManager::deletePointers() {
}

void SelfPairManager::findVariablePositions() {

	/*************************************************************************
	 *
	 *    vector<vector<vector<vector<map<string, vector<Interaction*> > > > > > subdividedInteractions:
	 *
	 *    4D vector containing a map of interactions by string (i.e. "CHARMM_VDW", "CHARMM_ANGL" ...
	 *    just like in the EnergySet, by positions i and j and residue type (identity) ii and jj
	 *
	 *    subdividedInteractions[i][ii][j][jj]
	 *
	 *    The table is subdivided by i and j in 4 sections
	 *    		1) fixed energies: interactions of atoms all belonging to a fixed position
	 *    		2) self energies: interactions of atoms all belonging to a single variable position
	 *    		3) variable-fixed energies: interactions with some atoms belonging on fixed and some to 
	 *    		   a single variable position
	 *    		4) pair energies: interactions that belong to two variable positions (and possibly some
	 *    		   fixed position atoms as well)
	 *
	 *        _____    -- j --
	 *        |0,0|<------------------ [0,0] fixed energies
	 *     |  -----
	 *     |  |1,0| 1,1          \
	 *        |   |_____          \
	 *     i  |2,0||2,1| 2,2       \
	 *        |   ||   |----|       \---- diagonal cells ([1,1]...[4,4] self energies
	 *     |  |3,0||3,1  3,2| 3,3    \
	 *     |  |   ||        |----|    \
	 *        |4,0||4,1  4,2  4,3| 4,4 \
	 *        ---|-------|  
	 *         |     |        
	 *         |     |
	 *         |     inner cells ([2,1], [3,1], [3,2] ... [4,2], [4,3]) pair energies
	 *         |
	 *         first column ([1,0]...[4,0]), variable-fixed energies
	 *
	 *    Each cell is further divided by identity (all 2D even for the fixed and
	 *    self, the fixed uses only [0,0] and the self uses only [i,0]
	 *
	 *    So for example 
	 *    
	 *    subdividedInteractions[0][0][0][0] is a map<string, vector<Interaction*> containing
	 *       all interactions of the fixed
	 *
	 *    subdividedInteractions[3][5][0][0] is a map<string, vector<Interaction*> containing
	 *       all the variable/fixed interactions of the 4th variable position in its 6th identity 
	 *
	 *    subdividedInteractions[3][5][3][0] is a map<string, vector<Interaction*> containing
	 *       all the self interactions of the 4th variable position in its 6th identity 
	 *
	 *    subdividedInteractions[3][5][1][3] is a map<string, vector<Interaction*> containing
	 *       all the pair interactions of the 4th variable position in its 6th identity and
	 *       the 2nd variable position in its 4th identity 
	 *
	 *    subdividedInteractions[3][5][1][3]["CHARMM_VDW"] is a vector<Interaction*> of all the
	 *       vdw interactions of this pair of identities
	 *
	 *    subdividedInteractions[3][5][1][3]["CHARMM_VDW"][4] is a specific Interaction pointer
	 *       for the vdw interaction of two atoms in this pair
	 *
	 *
	 ***************************************************************************/

	/**************************************
	 *  TO DO: ADD ALL RESETS!!!!!
	 **************************************/
	variableCount.clear();
	subdividedInteractions.clear();
	variableIdentities.clear();
	variablePositions.clear();
	variablePosIndex.clear();

	if (pSys != NULL) {
		vector<Position*> & positions = pSys->getPositions();

		// add the first box for the fixed interactions
		subdividedInteractions.push_back(vector<vector<vector<map<string, vector<Interaction*> > > > >(1, vector<vector<map<string, vector<Interaction*> > > >(1, vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()))));
		variablePositions.push_back(NULL); // an NULL position for the fixed (just a trick)
		variableIdentities.push_back(vector<Residue*>(1, (Residue*)NULL)); // an NULL identity for the fixed (just a trick)

		unsigned int varCounter = 0;
		for (unsigned int i=0; i<positions.size(); i++) {
			unsigned int totalRots = positions[i]->getTotalNumberOfRotamers();
			if (totalRots > 1) {
				varCounter++;
				//variablePosMap[positions[i]] = true;
				variablePosIndex[positions[i]] = varCounter;
				variableCount.push_back(totalRots);
				// add a new entry for each identity of the position
				subdividedInteractions.push_back(vector<vector<vector<map<string, vector<Interaction*> > > > >(positions[i]->size(), vector<vector<map<string, vector<Interaction*> > > >()));
				variablePositions.push_back(positions[i]);
				variableIdentities.push_back(vector<Residue*>());

				// to each identity, add an entry for each previous variable position
				for (unsigned int ii=0; ii<subdividedInteractions.back().size(); ii++) {
					subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()));
					for (unsigned int j=1; j<subdividedInteractions.size()-1; j++) {
						subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(subdividedInteractions[j].size(), map<string, vector<Interaction*> >()));

					}
					subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()));

					// save a pointer to the ii-th identity of this variable position
					variableIdentities.back().push_back(&(positions[i]->getIdentity(ii)));
				}
			} else {
				//variablePosMap[positions[i]] = false;
				variablePosIndex[positions[i]] = 0;  // index 0 is associated with fixed positions
			}
		}
	}
}

void SelfPairManager::subdivideInteractions() {
	for (map<string, vector<Interaction*> >::iterator k=pEnergyTerms->begin(); k!=pEnergyTerms->end(); k++) {
		for (vector<Interaction*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			vector<Atom*> & atoms = (*l)->getAtoms();
			unsigned int variableCounter = 0;
			unsigned int positionOne = 0;
			unsigned int positionTwo = 0;
			unsigned int identityOne = 0;
			unsigned int identityTwo = 0;
			bool isFixed = false;
			bool isVariable = false;
			for (vector<Atom*>::iterator m=atoms.begin(); m!=atoms.end(); m++) {
				Position * pPos = (*m)->getParentPosition();
				if (variablePosIndex[pPos] == 0) {
					isFixed = true;
					if (variableCounter == 1) {
						positionTwo = 0;
						identityTwo = 0;
					}
				} else {
					isVariable = true;
					unsigned int identityIndex = (*m)->getIdentityIndex();
					if (variableCounter == 0) {
						variableCounter++;
						positionOne = variablePosIndex[pPos];
						identityOne = identityIndex;
						if (!isFixed) {
							positionTwo = positionOne;
							identityTwo = 0;
						}
					} else if (variableCounter == 1) {
						if (variablePosIndex[pPos] != positionOne) {
							variableCounter++;
							if (variablePosIndex[pPos] < positionOne) {
								positionTwo = variablePosIndex[pPos];
								identityTwo = identityIndex;
							} else {
								positionTwo = positionOne;
								positionOne = variablePosIndex[pPos];
								identityTwo = identityOne;
								identityOne = identityIndex;
							}
						} else {
							if (identityOne != identityIndex) {
								cerr << "ERROR 54907: mismatching identities in " << (*l)->toString() << " refers to more than two variable positions in void SelfPairManager::subdivideInteractions()" << endl;
								exit(54907);
							}
						}
					} else if (variableCounter == 2 && variablePosIndex[pPos] != positionOne && variablePosIndex[pPos] != positionTwo) {
						cerr << "ERROR 54912: Interaction " << (*l)->toString() << " refers to more than two variable positions in void SelfPairManager::subdivideInteractions()" << endl;
						exit(54912);
					}
				}
			}

			subdividedInteractions[positionOne][identityOne][positionTwo][identityTwo][k->first].push_back(*l);

		}
	}

}

void SelfPairManager::calculateEnergies() {

	fixE = 0.0;
	selfE.clear();
	pairE.clear();

	fixEbyTerm.clear();
	selfEbyTerm.clear();
	pairEbyTerm.clear();

	fixCount = 0;
	selfCount.clear();
	pairCount.clear();

	fixCountByTerm.clear();
	selfCountByTerm.clear();
	pairCountByTerm.clear();
	
	IdRotAbsIndex.clear();
	vector<unsigned int> rotCoor(3, 0);

	rotamerDescriptors.clear();
	rotamerPos_Id_Rot.clear();

	/*********************************************************
	 *  HERE!!!!
	 *    invert [i][j][ii][jj]
	 *
	 *    to [i][ii][j][jj]
	 *
	 *
	 *    DONE?
	 *
	 *********************************************************/
	for (unsigned int i=0; i<subdividedInteractions.size(); i++) {
		// LOOP LEVEL 1: for each position i

		if (i>0) {
			selfE.push_back(vector<double>());
			pairE.push_back(vector<vector<vector<double> > >());

			selfEbyTerm.push_back(vector<map<string, double> >());
			pairEbyTerm.push_back(vector<vector<vector<map<string, double> > > >());

			selfCount.push_back(vector<unsigned int>());
			pairCount.push_back(vector<vector<vector<unsigned int> > >());

			selfCountByTerm.push_back(vector<map<string, unsigned int> >());
			pairCountByTerm.push_back(vector<vector<vector<map<string, unsigned int> > > >());

			rotamerDescriptors.push_back(vector<string>());
			rotamerPos_Id_Rot.push_back(vector<vector<unsigned int> >());
		}
		
		unsigned int rotI = -1;
		unsigned int overallConfI = 0;
		for (unsigned int ii=0; ii<subdividedInteractions[i].size(); ii++) {
			// LOOP LEVEL 2: for each identity ii at position i
			

			// get the number of rotamer for the position and their particular identity
			unsigned int totalConfI = 1;
			if (i != 0) {
				// a self or a pair interaction
				totalConfI = variableIdentities[i][ii]->getNumberOfAltConformations();
				// set the i/ii-th residue the initial rotamer
				variableIdentities[i][ii]->setActiveConformation(0);

				string chain = variableIdentities[i][ii]->getChainId();
				string resName = variableIdentities[i][ii]->getResidueName();
				int resNum = variableIdentities[i][ii]->getResidueNumber();
				string iCode = variableIdentities[i][ii]->getResidueIcode();
				unsigned int posIndex = variablePosIndex[variablePositions[i]];
				overallConfI += totalConfI;

				for (unsigned int cI=0; cI<totalConfI; cI++) {
					char c [1000];
					sprintf(c, "Position %1s %4d%1s (%4u), identity %-4s (%2u), rotamer %3u (%3u)", chain.c_str(), resNum, iCode.c_str(), posIndex, resName.c_str(), ii, cI, overallConfI - totalConfI + cI);
					rotamerDescriptors[i-1].push_back((string)c);
					rotamerPos_Id_Rot[i-1].push_back(vector<unsigned int>());
					rotamerPos_Id_Rot[i-1].back().push_back(posIndex);
					rotamerPos_Id_Rot[i-1].back().push_back(ii);
					rotamerPos_Id_Rot[i-1].back().push_back(cI);
				}
			}

			for (unsigned int cI=0; cI<totalConfI; cI++) {
				//  LOOP LEVEL 3: for each rotamer of pos/identity i/ii 

				rotI++;
				if (i>0) {
					selfE[i-1].push_back(0.0);
					pairE[i-1].push_back(vector<vector<double> >());

					selfEbyTerm[i-1].push_back(map<string, double>());
					pairEbyTerm[i-1].push_back(vector<vector<map<string, double> > >());
					
					selfCount[i-1].push_back(0);
					pairCount[i-1].push_back(vector<vector<unsigned int> >());

					selfCountByTerm[i-1].push_back(map<string, unsigned int>());
					pairCountByTerm[i-1].push_back(vector<vector<map<string, unsigned int> > >());
					
					// update the look up table from identity/rot to absolute rot
					rotCoor[0] = i;
					rotCoor[1] = ii;
					rotCoor[2] = cI;
					IdRotAbsIndex[rotCoor] = rotI;
				}

				if (cI > 0) {
					// change the rotamer of i/ii
					variableIdentities[i][ii]->setActiveConformation(cI);
				}

				for (unsigned int j=0; j<subdividedInteractions[i][ii].size(); j++) {
					// LOOP LEVEL 4: for each position j

					unsigned int type = 0; // fixed
					if (i>0) {
						if (j == i || j == 0) {
							type = 1; // self
						} else {
							type = 2; //pair
						}
					}

					if (type == 2) {
						pairE[i-1][rotI].push_back(vector<double>());

						pairEbyTerm[i-1][rotI].push_back(vector<map<string, double> >());

						pairCount[i-1][rotI].push_back(vector<unsigned int>());

						pairCountByTerm[i-1][rotI].push_back(vector<map<string, unsigned int> >());
					}

					unsigned int rotJ = -1;
					for (unsigned int jj=0; jj<subdividedInteractions[i][ii][j].size(); jj++) {
						// LOOP LEVEL 5: for each identity jj

						// get the number of rotamer for the position and their particular identity
						unsigned int totalConfJ = 1;
						if (j != 0 && j != i) {
							// a pair interaction
							totalConfJ = variableIdentities[j][jj]->getNumberOfAltConformations();
							// set the j/jj-th residue the initial rotamer
							variableIdentities[j][jj]->setActiveConformation(0);
						}

						for (unsigned int cJ=0; cJ<totalConfJ; cJ++) {
							//  LOOP LEVEL 6: for each rotamer of pos/identity j/jj 


							rotJ++;
							if (type == 2) {
								pairE[i-1][rotI][j-1].push_back(0.0);

								pairEbyTerm[i-1][rotI][j-1].push_back(map<string, double>());

								pairCount[i-1][rotI][j-1].push_back(0);

								pairCountByTerm[i-1][rotI][j-1].push_back(map<string, unsigned int>());

							}
							if (cJ > 0) {
								// change the rotamer of j/jj
								variableIdentities[j][jj]->setActiveConformation(cJ);
							}

							// finally calculate the energies
							for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[i][ii][j][jj].begin(); k!= subdividedInteractions[i][ii][j][jj].end(); k++) {
								if (!pESet->isTermActive(k->first)) {
									// inactive term
									continue;
								}
								for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {

									double E = (*l)->getEnergy();
									if (type == 0) {
										fixE += E;
										fixEbyTerm[k->first] += E;

										fixCount++;
										fixCountByTerm[k->first]++;
									} else if (type == 1) {
										selfE[i-1][rotI] += E;
										selfEbyTerm[i-1][rotI][k->first] += E;

										selfCount[i-1][rotI]++;
										selfCountByTerm[i-1][rotI][k->first]++;
									} else {
										pairE[i-1][rotI][j-1][rotJ] += E;
										pairEbyTerm[i-1][rotI][j-1][rotJ][k->first] += E;

										pairCount[i-1][rotI][j-1][rotJ]++;
										pairCountByTerm[i-1][rotI][j-1][rotJ][k->first]++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	/*
	cout << "UUU Fixed Energy: " << fixE << " (" << fixCount << ")" << endl;
	for (unsigned int i=0; i<pairE.size(); i++) {
		for (unsigned int j=0; j<pairE[i].size(); j++) {
			cout << "UUU Self Energy " << i << "/" << j << ": " << selfE[i][j] << " (" << selfCount[i][j] << ")" << endl;
			for (unsigned int k=0; k<pairE[i][j].size(); k++) {
				for (unsigned int l=0; l<pairE[i][j][k].size(); l++) {
					cout << "UUU Pair Energy " << i << "/" << j << " - " << k << "/" << l << ": " << pairE[i][j][k][l] << " (" << pairCount[i][j][k][l] << ")" << endl;
				}
			}
		}
	}
	*/


}

double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term) {
	if (_overallRotamerStates.size() != pairE.size()) {
		cerr << "ERROR 54917: incorrect number of positions in input (" << _overallRotamerStates.size() << " != " << pairE.size() << " in double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
		exit(54917);
	}

	
	double out = 0.0;
	if (_term == "") {
		out += fixE;
		for (unsigned int i=0; i<pairE.size(); i++) {
			if (_overallRotamerStates[i] >= pairE[i].size()) {
				cerr << "ERROR 54922: incorrect number of rotamer in variable position " << i << " in input (" << _overallRotamerStates[i] << " >= " << pairE[i].size() << " in double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
				exit(54922);
			}
			out += selfE[i][_overallRotamerStates[i]];
			for (unsigned int j=0; j<i; j++) {
				out += pairE[i][_overallRotamerStates[i]][j][_overallRotamerStates[j]];
			}
		}
	} else {
		if (fixEbyTerm.find(_term) == fixEbyTerm.end()) {
			cerr << "WARNING 54922: unrecognized term " << _term << " in double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
			return 0.0;
		}
		out += fixEbyTerm[_term];
		for (unsigned int i=0; i<pairEbyTerm.size(); i++) {
			if (_overallRotamerStates[i] >= pairEbyTerm[i].size()) {
				cerr << "ERROR 54927: incorrect number of rotamer in variable position " << i << " in input (" << _overallRotamerStates[i] << " >= " << pairEbyTerm[i].size() << " in double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
				exit(54927);
			}
			out += selfEbyTerm[i][_overallRotamerStates[i]][_term];
			for (unsigned int j=0; j<i; j++) {
				out += pairEbyTerm[i][_overallRotamerStates[i]][j][_overallRotamerStates[j]][_term];
			}
		}
	}
	return out;
}

unsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term) {
	if (_overallRotamerStates.size() != pairCount.size()) {
		cerr << "ERROR 54932: incorrect number of positions in input (" << _overallRotamerStates.size() << " != " << pairCount.size() << " in unsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
		exit(54932);
	}

	
	unsigned int out = 0;
	if (_term == "") {
		out += fixCount;
		for (unsigned int i=0; i<pairCount.size(); i++) {
			if (_overallRotamerStates[i] >= pairCount[i].size()) {
				cerr << "ERROR 54937: incorrect number of rotamer in variable position " << i << " in input (" << _overallRotamerStates[i] << " >= " << pairCount[i].size() << " in unsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
				exit(54937);
			}
			out += selfCount[i][_overallRotamerStates[i]];
			for (unsigned int j=0; j<i; j++) {
				out += pairCount[i][_overallRotamerStates[i]][j][_overallRotamerStates[j]];
			}
		}
	} else {
		if (fixCountByTerm.find(_term) == fixCountByTerm.end()) {
			cerr << "WARNING 54922: unrecognized term " << _term << " inunsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
			return 0.0;
		}
		out += fixCountByTerm[_term];
		for (unsigned int i=0; i<pairCountByTerm.size(); i++) {
			if (_overallRotamerStates[i] >= pairCountByTerm[i].size()) {
				cerr << "ERROR 54927: incorrect number of rotamer in variable position " << i << " in input (" << _overallRotamerStates[i] << " >= " << pairCountByTerm[i].size() << " in unsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
				exit(54927);
			}
			out += selfCountByTerm[i][_overallRotamerStates[i]][_term];
			for (unsigned int j=0; j<i; j++) {
				out += pairCountByTerm[i][_overallRotamerStates[i]][j][_overallRotamerStates[j]][_term];
			}
		}
	}
	return out;
}

//double SelfPairManager::getStateEnergy(vector<unsigned int> _residueStates, vector<unsigned int> _rotamerStates) {
//}
//double SelfPairManager::getStateEnergy(vector<string> _residueNames, vector<unsigned int> _rotamerStates) {
//}

string SelfPairManager::getSummary(vector<unsigned int> _overallRotamerStates) {
	ostringstream os;	
	os << setiosflags(ios::left);
	os << "================  ======================  ===============" << endl;
	os << setw(20) <<"Interaction Type"<< setw(22) <<"Energy" << setw(15) << "Number of Terms" << endl;
	os << "================  ======================  ===============" << endl;
	for (map<string, double>::const_iterator l = fixEbyTerm.begin(); l!=fixEbyTerm.end(); l++) {
		double E = getStateEnergy(_overallRotamerStates, l->first);
		if (E<1E+14 && E>-1E+14) {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getStateInteractionCount(_overallRotamerStates, l->first) << endl;
		} else {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getStateInteractionCount(_overallRotamerStates, l->first) << endl;
		}
	}
	os << "================  ======================  ===============" << endl;
	double E = getStateEnergy(_overallRotamerStates);
	if (E<1E+14 && E>-1E+14) {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getStateInteractionCount(_overallRotamerStates) << endl;
	} else {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getStateInteractionCount(_overallRotamerStates) << endl;
	}
	os << "================  ======================  ===============" << endl;
	return (os.str());

}


