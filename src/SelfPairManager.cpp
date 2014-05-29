/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

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

#include "MslOut.h"
static MslOut MSLOUT("SelfPairManager");

SelfPairManager::SelfPairManager() {
	setup();
}

SelfPairManager::SelfPairManager(System * _pSystem) {
	setup();
	setSystem(_pSystem);
}

SelfPairManager::~SelfPairManager() {
	deletePointers();
}

void SelfPairManager::setup() {
	pSys = NULL;
	pESet = NULL;

	deleteRng = true;
	pRng = new RandomNumberGenerator;

	runDEE = true;
	runSCMFBiasedMC = true;
	runUnbiasedMC = true;
	runSCMF = true;
	runEnum = true;

	verbose = true;

	// DEE Options
	DEEenergyOffset = 0.0;
	DEEdoSimpleGoldsteinSingle = true;
	DEEdoSimpleGoldsteinPair = false;
	runDEE = true;

	// Enumeration Options
	enumerationLimit = 50000;

	// SCMF Options
	maxSavedResults = 100;
	SCMFtemperature = 300;
	SCMFcycles = 100;

	// save energies by term
	saveEbyTerm = false;
	saveInteractionCount = false;

	onTheFly = false; // precompute pair energies by default

	// MCO Options
	mcStartT = 1000.0;
	mcEndT = 0.5;
	mcCycles = 50000;
	mcShape  = MonteCarloManager::EXPONENTIAL;
	mcMaxReject = 2000;
	mcDeltaSteps = 100;
	mcMinDeltaE = 0.01;
}

void SelfPairManager::copy(const SelfPairManager & _sysBuild) {
}

void SelfPairManager::deletePointers() {
	if (deleteRng == true) {
		delete pRng;
	}
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
	 *       all the variable/fixed interactions of the 3th variable position in its 6th identity 
	 *
	 *    subdividedInteractions[3][5][3][0] is a map<string, vector<Interaction*> containing
	 *       all the self interactions of the 3th variable position in its 6th identity 
	 *
	 *    subdividedInteractions[3][5][1][3] is a map<string, vector<Interaction*> containing
	 *       all the pair interactions of the 3th variable position in its 6th identity and
	 *       the 1nd variable position in its 4th identity 
	 *
	 *    subdividedInteractions[3][5][1][3]["CHARMM_VDW"] is a vector<Interaction*> of all the
	 *       vdw interactions of this pair of identities
	 *
	 *    subdividedInteractions[3][5][1][3]["CHARMM_VDW"][4] is a specific Interaction pointer
	 *       for the vdw interaction of two atoms in this pair
	 *
	 *
	 ***************************************************************************/
	/*
	bool autoFindPositions = false;
	vector<bool> variablePositionFound_flag(_variablePositions.size(), false); // keep track of what position was found
	if (_variablePositions.size() == 0) {
		 autoFindPositions = true;
	}
	*/

	/**************************************
	 *  TO DO: ADD ALL RESETS!!!!! Done???
	 **************************************/
	variableCount.clear();
	subdividedInteractions.clear();
	variableIdentities.clear();
	slaveIdentities.clear();
	variablePositions.clear();
//	slavePositions.clear();
	variablePosIndex.clear();
	//slavePosIndex.clear();

	if (pSys != NULL) {
		vector<Position*> & positions = pSys->getPositions();

		// add the first box for the fixed interactions
		subdividedInteractions.push_back(vector<vector<vector<map<string, vector<Interaction*> > > > >(1, vector<vector<map<string, vector<Interaction*> > > >(1, vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()))));
		variablePositions.push_back(NULL); // an NULL position for the fixed (just a trick)
		variableIdentities.push_back(vector<Residue*>(1, (Residue*)NULL)); // an NULL identity for the fixed (just a trick)
		slaveIdentities.push_back(vector<vector<Residue*> >(1,vector<Residue*>(0, (Residue*)NULL))); // an NULL identity for the fixed (just a trick)

		unsigned int varCounter = 0;
		for (unsigned int i=0; i<positions.size(); i++) {
			unsigned int totalRots = positions[i]->getTotalNumberOfRotamers();

			/*
			// determine if this is a variable position, either because it was given, or if autofind is on, because
			// it has more than one rotamer
			bool isVariable = false;
			if (!autoFindPositions) {
				for (unsigned int j=0; j<_variablePositions.size(); j++) {
					if (MslTools::comparePositionIds(_variablePositions[j], positions[i]->getPositionId())) {
						isVariable = true;
						variablePositionFound_flag[j] = true;
						break;
					}
				}
			} else {
				if (totalRots > 1) {
					isVariable = true;
				}
			}
			*/
				
			if (pSys->isPositionVariable(i)) {
				// if it is UNLINKED or MASTER
				if (positions[i]->getLinkedPositionType() != Position::SLAVE) {
					varCounter++;
					//variablePosMap[positions[i]] = true;
					//cout << "VariablePosIndex["<<positions[i]<<"]: = "<<varCounter<<endl;
					variablePosIndex[positions[i]] = varCounter;
					variableCount.push_back(totalRots);
					// add a new entry for each identity of the position
					subdividedInteractions.push_back(vector<vector<vector<map<string, vector<Interaction*> > > > >(positions[i]->identitySize(), vector<vector<map<string, vector<Interaction*> > > >()));
					variablePositions.push_back(positions[i]);
				//	slavePositions.push_back(vector<Position*>());
					variableIdentities.push_back(vector<Residue*>());
					slaveIdentities.push_back(vector<vector<Residue*> >());

					vector<Position*> linkedPos = positions[i]->getLinkedPositions();
					unsigned int numberOfSlaves = linkedPos.size();

					// to each identity, add an entry for each previous variable position
					for (unsigned int ii=0; ii<subdividedInteractions.back().size(); ii++) {
						subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()));
						for (unsigned int j=1; j<subdividedInteractions.size()-1; j++) {
							subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(subdividedInteractions[j].size(), map<string, vector<Interaction*> >()));

						}
						subdividedInteractions.back()[ii].push_back(vector<map<string, vector<Interaction*> > >(1, map<string, vector<Interaction*> >()));

						// save a pointer to the ii-th identity of this variable position
						variableIdentities.back().push_back(&(positions[i]->getIdentity(ii)));
						slaveIdentities.back().push_back(vector<Residue*>(numberOfSlaves, (Residue*)NULL));
					}
					if (positions[i]->getLinkedPositionType() == Position::MASTER) {
						// if MASTER add the linked SLAVES
						//vector<Position*> linkedPos = pSys->getPosition(i).getLinkedPositions();

						for (unsigned int j=0; j<linkedPos.size(); j++) {

							bool matching = true;
							string msg;
							if (positions[i]->identitySize() == linkedPos[j]->identitySize()) {
								for (unsigned int id=0; id<positions[i]->identitySize(); id++) {
									if (positions[i]->getIdentity(id).getResidueName() == linkedPos[j]->getIdentity(id).getResidueName()) {
										if (positions[i]->getTotalNumberOfRotamers(id) != linkedPos[j]->getTotalNumberOfRotamers(id)) {
											// different number of rotamers for the id-th identity
											matching = false;
											msg = "Different number of rotamers for identity " + positions[i]->getIdentity(id).getResidueName() + " at linked positions " + positions[i]->getPositionId() + "-" + linkedPos[j]->getPositionId() + " (" + MslTools::intToString(positions[i]->getTotalNumberOfRotamers(id)) + " != " + MslTools::intToString(linkedPos[j]->getTotalNumberOfRotamers(id)) + ")";
											break;
										}
									} else {
										// different residue name for the id-th identity
										matching = false;
										msg = "Different residue names for the " + MslTools::intToString(j) + "-th identity at linked positions " + positions[i]->getPositionId() + "-" + linkedPos[j]->getPositionId() + " (" + positions[i]->getIdentity(id).getResidueName() + " != " + linkedPos[j]->getIdentity(id).getResidueName() + ")";
										break;
									}
								}
							} else {
								// different number of identities
								matching = false;
								msg = "Different number of identities at linked positions " + positions[i]->getPositionId() + "-" + linkedPos[j]->getPositionId() + " (" + MslTools::intToString(positions[i]->identitySize()) + " != " + MslTools::intToString(linkedPos[j]->identitySize()) + ")";
							}

							if (!matching) {
								cerr << "ERROR 48923: linked positions do not match: " << msg << " in void SelfPairManager::findVariablePositions()" << endl;
								exit(48923);
							}

							variablePosIndex[linkedPos[j]] = varCounter;
							//slavePosIndex[linkedPos[j]] = varCounter;
							// add a new entry for each identity of the position
					//		slavePositions[i].push_back(linkedPos[j]);
							
							/************************************************************
							 *     TODO:
							 *
							 *     Here we need to check that the linked positions have
							 *     the same number of identities and each identity the
							 *     same number of rotamers
							 ************************************************************/
							for (unsigned int ii=0; ii<subdividedInteractions.back().size(); ii++) {
								// save a pointer to the ii-th identity of this slave position
								slaveIdentities.back()[ii][j] = &(linkedPos[j]->getIdentity(ii));
							}
						}

					}
				}
			} else {
				//variablePosMap[positions[i]] = false;
				variablePosIndex[positions[i]] = 0;  // index 0 is associated with fixed positions
			}
		}
	}
	/*
	for (unsigned int i=0; i<subdividedInteractions.size(); i++) {
		for (unsigned int ii=0; ii<subdividedInteractions[i].size(); ii++) {
			for (unsigned int j=0; j<subdividedInteractions[i][ii].size(); j++) {
				for (unsigned int jj=0; jj<subdividedInteractions[i][ii][j].size(); jj++) {
					cout << i << "," << ii << "," << j << "," << jj << endl;
				}
			}
		}
	}
	*/

	/*
	for (unsigned int i=0; i<_variablePositions.size(); i++) {
		if (!variablePositionFound_flag[i]) {
			cerr << "WARNING 81145: variable position " << _variablePositions[i] << " not found in System at void SelfPairManager::findVariablePositions(vector<string> _variablePositions)" << endl;
		}
	}
	*/
}

void SelfPairManager::subdivideInteractions() {
	for (map<string, vector<Interaction*> >::iterator k=pEnergyTerms->begin(); k!=pEnergyTerms->end(); k++) {
		// for each term

		for (vector<Interaction*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			// for each interaction
			vector<Atom*> & atoms = (*l)->getAtomPointers();
			unsigned int variableCounter = 0;
			unsigned int positionOne = 0;
			unsigned int positionTwo = 0;
			unsigned int identityOne = 0;
			unsigned int identityTwo = 0;
			bool isFixed = false;
		//	bool isVariable = false;
			bool skipInteraction = false;
			for (vector<Atom*>::iterator m=atoms.begin(); m!=atoms.end(); m++) {
				Position * pPos = (*m)->getParentPosition();
				if((*m)->getHidden()) {
					skipInteraction = true;
					break;
				}
				if (variablePosIndex[pPos] == 0) {
					// a fixed position, the interction could be fixed or self
					isFixed = true;
					if (variableCounter == 1) {
						// it is a self, pos 2 index 0 0 for self interactions
						positionTwo = 0;
						identityTwo = 0;
					}
					if (pPos->getNumberOfIdentities() > 1) {
						// this is a fixed position even if it has multiple identitities,
						// it must have been set manually
						// skip all the interactions of the residues that are not the current
						if ((*m)->getParentResidue() != &(pPos->getCurrentIdentity())) {
							//cout << "UUU skipping " << *((*m)->getParentResidue()) << " for " << *pPos << endl;
							skipInteraction = true;
						}
					}
				} else {
					// a variable position, interaction could be self or pair
			//		isVariable = true;
					unsigned int identityIndex = (*m)->getIdentityIndex();
					if (variableCounter == 0) {
						// first variable atom found so far
						variableCounter++;
						positionOne = variablePosIndex[pPos];
						identityOne = identityIndex;
						if (!isFixed) {
							// if no fixed atom was found before we set it as self, with pos 2 index = pos 1 index and the indentity set to 0
							positionTwo = positionOne;
							identityTwo = 0;
						}
					} else if (variableCounter == 1) {
						// we found a variable atom before
						if (variablePosIndex[pPos] != positionOne) {
							// this is a different position, must be a pair interaction
							variableCounter++;
							// the next is to maintaing a lower diagonal table
							if (variablePosIndex[pPos] < positionOne) {
								positionTwo = variablePosIndex[pPos];
								identityTwo = identityIndex;
							} else {
								// invert the indeces
								positionTwo = positionOne;
								positionOne = variablePosIndex[pPos];
								identityTwo = identityOne;
								identityOne = identityIndex;
							}
						} else {
							// it is the same position, still a self
							if (identityOne != identityIndex) {
								/************************************************************************************************************************* 
								*  If we get here it means 
								*  1) the two positions are the same and the identities are different or
								*  2) the two positions are master and slave and the identities are different
								* in case 2) we can be sure that the two are master and slave because  variablePosIndex[pPos] == positionOne in this block
								* in both cases the right thing to do is to skip the interaction
								*************************************************************************************************************************/

								skipInteraction = true;
								
								/*
								if (pPos->getLinkedPositionType() != Position::SLAVE) {
									// different identities of linked positions, not to be calculated, skip it
									skipInteraction = true;
								} else {
									// a problem
									cerr << "ERROR 54907: mismatching identities in " << (*l)->toString() << " refers to more than two variable positions in void SelfPairManager::subdivideInteractions()" << endl;
									exit(54907);
								}
								*/
							}
						}
					} else if (variableCounter == 2 && variablePosIndex[pPos] != positionOne && variablePosIndex[pPos] != positionTwo) {
						cerr << "ERROR 54912: Interaction " << (*l)->toString() << " refers to more than two variable positions in void SelfPairManager::subdivideInteractions()" << endl;
						exit(54912);
					}
				}
			}
			if (!skipInteraction) {
				subdividedInteractions[positionOne][identityOne][positionTwo][identityTwo][k->first].push_back(*l);
			}

		}
	}

}

void SelfPairManager::saveMin(double _boundE, vector<unsigned int> _stateVec, int _maxSaved) {
	// case the list is empty

	if (minBound.size() == 0) {
		minBound.push_back(_boundE);
		minStates.push_back(_stateVec);
		return;
	} else if (minBound.size() == _maxSaved) {
		if (_boundE >= minBound.back()) {
			return;
		}
	}
	
	// make sure that we have not yet saved the state
	vector<vector<unsigned int> >::iterator v;

	// ... else try to fit it in the middle
	v = minStates.begin();
	vector<double>::iterator b;
	int counter = 0;
	for (b = minBound.begin(); b != minBound.end(); b++, v++, counter++) {
		double Enew = _boundE;
		double Eold = *b;
		if (Enew <= Eold) {
			if (Enew == Eold) {
				// make sure that we have not yet saved the state
				bool identical = true;
				for (unsigned int i=0; i<_stateVec.size(); i++) {
					if (_stateVec[i] != (*v)[i]) {
						identical = false;
						break;
					}
				}

				if (identical) {
					return;
				}
			}
			minBound.insert(b, _boundE);
			minStates.insert(v, _stateVec);
			// trim the list if oversized
			if (minBound.size() > _maxSaved) {
				b = minBound.end()-1;
				v = minStates.end()-1;
				minBound.erase(b);
				minStates.erase(v);
			}
			return;
		}
	}

	// ... or stick it at the end if the list is short
	minBound.push_back(_boundE);
	minStates.push_back(_stateVec);

}

void SelfPairManager::calculateEnergies() {
	calculateFixedEnergies();
	calculateSelfEnergies();
	calculatePairEnergies();
}
void SelfPairManager::recalculateNonSavedEnergies(vector<vector<vector<vector<bool> > > > savedPairEnergies) {
	calculateFixedEnergies();
	calculateSelfEnergies();
	recalculateNonSavedPairEnergies(savedPairEnergies);
}
void SelfPairManager::calculateFixedEnergies() {
	fixE = 0.0;
	fixEbyTerm.clear();
	fixCount = 0;
	fixCountByTerm.clear();
	// fixed energy

	for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[0][0][0][0].begin(); k!= subdividedInteractions[0][0][0][0].end(); k++) {
		if (!pESet->isTermActive(k->first)) {
			continue;
		}
		if(saveEbyTerm) {
			fixEbyTerm[k->first] = 0.0;
			fixCountByTerm[k->first] = k->second.size();
		}
		if(saveInteractionCount) {
			fixCount += k->second.size();
		}
		double E = 0.0;
		for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
			E += (*l)->getEnergy();
		}
		E *= weights[k->first];
		fixE += E;
		if(saveEbyTerm) {
			fixEbyTerm[k->first] += E; 
		}

	}
}
void SelfPairManager::calculateSelfEnergies() {
	selfE.clear();
	selfEbyTerm.clear();
	selfCount.clear();
	selfCountByTerm.clear();

	//IdRotAbsIndex.clear();
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
	//cout << "SUBDIV: "<<subdividedInteractions.size()<<endl;
	// interaction with the fixed atoms
	for (unsigned int i=1; i<subdividedInteractions.size(); i++) {
		// LOOP LEVEL 1: for each position i
		selfE.push_back(vector<double>());
		if(saveInteractionCount) {
			selfCount.push_back(vector<unsigned int>());
		}
		if(saveEbyTerm) {
			selfEbyTerm.push_back(vector<map<string, double> >());
			selfCountByTerm.push_back(vector<map<string, unsigned int> >());
		}

		rotamerDescriptors.push_back(vector<string>());
		rotamerPos_Id_Rot.push_back(vector<vector<unsigned int> >());

		unsigned overallConfI = 0;
		//cout << "\tSUBDIV: "<<subdividedInteractions[i].size()<<endl;
		for (unsigned int ii=0; ii<subdividedInteractions[i].size(); ii++) {
			// LOOP LEVEL 2: for each identity ii at position i

			unsigned int totalConfI = variableIdentities[i][ii]->getNumberOfRotamers();
			string chain = variableIdentities[i][ii]->getChainId();
			string resName = variableIdentities[i][ii]->getResidueName();
			int resNum = variableIdentities[i][ii]->getResidueNumber();
			string iCode = variableIdentities[i][ii]->getResidueIcode();
			unsigned int posIndex = variablePosIndex[variablePositions[i]];
			//cout << "Position: "<<resNum<<" "<<resName<<endl;
			//cout << "\t\tRotamers: "<<totalConfI<<endl;
			for (unsigned int cI=0; cI<totalConfI; cI++) {
				char c [1000];
				// TO DO: adde the linked to the string
				sprintf(c, "Position %1s %4d%1s (%4u), identity %-4s (%2u), rotamer %3u (%3u)", chain.c_str(), resNum, iCode.c_str(), posIndex, resName.c_str(), ii, cI, overallConfI + cI);
				rotamerDescriptors[i-1].push_back((string)c);
				rotamerPos_Id_Rot[i-1].push_back(vector<unsigned int>());
				rotamerPos_Id_Rot[i-1].back().push_back(posIndex);
				rotamerPos_Id_Rot[i-1].back().push_back(ii);
				rotamerPos_Id_Rot[i-1].back().push_back(cI);

			}

			for(unsigned int cI = 0; cI < totalConfI; cI++) {
				// LOOP LEVEL 3: for each conformation  compute selfE and if greater than threshold discard it
				variableIdentities[i][ii]->setActiveConformation(cI);
			//	cout << "UUU " << i << " " << slaveIdentities.size() << endl;
			//	cout << "  UUU " << ii << " " << slaveIdentities[i].size() << endl;
				for (unsigned int iii=0; iii<slaveIdentities[i][ii].size(); iii++) {
					slaveIdentities[i][ii][iii]->setActiveConformation(cI);
				}
				double energy = 0;
				unsigned count = 0;
				map<string,unsigned int> countByTerm;
				map<string,double> energyByTerm;
				for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[i][ii][0][0].begin(); k!= subdividedInteractions[i][ii][0][0].end(); k++) {

					// LOOP LEVEL 4 for each energy term
					if (!pESet->isTermActive(k->first)) {
						// inactive term
						continue;
					}
					if(saveEbyTerm) {
						countByTerm[k->first] = k->second.size(); 
					}
					count += k->second.size();
					double E = 0.0;
					for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
						E += (*l)->getEnergy();
						//double e = (*l)->getEnergy();
						//energy += e;
						//if(saveEbyTerm) {
						//	energyByTerm[k->first] += e; 
						//}
					}
					E *= weights[k->first];
					energy += E;
					if(saveEbyTerm) {
						energyByTerm[k->first] += E; 
					}
				}

				// compute energies with self
				for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[i][ii][i][0].begin(); k!= subdividedInteractions[i][ii][i][0].end(); k++) {
					// LOOP LEVEL 4 for each energy term
					if (!pESet->isTermActive(k->first)) {
						// inactive term
						continue;
					}
					if(saveEbyTerm) {
						countByTerm[k->first] += k->second.size(); 
					}
					if(saveInteractionCount) {
						count += k->second.size();
					}
					double E = 0.0;
					for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
						E += (*l)->getEnergy();
						//double e = (*l)->getEnergy();
						//energy += e;
						//if(saveEbyTerm) {
						//	energyByTerm[k->first] += e; 
						//}
					}
					E *= weights[k->first];
					energy += E;
					if(saveEbyTerm) {
						energyByTerm[k->first] += E; 
					}
				}
				//cout << "\t\t\tENERGY["<<selfE.size()-1<<"]["<<selfE.back().size()+1<<"]: "<<energy<<endl;
				selfE.back().push_back(energy);
				if(saveInteractionCount) {
					selfCount.back().push_back(count);
				}
				if(saveEbyTerm) {
					selfEbyTerm.back().push_back(energyByTerm);
					selfCountByTerm.back().push_back(countByTerm);
				}
			}
			overallConfI += totalConfI;
		}
	}

}
void SelfPairManager::recalculateNonSavedPairEnergies(vector<vector<vector<vector<bool> > > > savedPairEnergies) {
	pairEFlag = savedPairEnergies;

	for (unsigned int i=1; i<subdividedInteractions.size(); i++) {
		// LOOP LEVEL 1: for each position i
		unsigned int rotI = -1;

		for (unsigned int ii=0; ii<subdividedInteractions[i].size(); ii++) {
			// LOOP LEVEL 2: for each identity ii at position i

			// get the number of rotamer for the position and their particular identity
			unsigned int totalConfI = 1;
			// a self or a pair interaction
			totalConfI = variableIdentities[i][ii]->getNumberOfRotamers();
			// set the i/ii-th residue the initial rotamer
			variableIdentities[i][ii]->setActiveConformation(0);
			for (unsigned int iii=0; iii<slaveIdentities[i][ii].size(); iii++) {
				slaveIdentities[i][ii][iii]->setActiveConformation(0);
			}

			string chain = variableIdentities[i][ii]->getChainId();
			string resName = variableIdentities[i][ii]->getResidueName();
			string iCode = variableIdentities[i][ii]->getResidueIcode();

			for (unsigned int cI=0; cI<totalConfI; cI++) {
				//  LOOP LEVEL 3: for each rotamer of pos/identity i/ii 

				rotI++;
				if (cI > 0) {
					// change the rotamer of i/ii
					variableIdentities[i][ii]->setActiveConformation(cI);
					for (unsigned int iii=0; iii<slaveIdentities[i][ii].size(); iii++) {
						slaveIdentities[i][ii][iii]->setActiveConformation(cI);
					}
				}

				for (unsigned int j=1; j<subdividedInteractions[i][ii].size(); j++) {
					// LOOP LEVEL 4: for each position j

					if (j == i ) {
						continue; // self
					}

					unsigned int rotJ = -1;
					for (unsigned int jj=0; jj<subdividedInteractions[i][ii][j].size(); jj++) {
						// LOOP LEVEL 5: for each identity jj

						// get the number of rotamer for the position and their particular identity
						unsigned int totalConfJ = variableIdentities[j][jj]->getNumberOfRotamers();
						// set the j/jj-th residue the initial rotamer
						variableIdentities[j][jj]->setActiveConformation(0);
						for (unsigned int jjj=0; jjj<slaveIdentities[j][jj].size(); jjj++) {
							slaveIdentities[j][jj][jjj]->setActiveConformation(0);
						}

						for (unsigned int cJ=0; cJ<totalConfJ; cJ++) {
							//  LOOP LEVEL 6: for each rotamer of pos/identity j/jj 

							rotJ++;

							if (cJ > 0) {
								// change the rotamer of j/jj
								variableIdentities[j][jj]->setActiveConformation(cJ);
								for (unsigned int jjj=0; jjj<slaveIdentities[j][jj].size(); jjj++) {
									slaveIdentities[j][jj][jjj]->setActiveConformation(cJ);
								}
							}

							if((saveEbyTerm || saveInteractionCount || !onTheFly) && !pairEFlag[i-1][rotI][j-1][rotJ]) {	
								// finally calculate the energies
								pairE[i-1][rotI][j-1][rotJ] = 0.0;

								for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[i][ii][j][jj].begin(); k!= subdividedInteractions[i][ii][j][jj].end(); k++) {
									if (!pESet->isTermActive(k->first)) {
										// inactive term
										continue;
									}
									if (saveEbyTerm) {
										pairEbyTerm[i-1][rotI][j-1][rotJ][k->first] = 0.0;
									}

									if(!onTheFly) {
										double E = 0.0;

										for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
											E += (*l)->getEnergy();
										}
										E *= weights[k->first];
										pairE[i-1][rotI][j-1][rotJ] += E;
										if(saveEbyTerm) {
											pairEbyTerm[i-1][rotI][j-1][rotJ][k->first] += E;
										}
									}
								}
								//cout << "N " << i-1 << " " << rotI << " " << j-1 << " " << rotJ << " : " << pairE[i-1][rotI][j-1][rotJ] << endl;
							}
							//else { 
							//	cout << "S " << i-1 << " " << rotI << " " << j-1 << " " << rotJ << " : " << pairE[i-1][rotI][j-1][rotJ] << endl;
							//}
						}
					}
				}
			}
		}
	}
	
	/*
	//cout << "UUU Fixed Energy: " << fixE << " (" << fixCount << ") " <<pairE.size()<<" "<<selfE.size()<<" "<<selfCount.size()<< endl;
	for (unsigned int i=0; i<pairE.size(); i++) {
		for (unsigned int j=0; j<pairE[i].size(); j++) {
		        //cout << "UUU Self Energy " << i << "/" << j << ": " << selfE[i][j] << endl;//" (" << selfCount[i][j] << ")" << endl;
			for (unsigned int k=0; k<pairE[i][j].size(); k++) {
				for (unsigned int l=0; l<pairE[i][j][k].size(); l++) {
				  cout << "UUU Pair Energy " << i << "/" << j << " - " << k << "/" << l << ": " << pairE[i][j][k][l] <<endl;  //<< " (" << pairCount[i][j][k][l] << ")" << endl;
				  energySum += pairE[i][j][k][l];
				}
			}
		}
	}
	*/

}

void SelfPairManager::calculatePairEnergies() {

	pairE.clear();
	pairEFlag.clear();
	pairEbyTerm.clear();
	pairCount.clear();
	pairCountByTerm.clear();

	for (unsigned int i=1; i<subdividedInteractions.size(); i++) {
		// LOOP LEVEL 1: for each position i

		// not a fixed energy, increment the tables
		pairE.push_back(vector<vector<vector<double> > >());
		if(onTheFly) {
			pairEFlag.push_back(vector<vector<vector<bool> > >());
		}

		if(saveInteractionCount) {
			pairCount.push_back(vector<vector<vector<unsigned int> > >());
		}

		if(saveEbyTerm) {
			pairEbyTerm.push_back(vector<vector<vector<map<string, double> > > >());
			pairCountByTerm.push_back(vector<vector<vector<map<string, unsigned int> > > >());
		}
		
		unsigned int rotI = -1;
		unsigned int overallConfI = 0;
		for (unsigned int ii=0; ii<subdividedInteractions[i].size(); ii++) {
			// LOOP LEVEL 2: for each identity ii at position i

			// get the number of rotamer for the position and their particular identity
			unsigned int totalConfI = 1;
			// a self or a pair interaction
			totalConfI = variableIdentities[i][ii]->getNumberOfRotamers();
			// set the i/ii-th residue the initial rotamer
			variableIdentities[i][ii]->setActiveConformation(0);
			for (unsigned int iii=0; iii<slaveIdentities[i][ii].size(); iii++) {
				slaveIdentities[i][ii][iii]->setActiveConformation(0);
			}

			string chain = variableIdentities[i][ii]->getChainId();
			string resName = variableIdentities[i][ii]->getResidueName();
			string iCode = variableIdentities[i][ii]->getResidueIcode();
			overallConfI += totalConfI;

			for (unsigned int cI=0; cI<totalConfI; cI++) {
				//  LOOP LEVEL 3: for each rotamer of pos/identity i/ii 

				rotI++;
				pairE[i-1].push_back(vector<vector<double> >());
				if(onTheFly) {
					pairEFlag[i-1].push_back(vector<vector<bool> >());
				}

				if(saveEbyTerm) {
					pairEbyTerm[i-1].push_back(vector<vector<map<string, double> > >());
					pairCountByTerm[i-1].push_back(vector<vector<map<string, unsigned int> > >());
				}
				
				if(saveInteractionCount) {
					pairCount[i-1].push_back(vector<vector<unsigned int> >());
				}
				
				
				if (cI > 0) {
					// change the rotamer of i/ii
					variableIdentities[i][ii]->setActiveConformation(cI);
					for (unsigned int iii=0; iii<slaveIdentities[i][ii].size(); iii++) {
						slaveIdentities[i][ii][iii]->setActiveConformation(cI);
					}
				}

				for (unsigned int j=1; j<subdividedInteractions[i][ii].size(); j++) {
					// LOOP LEVEL 4: for each position j

					if (j == i ) {
						continue; // self
					}

					pairE[i-1][rotI].push_back(vector<double>());
					if(onTheFly) {
						pairEFlag[i-1][rotI].push_back(vector<bool>());
					}
					if(saveInteractionCount) {
						pairCount[i-1][rotI].push_back(vector<unsigned int>());
					}

					if(saveEbyTerm) {
						pairEbyTerm[i-1][rotI].push_back(vector<map<string, double> >());
						pairCountByTerm[i-1][rotI].push_back(vector<map<string, unsigned int> >());
					}

					unsigned int rotJ = -1;
					for (unsigned int jj=0; jj<subdividedInteractions[i][ii][j].size(); jj++) {
						// LOOP LEVEL 5: for each identity jj

						// get the number of rotamer for the position and their particular identity
						unsigned int totalConfJ = variableIdentities[j][jj]->getNumberOfRotamers();
						// set the j/jj-th residue the initial rotamer
						variableIdentities[j][jj]->setActiveConformation(0);
						for (unsigned int jjj=0; jjj<slaveIdentities[j][jj].size(); jjj++) {
							slaveIdentities[j][jj][jjj]->setActiveConformation(0);
						}

						for (unsigned int cJ=0; cJ<totalConfJ; cJ++) {
							//  LOOP LEVEL 6: for each rotamer of pos/identity j/jj 

							rotJ++;

							pairE[i-1][rotI][j-1].push_back(0.0);
							if(onTheFly) {
								pairEFlag[i-1][rotI][j-1].push_back(false);
							}

							if(saveEbyTerm) {
								pairEbyTerm[i-1][rotI][j-1].push_back(map<string, double>());
								pairCountByTerm[i-1][rotI][j-1].push_back(map<string, unsigned int>());
							}
							if(saveInteractionCount) {
								pairCount[i-1][rotI][j-1].push_back(0);
							}

							if (cJ > 0) {
								// change the rotamer of j/jj
								variableIdentities[j][jj]->setActiveConformation(cJ);
								for (unsigned int jjj=0; jjj<slaveIdentities[j][jj].size(); jjj++) {
									slaveIdentities[j][jj][jjj]->setActiveConformation(cJ);
								}
							}

							if(saveEbyTerm || saveInteractionCount || !onTheFly) {	
								// finally calculate the energies
								for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[i][ii][j][jj].begin(); k!= subdividedInteractions[i][ii][j][jj].end(); k++) {
									if (!pESet->isTermActive(k->first)) {
										// inactive term
										continue;
									}
									if(saveInteractionCount) {
										pairCount[i-1][rotI][j-1][rotJ] += k->second.size();
									}
									if(saveEbyTerm) {
										pairEbyTerm[i-1][rotI][j-1][rotJ][k->first] = 0;
										pairCountByTerm[i-1][rotI][j-1][rotJ][k->first] = k->second.size();
									}

									if(!onTheFly) {
										double E = 0.0;
										for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
											E += (*l)->getEnergy();
										}
										E *= weights[k->first];
										pairE[i-1][rotI][j-1][rotJ] += E;
										if(saveEbyTerm) {
											pairEbyTerm[i-1][rotI][j-1][rotJ][k->first] += E;
										}
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
	cout << "UUU Fixed Energy: " << fixE << " (" << fixCount << ") " <<pairE.size()<<" "<<selfE.size()<<" "<<selfCount.size()<< endl;
	for (unsigned int i=0; i<pairE.size(); i++) {
		for (unsigned int j=0; j<pairE[i].size(); j++) {
		        cout << "UUU Self Energy " << i << "/" << j << ": " << selfE[i][j] << endl;//" (" << selfCount[i][j] << ")" << endl;
			for (unsigned int k=0; k<pairE[i][j].size(); k++) {
				for (unsigned int l=0; l<pairE[i][j][k].size(); l++) {
				  cout << "UUU Pair Energy " << i << "/" << j << " - " << k << "/" << l << ": " << pairE[i][j][k][l] <<endl;  //<< " (" << pairCount[i][j][k][l] << ")" << endl;
				}
			}
		}
	}
	*/


}

double SelfPairManager::getAliveCombinations() const {
	double total = 1;
	for (unsigned int i=0; i<aliveMask.size(); i++) {
		unsigned int count = 0;
		for (unsigned int j=0; j<aliveMask[i].size(); j++) {
			if (aliveMask[i][j]) {
				count++;
			}
		}
		total *= count;
	}
	return total;
}


double SelfPairManager::computeSelfE(unsigned pos, unsigned rot) {
	// when onTheFly is implemented for selfE this will need to be updated
	return selfE[pos][rot];

}
double SelfPairManager::computePairE(unsigned pos1, unsigned rot1, unsigned pos2, unsigned rot2, string _term) {
	if (_term != "" && weights.find(_term) == weights.end()) {
		// term does not exist
		return 0.0;
	}
	if(!onTheFly) {
	  return pairE[pos1][rot1][pos2][rot2]; 
	}
	if(!pairEFlag[pos1][rot1][pos2][rot2]) {
	  // cout << pos1 << "," << rot1 << "," << pos2 << "," <<  rot2 << " " << endl;
		// We need to compute the identity number for pos1 and pos2 and set the correct rotamer number
		unsigned id1 = rotamerPos_Id_Rot[pos1][rot1][1]; // corresponding identity

		unsigned id2 = rotamerPos_Id_Rot[pos2][rot2][1]; // corresponding identity

		variablePositions[pos1 + 1]->setActiveRotamer(rot1);
		variablePositions[pos2 + 1]->setActiveRotamer(rot2);

		pairE[pos1][rot1][pos2][rot2] = 0.0;
		for (map<string, vector<Interaction*> >::iterator k=subdividedInteractions[pos1+1][id1][pos2+1][id2].begin(); k!= subdividedInteractions[pos1 + 1][id1][pos2 + 1][id2].end(); k++) {
			if (!pESet->isTermActive(k->first)) {
				// inactive term
				continue;
			}
			double E = 0.0;
			for (vector<Interaction*>::iterator l=k->second.begin(); l!= k->second.end(); l++) {
				E += (*l)->getEnergy();
			}
			E *= weights[k->first];
			pairE[pos1][rot1][pos2][rot2] += E;
			if(saveEbyTerm) {
				pairEbyTerm[pos1][rot1][pos2][rot2][k->first] += E;
			}
		}
		// Assume pairEFlag array exists 
		pairEFlag[pos1][rot1][pos2][rot2] = true;
	      }


	
	if (_term == "") {
		return pairE[pos1][rot1][pos2][rot2]; 
	} else {
		return pairEbyTerm[pos1][rot1][pos2][rot2][_term]; 
	}
}

/*
double SelfPairManager::computePairEbyTerm(unsigned pos1, unsigned rot1, unsigned pos2, unsigned rot2, string _term) {
	if(!saveEbyTerm) {
		return 0.0;
	}
	if(!onTheFly) {
		return pairEbyTerm[pos1][rot1][pos2][rot2][_term]; 
	}
	// We need to compute the identity number for pos1 and pos2 and set the correct rotamer number
	unsigned id1 = rotamerPos_Id_Rot[pos1][rot1][1]; // corresponding identity
	unsigned rotamer1 = rotamerPos_Id_Rot[pos1][rot1][2]; // rotamer number

	unsigned id2 = rotamerPos_Id_Rot[pos2][rot2][1]; // corresponding identity
	unsigned rotamer2 = rotamerPos_Id_Rot[pos2][rot2][2]; // rotamer number

	variableIdentities[pos1 + 1][id1]->setActiveConformation(rotamer1);
	variableIdentities[pos2 + 1][id2]->setActiveConformation(rotamer2);

	pairEbyTerm[pos1][rot1][pos2][rot2][_term] = 0.0;
	for (vector<Interaction*>::iterator k=subdividedInteractions[pos1 + 1][id1][pos2 + 1][id2][_term].begin(); k!= subdividedInteractions[pos1 + 1][id1][pos2 + 1][id2][_term].end(); k++) {
		pairEbyTerm[pos1][rot1][pos2][rot2][_term] += (*k)->getEnergy();
	}
	return pairEbyTerm[pos1][rot1][pos2][rot2][_term];
}
*/

double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term) {
	if (_term != "" && weights.find(_term) == weights.end()) {
		// term does not exist
		return 0.0;
	}
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
				out += computePairE(i,_overallRotamerStates[i],j,_overallRotamerStates[j]);
			}
		}
	} else {
		if (fixEbyTerm.find(_term) != fixEbyTerm.end()) {
			out += fixEbyTerm[_term];
		}
		for (unsigned int i=0; i<pairEbyTerm.size(); i++) {
			if (_overallRotamerStates[i] >= pairEbyTerm[i].size()) {
				cerr << "ERROR 54927: incorrect number of rotamer in variable position " << i << " in input (" << _overallRotamerStates[i] << " >= " << pairEbyTerm[i].size() << " in double SelfPairManager::getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term)" << endl;
				exit(54927);
			}
			out += selfEbyTerm[i][_overallRotamerStates[i]][_term];
			for (unsigned int j=0; j<i; j++) {
				out += computePairE(i,_overallRotamerStates[i],j,_overallRotamerStates[j],_term);
			}
		}
	}
	return out;
}

unsigned int SelfPairManager::getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term) {
	if (_term != "" && weights.find(_term) == weights.end()) {
		// term does not exist
		return 0;
	}
	if(!saveInteractionCount) {
		return 0;
	}
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
		if (fixCountByTerm.find(_term) != fixCountByTerm.end()) {
			out += fixCountByTerm[_term];
		}
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

string SelfPairManager::getSummary(vector<unsigned int> _overallRotamerStates, unsigned int _precision) {
	ostringstream os;	
	os << setiosflags(ios::left);
	if (saveInteractionCount) {
		os << "================  ======================  ===============" << endl;
		os << setw(20) <<"Interaction Type"<< setw(22) <<"Energy" << setw(15) << "Number of Terms" << endl;
		os << "================  ======================  ===============" << endl;
	} else {
		os << "================  ======================" << endl;
		os << setw(20) <<"Interaction Type"<< setw(22) <<"Energy" << endl;
		os << "================  ======================" << endl;
	}
	
	if (saveEbyTerm) {
		map<string,vector<Interaction*> > * eTerms = pESet->getEnergyTerms();
		for (map<string, vector<Interaction*> >::const_iterator l = eTerms->begin(); l!=eTerms->end(); l++) {
			if(!pESet->isTermActive(l->first)) {
				continue;
			}
			double E = getStateEnergy(_overallRotamerStates, l->first);
			if (E<1E+14 && E>-1E+14) {
				//os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getStateInteractionCount(_overallRotamerStates, l->first) << endl;
				os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(_precision) << E;
			//	os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << E;
			} else {
				os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(_precision) << "********************";
			//	os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************";
			}
			if (saveInteractionCount) {
				os << setw(15) << getStateInteractionCount(_overallRotamerStates, l->first) << endl;
			} else {
				os << endl;
			}
		}
		if (saveInteractionCount) {
			os << "================  ======================  ===============" << endl;
		} else {
			os << "================  ======================" << endl;
		}
	}
	double E = getStateEnergy(_overallRotamerStates);
	if (E<1E+14 && E>-1E+14) {
		//os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getStateInteractionCount(_overallRotamerStates) << endl;
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(_precision) << E;
	//	os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << E;
	} else {
		//os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getStateInteractionCount(_overallRotamerStates) << endl;
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(_precision) << "********************";
	//	os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << "********************";
	}
	if (saveInteractionCount) {
		os << setw(15) << getStateInteractionCount(_overallRotamerStates) << endl;
		os << "================  ======================  ===============" << endl;
	} else {
		os << endl;
		os << "================  ======================" << endl;
	}
	return (os.str());

}

double SelfPairManager::getFixedEnergy() const {
	return fixE;
}

std::vector<std::vector<double> > & SelfPairManager::getSelfEnergy() {
	return selfE;	
}

std::vector<std::vector<std::vector<std::vector<double> > > > & SelfPairManager::getPairEnergy() {
	return pairE;	
}

void SelfPairManager::setRunDEE(bool _singles, bool _pairs) {
	runDEE = _singles || _pairs;
	DEEdoSimpleGoldsteinSingle = _singles;
	DEEdoSimpleGoldsteinPair = _pairs;
}
void SelfPairManager::setRunUnbiasedMC(bool _toogle) {
	runUnbiasedMC = _toogle;
}
void SelfPairManager::setRunSCMFBiasedMC(bool _toogle) {
	runSCMFBiasedMC = _toogle;
}
void SelfPairManager::setRunSCMF(bool _toogle) {
	runSCMF = _toogle;
}
void SelfPairManager::setRunEnum(bool _toogle) {
	runEnum = _toogle;
}

void SelfPairManager::setVerbose(bool _toogle) {
	verbose = _toogle;
}

vector<unsigned int> SelfPairManager::getSCMFstate() {
	return mostProbableSCMFstate;
}

vector<unsigned int> SelfPairManager::getBestUnbiasedMCState() {
	return bestUnbiasedMCstate;
}

vector<unsigned int> SelfPairManager::getBestSCMFBiasedMCState() {
	return bestSCMFBiasedMCstate;
}

vector<vector<unsigned int> > SelfPairManager::getDEEAliveRotamers() {
	return aliveRotamers;
}

vector<vector<bool> > SelfPairManager::getDEEAliveMask() {
	return aliveMask;
}

vector<double> SelfPairManager::getMinBound() {
	return minBound;
}
vector<vector<unsigned int> > SelfPairManager::getMinStates() {
	return minStates;
}


void SelfPairManager::runOptimizer() {
	

	minBound.clear();
	minStates.clear();

	aliveRotamers.clear();
	aliveMask.clear();

	/*****************************************
	 *  Create the masks for SCMF
	 *  
	 *  A mask reduces the rotamer space at each position so that there is a
	 *  residue only, vs the whole rotameric space at all other positions
	 *****************************************/

	for (unsigned int i=0; i < selfE.size(); i++) {
		unsigned int rots = selfE[i].size();
		aliveMask.push_back(vector<bool>(rots, true));
		aliveRotamers.push_back(vector<unsigned int>());
		for (unsigned int j=0; j<rots; j++) {
			aliveRotamers.back().push_back(j);
		}
	}

	if(onTheFly && (runDEE || runSCMF)) {
		onTheFly = false;
		cerr << "WARNING 12324: DEE and/or SCMF need to be run, so precomputing all energies " << endl; 
		calculateEnergies();
	}

	double finalCombinations= 0;
	if (runDEE) {
		finalCombinations = runDeadEndElimination();
	} else {
		finalCombinations = getAliveCombinations();
	}

	if(runEnum) {
		if (finalCombinations > enumerationLimit) {
			if(verbose) {
				cout << "The number of combinations " << finalCombinations << " exceeds the limit (" << enumerationLimit << ") provided by the user - not enumerating." << endl;
			}
		} else {
			runEnumeration();
			return;
		}
	}
	 
	if(runSCMF) {
		runSelfConsistentMeanField();
	}
	if(runUnbiasedMC) {
		runUnbiasedMonteCarlo();
	}
}



double SelfPairManager::runDeadEndElimination() {
	
	/******************************************************************************
	 *                        === DEAD END ELIMINATION ===
	 ******************************************************************************/
	time_t startDEEtime, endDEEtime;
	double DEETime;

	time (&startDEEtime);

	//bool singleSolution = false;

	vector<vector<double> >& oligomersSelf = getSelfEnergy();
	vector<vector<vector<vector<double> > > >& oligomersPair = getPairEnergy();

	DeadEndElimination DEE(oligomersSelf, oligomersPair);
	double finalCombinations = DEE.getTotalCombinations();

	if (verbose) {
		DEE.setVerbose(true, 2);
	}
	else {
		DEE.setVerbose(false, 0);
	}
	DEE.setEnergyOffset(DEEenergyOffset);
	if (verbose) {
		cout << "===================================" << endl;
		cout << "Run Dead End Elimination" << endl;
		cout << "Initial combinations: " << DEE.getTotalCombinations() << endl;
	}
	while (true) {
		if (DEEdoSimpleGoldsteinSingle) {
			if (!DEE.runSimpleGoldsteinSingles() ) {
				break;
			}
		}
		unsigned int combinations = int(DEE.getTotalCombinations());
		if ((runEnum && combinations < enumerationLimit) || combinations == 1) {
			break;	
		}
		if (DEEdoSimpleGoldsteinPair) {
			if (!DEE.runSimpleGoldsteinPairsOnce()) {
				break;
			}
		}
		if (verbose) cout << "Current combinations: " << DEE.getTotalCombinations() << endl;
	}
	finalCombinations = DEE.getTotalCombinations();
	if (finalCombinations < enumerationLimit) {
		if (verbose) cout << "Alive combinations = " << finalCombinations << ": enumerate" << endl;
		//singleSolution = false;
		aliveRotamers = DEE.getAliveStates();
		aliveMask = DEE.getMask();
	} else {
		if (verbose) cout << "Alive combinations = " << finalCombinations << ": run SCMF/MC" << endl;
		//singleSolution = false;
		aliveRotamers = DEE.getAliveStates();
		aliveMask = DEE.getMask();
	}

	time (&endDEEtime);
	DEETime = difftime (endDEEtime, startDEEtime);

	if (verbose) {
		cout << endl << "DEE Time: " << DEETime << " seconds" << endl;
		cout << "===================================" << endl;
	}
	return finalCombinations;

}


void SelfPairManager::runEnumeration() {
	/******************************************************************************
	 *                     === ENUMERATION ===
	 ******************************************************************************/
	if (verbose) {
		cout << "===================================" << endl;
		cout << "Enumerate the states and find the mins" << endl;
	}

	Enumerator aliveEnum(aliveRotamers);

	for (int i=0; i<aliveEnum.size(); i++) {
		if (verbose) {
			cout << "State " << i << ":" << endl;
			//printIdentitiesAndRotamer(aliveEnum[i], opt);
			for (int j=0; j < aliveEnum[i].size(); j++){
				cout << aliveEnum[i][j] << ",";
			}
			cout << endl;
		}
		//double oligomerEnergy = oligomer.getTotalEnergyFromPairwise(aliveEnum[i]);
		double oligomerEnergy = getStateEnergy(aliveEnum[i]);
		saveMin(oligomerEnergy, aliveEnum[i], maxSavedResults);
	}
	if (verbose) {
		cout << "===================================" << endl;
	}
	return;

}


void SelfPairManager::runSelfConsistentMeanField() {

	/******************************************************************************
	 *                     === SELF CONSISTENT MEAN FIELD ===
	 ******************************************************************************/
	time_t startSCMFtime, endSCMFtime;
	double SCMFTime;

	time (&startSCMFtime);

	double oligomerFixed = getFixedEnergy();
	vector<vector<double> > & oligomersSelf = getSelfEnergy();
	vector<vector<vector<vector<double> > > >& oligomersPair = getPairEnergy();

	SelfConsistentMeanField SCMF;
	SCMF.setRandomNumberGenerator(pRng);
	SCMF.setEnergyTables(oligomerFixed,oligomersSelf,oligomersPair);
	SCMF.setT(SCMFtemperature);
	SCMF.setVerbose(verbose);
	if (runDEE) {
		SCMF.setMask(aliveMask);
	}

	if (verbose) {
		cout << "===================================" << endl;
		cout << "Run Self Consistent Mean Field" << endl;
		cout << endl;
	}
	for (int i=0; i < SCMFcycles; i++) {
		SCMF.cycle();
		if (verbose) {
			cout << "Cycle " << SCMF.getNumberOfCycles() << " p variation " << SCMF.getPvariation() << endl;
		}
	}
	time (&endSCMFtime);
	SCMFTime = difftime (endSCMFtime, startSCMFtime);
	mostProbableSCMFstate = SCMF.getMostProbableState();
	if (verbose) {
		cout << endl;
		cout << "Final SCMF probability variation: " << SCMF.getPvariation() << endl;
		cout << "Most probable state: ";
		for (int j=0; j < mostProbableSCMFstate.size(); j++){
			cout << mostProbableSCMFstate[j] << ",";
		}
		cout << endl;
		cout << "Most Probable State Energy: " << SCMF.getStateEnergy(mostProbableSCMFstate) << endl;
		cout << "SCMF Time: " << SCMFTime << " seconds" << endl;
		cout << "===================================" << endl;
	}
	saveMin(SCMF.getStateEnergy(mostProbableSCMFstate),mostProbableSCMFstate,maxSavedResults);
	if(runSCMFBiasedMC) {
		bestSCMFBiasedMCstate  = SCMF.runMC(mcStartT, mcEndT, mcCycles, mcShape, mcMaxReject, mcDeltaSteps, mcMinDeltaE);
		if(verbose) {
			cout << "Best SCMF biased state: ";
			for (int j=0; j < bestSCMFBiasedMCstate.size(); j++){
				cout << bestSCMFBiasedMCstate[j] << ",";
			}
			cout << endl;
		}

		saveMin(SCMF.getStateEnergy(bestSCMFBiasedMCstate),bestSCMFBiasedMCstate,maxSavedResults);
	}
}

void SelfPairManager::runUnbiasedMonteCarlo() {

	/******************************************************************************
	 *                     === Unbiased MONTE CARLO OPTIMIZATION ===
	 ******************************************************************************/

	vector<vector<double> > & oligomersSelf = getSelfEnergy();
	vector<vector<vector<vector<double> > > >& oligomersPair = getPairEnergy();

		
//	an unbiased monte carlo method using the most probable SCMF state as the start
	MonteCarloOptimization MCO;
	if(onTheFly) {
		MCO.setSelfPairManager(this);
	} else {
		MCO.addEnergyTable(oligomersSelf, oligomersPair);
	}
	MCO.setRandomNumberGenerator(pRng);
	if(runSCMF) {
		// set the starting state appropriately
		//MCO.setInitializationState(MSL::MonteCarloOptimization::RANDOM);
		//MCO.setInitializationState(MSL::MonteCarloOptimization::LOWESTSELF);
		//MCO.setInitializationState(MSL::MonteCarloOptimization::QUICKSCAN);
		if(runSCMFBiasedMC) {
			MCO.setInitializationState(bestSCMFBiasedMCstate);
		} else {
			MCO.setInitializationState(mostProbableSCMFstate);
		}
	}

	//MCO.setVerbose(verbose);
	if(runDEE) {
		MCO.setInputRotamerMasks(aliveMask);
	}
	bestUnbiasedMCstate  = MCO.runMC(mcStartT, mcEndT, mcCycles, mcShape, mcMaxReject, mcDeltaSteps, mcMinDeltaE);
	double bestEnergy =  getStateEnergy(bestUnbiasedMCstate); 
	saveMin(bestEnergy,bestUnbiasedMCstate,maxSavedResults);
	if (verbose) {
		cout << endl;
		cout << "Best Unbiased MC state: ";
		for (int j=0; j < bestUnbiasedMCstate.size(); j++){
			cout << bestUnbiasedMCstate[j] << ",";
		}
		cout << endl;
		cout << "Best Unbiased MC state Energy: " << bestEnergy << endl;
		cout << "===================================" << endl;
	}
}



double SelfPairManager::getInteractionEnergy(int _pos, int _rot, vector<unsigned int>& _currentState){

	double energy = 0.0;
	energy += computeSelfE(_pos,_rot); 
	for (uint i = 0;i <_currentState.size();i++){
		int pos2 = i;
		int rot2 = _currentState[i];
		if (_pos == pos2) continue;
		// IF _pos is LINKED THEN WHAT?
		// Then if pos2 is linked to _pos, we need to use _rot instead of rot2?
		if (_pos > pos2){
			energy += computePairE(_pos,_rot,pos2,rot2);
		} else {
			energy += computePairE(pos2,rot2,_pos,_rot);
		}
	}
	return energy;
}
int SelfPairManager::selectRandomStateAtPosition(int _position,vector<unsigned int>& _currentState) const {
	if (_position >= _currentState.size()) {
		cerr << "ERROR 72960: position " << _position << " out of range in int SelfPairManager::selectRandomStateAtPosition(int _position,vector<unsigned int>& _currentState) const " << endl;
		return 0;
	}
	if (selfE[_position].size() == 1) {
		return 0;
	}

	vector<double> residualP;
	double sumP = 0.0;
	for (int i=0; i<selfE[_position].size(); i++) {
		if (i == _currentState[_position] || !aliveMask[_position][i]) {
			residualP.push_back(0.0);
		} else {
			residualP.push_back(1.0);
		}
		sumP = residualP.back();
	}

	if (sumP == 0.0) {
		// no move avalable, stay there
		return _currentState[_position];
	}

	pRng->setDiscreteProb(residualP);
	return pRng->getRandomDiscreteIndex();
}

void SelfPairManager::getRandomState(vector<unsigned int>& _currentState) {
	for (int i=0; i<selfE.size(); i++) {
		_currentState[i] = selectRandomStateAtPosition(i,_currentState);
	}
}


void SelfPairManager::runGreedyOptimizer(int _cycles) {
	//Set up the mask so that every rotamer will be checked
	std::vector< std::vector<bool> > mask(getNumPositions());
	for (unsigned int i = 0; i < mask.size(); i++) {
		mask[i] = vector<bool>(getNumRotamers(i), true);
	}
	return runGreedyOptimizer(_cycles, mask);


}

void SelfPairManager::runGreedyOptimizer(int _cycles, std::vector< std::vector<bool> > _mask) {

	minBound.clear();
	minStates.clear();

	int currCycle = 0;
	int conMaxCycles=200; // it usually takes < 10 traversals to converge each time

	int numPositions = getNumPositions();

	// Lets pick the first rotamer in all positions for the first state
	vector<unsigned int> state(numPositions,0);
	vector<unsigned int> prevState;

	// Reset the aliveMask, we dont really need it, but the getRandomState needs it
	aliveMask.clear();
	for (unsigned int i=0; i < selfE.size(); i++) {
		aliveMask.push_back(vector<bool>(selfE[i].size(), true));
	}

	double energy=MslTools::doubleMax;
	vector<unsigned int> positionOrder;
	for(int i = 0; i < numPositions; i++) {
		positionOrder.push_back(i);
	}
	while(currCycle < _cycles){
		
		int convCycle = 0;
		
		while (convCycle < conMaxCycles){
			//Assignment of prevState
			prevState=state;
		
			//get a random order of positions
			random_shuffle(positionOrder.begin(),positionOrder.end());
			//Actual loop to calculate energy and modify state
			for (int i=0; i<positionOrder.size(); i++) {
				energy = MslTools::doubleMax;
				double thisEnergy = MslTools::doubleMax;
				for(int j=0;j<getNumRotamers(positionOrder[i]);j++){
					//Check if the rotamer has been masked out-- if so, we won't check it
					if(_mask[positionOrder[i]][j] == 1) {		
						thisEnergy = getInteractionEnergy(positionOrder[i],j,state);
						if ( thisEnergy < energy){
							energy = thisEnergy;
							state[positionOrder[i]]=j;
						}
					}
				}
			}
	
			//Checking for convergence 
			bool converge = true;
			for (int l=0;l<numPositions;l++) {
				if (state[l]!=prevState[l]) converge=false;
			}
		
			if (converge) break;	
			convCycle++;
		}
		double tmpE = getStateEnergy(state);
		saveMin(tmpE,state,maxSavedResults);
		if(verbose) {
			cout << "Cycle " << currCycle << " converged in " << convCycle << " cycles with Energy " << tmpE << endl;
			cout << "State " << currCycle << " ";
			for(int i = 0; i < state.size(); i++) {
				cout << state[i] << "," ;
			}
			cout << endl;
			
		}
		currCycle++;
		getRandomState(state);
	}
}

vector<unsigned int> SelfPairManager::runLP(bool _runMIP) {
#ifdef __GLPK__
		LinearProgrammingOptimization lpo;
		vector<vector<double> >& oligomersSelf = getSelfEnergy();
		vector<vector<vector<vector<double> > > >& oligomersPair = getPairEnergy();
		lpo.addEnergyTable(oligomersSelf,oligomersPair);
		lpo.setVerbose(verbose);
		return lpo.getSolution(_runMIP);
#else
		cerr << "GLPK library needs to be installed to run Linear Programming Optimization" << endl;
		return vector<unsigned int> ();
#endif
}

// In current state of system, get the energy for this term for this _posId
// E = selfE 
double SelfPairManager::computeSelfE(string _posId, string _resName, string _term){
  if (pSys == NULL || !pSys->positionExists(_posId)){
    cerr << "ERROR 343343 position: "<<_posId<<" does not exist in system or sytem is NULL.\n";
    exit(343343);
  }

  Position &pos = pSys->getPosition(_posId);
  int posIndex = variablePosIndex[&pos];

  //cout << "Variable index for position["<<_posId<<"] = "<<posIndex<<" number of rotamers: "<<variableCount[posIndex-1]<<" "<<variableCount.size()<<" "<<variableIdentities.size()<<" "<<variableIdentities[posIndex].size()<<" "<<selfE.size()<<" "<<selfE[0].size()<<" "<<selfEbyTerm.size()<<endl;

  // Iterate over rotamers at this position: variableCount[posIndex]
  double minE = MslTools::doubleMax;
  int minR = MslTools::intMax;
  for (uint i = 1; i < variableCount[posIndex-1];i++){

    // Only consider rotamers of the specific residue type
    int rotamerIndex = rotamerPos_Id_Rot[posIndex-1][i][1];

    if (variableIdentities[posIndex][rotamerIndex]->getResidueName() != _resName) continue;
    //cout << "Rotamer["<<i<<","<<rotamerIndex<<"]: is "<<variableIdentities[posIndex][rotamerIndex]->getResidueName()<<" looking for "<<_resName<<endl;

    double E = 0.0;
    E = selfE[posIndex-1][i];

    if (E < minE){
      minE = E;
      minR = i;
    }
  }

  if (_term != ""){
    minE = selfEbyTerm[posIndex-1][minR][_term];
  }

  minStates.clear();
  vector<unsigned int> state;
  state.push_back(minR);
  minStates.push_back(state);

  return minE;
  
  
}


// In current state of system, get the energy for this term for this _posId
// E = selfE 
double SelfPairManager::computeBestPairE(string _posId, string _resName, string _posId2, string _resName2, string _term){
  if (pSys == NULL || !pSys->positionExists(_posId)){
    cerr << "ERROR 343343 position: "<<_posId<<" does not exist in system or sytem is NULL.\n";
    exit(343343);
  }
  if (!pSys->positionExists(_posId2)){
    cerr << "ERROR 343343 position: "<<_posId2<<" does not exist in system or sytem is NULL.\n";
    exit(343344);
  }

  Position &pos1 = pSys->getPosition(_posId);
  for (uint a = 0; a < pos1.getAtomPointers().size();a++){
    pos1.getAtom(a).setSelectionFlag(pos1.getPositionId(),true);
  }

  Position &pos2 = pSys->getPosition(_posId2);

  for (uint a = 0; a < pos2.getAtomPointers().size();a++){
    pos2.getAtom(a).setSelectionFlag(pos2.getPositionId(),true);
  }
  if (_term != ""){
    pESet->setAllTermsInactive();
    pESet->setTermActive(_term);
  }

  // For each pair of rotamers
  double startE = pESet->getTotalEnergy();
  double minE = MslTools::doubleMax;
  int minR1   = MslTools::intMax;
  int minR2   = MslTools::intMax;
  for (uint r1 = 0; r1 < pos1.getTotalNumberOfRotamers(_resName);r1++){
    pos1.setActiveRotamer(_resName,r1);
    for (uint r2 = 0; r2 < pos2.getTotalNumberOfRotamers(_resName2);r2++){
      pos2.setActiveRotamer(_resName2,r2);
      double E = pESet->calcEnergy(pos1.getPositionId(), pos2.getPositionId());
      double totalE = pESet->getTotalEnergy();
      
      if (E < minE && (totalE - startE) < 50  ){
	minE = E;
	minR1 = r1;
	minR2 = r2;
      }
    }
  }

  pos1.setActiveRotamer(_resName,minR1);
  pos2.setActiveRotamer(_resName2,minR2);

  return minE;
}

