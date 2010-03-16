/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#include "CharmmSystemBuilder.h"

using namespace MSL;
using namespace std;


CharmmSystemBuilder::CharmmSystemBuilder() {
	setup();
}

CharmmSystemBuilder::CharmmSystemBuilder(string _topologyFile, string _parameterFile, string _solvationFile) {
	setup();
	readTopology(_topologyFile);
	readParameters(_parameterFile);
	if (_solvationFile != "") {
		readSolvation(_solvationFile);
	}
}

CharmmSystemBuilder::CharmmSystemBuilder(const CharmmSystemBuilder & _sysBuild) {
	setup();
	copy(_sysBuild);
}

CharmmSystemBuilder::~CharmmSystemBuilder() {
	delete pTopReader;
	delete pParReader;
	delete pEEF1ParReader;
}

void CharmmSystemBuilder::operator=(const CharmmSystemBuilder & _sysBuild) {
	copy(_sysBuild);
}


void CharmmSystemBuilder::setup() {
	pTopReader = new CharmmTopologyReader;
	pParReader = new CharmmParameterReader;
	pEEF1ParReader = new CharmmEEF1ParameterReader;
	vdwRescalingFactor = 1.0;
	buildNonBondedInteractions = true;
	elec14factor = 1;
	dielectricConstant = 1;
	useRdielectric = true;
	useSolvation = false;
	solvent = pEEF1ParReader->getDefaultSolvent();
	useSolvation = false;
	solvent = pEEF1ParReader->getDefaultSolvent();
}

void CharmmSystemBuilder::copy(const CharmmSystemBuilder & _sysBuild) {
	*pTopReader = *_sysBuild.pTopReader;
	*pParReader = *_sysBuild.pParReader;
	*pEEF1ParReader = *_sysBuild.pEEF1ParReader;
	buildNonBondedInteractions = _sysBuild.buildNonBondedInteractions;
	vdwRescalingFactor = _sysBuild.vdwRescalingFactor;
	elec14factor = _sysBuild.elec14factor;
	dielectricConstant = _sysBuild.dielectricConstant;
	useRdielectric = _sysBuild.useRdielectric;
}

void CharmmSystemBuilder::deletePointers() {
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator k=polymerDefi.begin(); k!=polymerDefi.end(); k++) {
		for (vector<vector<CharmmTopologyResidue*> >::iterator l=k->begin(); l!=k->end(); l++) {
			for (vector<CharmmTopologyResidue*>::iterator m=l->begin(); m!=l->end(); m++) {
				delete *m;
			}
		}
	}
	polymerDefi.clear();
	atomMap.clear();
}

bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence) {
	deletePointers();
	
	vector<vector<vector<string> > > seq = _sequence.getSequence();

	for (vector<vector<vector<string> > >::iterator k=seq.begin(); k!=seq.end(); k++) {
		// for each chain
		polymerDefi.push_back(vector<vector<CharmmTopologyResidue*> >());
		atomMap.push_back(vector<vector<map<string, Atom*> > >());

		for (vector<vector<string> >::iterator l=k->begin(); l!=k->end(); l++) {
			// for each position
			polymerDefi.back().push_back(vector<CharmmTopologyResidue*>());
			atomMap.back().push_back(vector<map<string, Atom*> >());

			for (vector<string>::iterator m=l->begin(); m!=l->end(); m++) {
				// each identity

				vector<string> split = MslTools::tokenize(*m, "-");
				// dangerously assuming that the initial string isn't blank!
				if (pTopReader->residueExists(split[0])) {
					polymerDefi.back().back().push_back(new CharmmTopologyResidue(pTopReader->getLastFoundResidue()));
					atomMap.back().back().push_back(map<string, Atom*>());
				} else {
					cerr << "WARNING 19129: residue " << split[0] << " does not exist in topology, in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
					return false;
				}
				// if the patches for the first or the last residue were not declared
				// get them from the residues
				if (split.size() == 1 && l==k->begin()) {
					// add default first patch to current identity
					split.push_back(polymerDefi.back().back().back()->getFirstDefaultPatch());
				}
				if (split.size() == 1 && l==k->end() - 1) {
					// add default last patch to current identity
					split.push_back(polymerDefi.back().back().back()->getLastDefaultPatch());
				}
				// apply patches (if any)
				for (vector<string>::iterator n=split.begin()+1; n!=split.end(); n++) {
					// patch current identity
					if (pTopReader->residueExists(*n)) {
						if (!polymerDefi.back().back().back()->applyPatch(pTopReader->getLastFoundResidue())) {
							cerr << "WARNING 19134: cannot apply patch " << *n << " to residue, in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
							return false;
						}
					}
				}
			//	cout << polymerDefi.back().back().back()->toString();

			}
		}
	}
	AtomPointerVector sysAtoms;
	string name;
	string type;
	double partialCharge;
	string element;
	int group;
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator i=polymerDefi.begin(); i!=polymerDefi.end(); i++) {
		string chainId = _sequence.getChainId(i-polymerDefi.begin());
	//	cout << "UUU chain " << chainId << endl;
		for (vector<vector<CharmmTopologyResidue*> >::iterator j=i->begin(); j!=i->end(); j++) {
			string resNumStr = _sequence.getResidueNumber(i-polymerDefi.begin(), j-i->begin());
			int resNum;
			string iCode;
			MslTools::splitIntAndString(resNumStr, resNum, iCode);
	//		cout << "UUU chain/resnum " << chainId << "/" << resNum << iCode << endl;
			for (vector<CharmmTopologyResidue*>::iterator k=j->begin(); k!=j->end(); k++) {
	//			cout << "UUUG " << (*k)->toString() << endl;
				string resName = (*k)->getName();
	//			cout << "UUU chain/resnum resname " << chainId << "/" << resNum << " " << resName << endl;
				for (unsigned int l=0; l<(*k)->atomSize(); l++) {
					(*k)->getTopolAtom(l, name, type, partialCharge, element, group);
	//				cout << "UUU chain/resnum resname " << chainId << "/" << resNum << " " << resName << " " << name << endl;
					Atom * a = new Atom(name);
					a->setChainId(chainId);
					a->setResidueName(resName);
					a->setResidueNumber(resNum);
					a->setResidueIcode(iCode);
					a->setCharge(partialCharge);
					a->setType(type);
					a->setElement(element);
					a->setGroupNumber(group);
					a->wipeCoordinates();
					sysAtoms.push_back(a);
	//				cout << "UUU " << *a << endl;
				};
				//unsigned int _i, string & _name, string & _type, string _partialCharge, string & _element, unsigned int & _group
			}
		}
	}
//	cout << "UUU add " << sysAtoms.size() << " to system" << endl;
//	for (AtomPointerVector::iterator k=sysAtoms.begin(); k!=sysAtoms.end(); k++) {
//		cout << **k << endl;
//	}
	_system.addAtoms(sysAtoms);
//	AtomPointerVector fromSys = _system.getAtoms();
//	cout << "UUU System created with " << _system.atomSize() << " atoms" << endl;
//	if (_system.atomSize() != sysAtoms.size()) {
//		// the number of atoms is not consistent
//		//cerr << "WARNING..." << endl;
//		//_system.reset();
//		return false;
//	}
	// garbage collection
	for (AtomPointerVector::iterator k=sysAtoms.begin(); k!=sysAtoms.end(); k++) {
		delete *k;
	}

	/************************************************************
	 *  Create the atoms map
	 ************************************************************/
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator k=polymerDefi.begin(); k!=polymerDefi.end(); k++) {
		string chainId = _sequence.getChainId(k-polymerDefi.begin());
		for (vector<vector<CharmmTopologyResidue*> >::iterator l=k->begin(); l!=k->end(); l++) {
			string resNumStr = _sequence.getResidueNumber(k-polymerDefi.begin(), l-k->begin());
			int resNum;
			string iCode;
			MslTools::splitIntAndString(resNumStr, resNum, iCode);
			for (vector<CharmmTopologyResidue*>::iterator m=l->begin(); m!=l->end(); m++) {
//				cout << (*m)->getName() << " " << resNumStr << " " << chainId << endl;
				string resName = (*m)->getName();
				for (unsigned int i=0; i<(*m)->atomSize(); i++) {
					(*m)->getTopolAtom(i, name, type, partialCharge, element, group);
		//			cout << name << endl;
					//if (_system.exists(chainId, resNumStr, name, resName)) {
					if (_system.atomExists(chainId, resNum, iCode, resName, name)) {
						// a vector/vector/vector/map(string, Atom*)
						atomMap[k-polymerDefi.begin()][l-k->begin()][m-l->begin()][name] = &(_system.getLastFoundAtom());
					} else {
						//cerr << "WARNING..." << endl;
						//_system.reset();
						return false;
					}
				}
			}
		}
	}

	/************************************************************
	 *  Create the atoms map
	 ************************************************************/
	vector<string> icAtoms;
	vector<double> icValues;
	bool improperFlag;
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator k=polymerDefi.begin(); k!=polymerDefi.end(); k++) {
		// for each chain
		for (vector<vector<CharmmTopologyResidue*> >::iterator l=k->begin(); l!=k->end(); l++) {
			// for each postion
			for (vector<CharmmTopologyResidue*>::iterator m=l->begin(); m!=l->end(); m++) {
				// for each identity
	//			cout << "UUUU " << (*m)->icSize() << endl;
				for (unsigned int i=0; i<(*m)->icSize(); i++) {
					(*m)->getIcLine(i, icAtoms, icValues, improperFlag);
					vector<vector<Atom*> > icAtomPointers(4, vector<Atom*>());
					for (vector<string>::iterator n=icAtoms.begin(); n!=icAtoms.end(); n++) {
						if (n->substr(0,1) == "-") {
							if (l > k->begin()) {
								// not the first residue of the chain
								vector<vector<CharmmTopologyResidue*> >::iterator ll=l-1;
								for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
									// for each identity get all the -X atoms

//									cout << "UUUX " << atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].size() << endl;
									map<string, Atom*>::iterator found = atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
									if (found == atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].end()) {
										icAtomPointers[n-icAtoms.begin()].push_back(NULL);
//										cout << "UUU atom -" << n->substr(1, n->size()-1) << " NOT found" << endl;
									} else {
										icAtomPointers[n-icAtoms.begin()].push_back(found->second);
//										cout << "UUU atom -" << n->substr(1, n->size()-1) << " found (" << icAtomPointers[n-icAtoms.begin()].size() << ")" << endl;
									}

								}
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
							}
							//cout << "UUU minus " << *n << endl;
						} else if (n->substr(0,1) == "+") {
							if (l<k->end()-1) {
								// not the last residue of the chain
								vector<vector<CharmmTopologyResidue*> >::iterator ll=l+1;
								for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
									// for each identity get all the -X atoms

//									cout << "UUUY " << atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].size() << endl;
									map<string, Atom*>::iterator found = atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
									if (found == atomMap[k-polymerDefi.begin()][ll-k->begin()][o-ll->begin()].end()) {
										icAtomPointers[n-icAtoms.begin()].push_back(NULL);
//										cout << "UUU atom +" << n->substr(1, n->size()-1) << " NOT found" << endl;
									} else {
										icAtomPointers[n-icAtoms.begin()].push_back(found->second);
//										cout << "UUU atom +" << n->substr(1, n->size()-1) << " found (" << icAtomPointers[n-icAtoms.begin()].size() << ")" << endl;
									}

								}
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
							}
							//cout << "UUU plus " << *n << endl;
						} else {
							map<string, Atom*>::iterator found = atomMap[k-polymerDefi.begin()][l-k->begin()][m-l->begin()].find(*n);
							if (found == atomMap[k-polymerDefi.begin()][l-k->begin()][m-l->begin()].end()) {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
//								cout << "UUU atom " << *n << " NOT found" << endl;
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(found->second);
//								cout << "UUU atom " << *n << " found (" << icAtomPointers[n-icAtoms.begin()].size() << ")" << endl;
							}
						}
					}
					// add the IC terms to the system
					//vector<vector<Atom*> > icAtomPointers(4, vector<Atom*>());
					for (vector<Atom*>::iterator k=icAtomPointers[0].begin(); k!=icAtomPointers[0].end(); k++) {
						for (vector<Atom*>::iterator l=icAtomPointers[1].begin(); l!=icAtomPointers[1].end(); l++) {
							for (vector<Atom*>::iterator m=icAtomPointers[2].begin(); m!=icAtomPointers[2].end(); m++) {
								for (vector<Atom*>::iterator n=icAtomPointers[3].begin(); n!=icAtomPointers[3].end(); n++) {
									if ((*k != NULL || *n != NULL) && *l != NULL && *m != NULL) {
										
										// If icValues are zero (usually zero for the PatchResidues' atoms, get the minimum values from the parameter file
										if (icValues[0] == 0.0) {
											if(improperFlag) {
												if(*k) {
													//icValues[0] = (pParReader->bondParam((*k)->getType(),(*m)->getType()))[1];
													vector<double> param;
													if (pParReader->bondParam(param, (*k)->getType(),(*m)->getType())) {
														icValues[0] = param[1];
													}
												}
											} else {
												if(*k) {
													//icValues[0] = (pParReader->bondParam((*k)->getType(),(*l)->getType()))[1];
													vector<double> param;
													if (pParReader->bondParam(param, (*k)->getType(),(*l)->getType())) {
														icValues[0] = param[1];
													}
												}
											}
										}	

										if (icValues[1] == 0.0) {
											if(improperFlag) {
												if(*k) {
													//icValues[1] = (pParReader->angleParam((*k)->getType(),(*m)->getType(),(*l)->getType()))[1];
													vector<double> param;
													if (pParReader->angleParam(param, (*k)->getType(),(*m)->getType(),(*l)->getType())) {
														icValues[1] = param[1];
													}
												}
											} else {
												if(*k) {
													//icValues[1] = (pParReader->angleParam((*k)->getType(),(*l)->getType(),(*m)->getType()))[1];
													vector<double> param;
													if (pParReader->angleParam(param, (*k)->getType(),(*l)->getType(),(*m)->getType())) {
														icValues[1] = param[1];
													}
												}
											}
										}
										
										if (icValues[3] == 0.0) {
											if(*n) {
												//icValues[3] = (pParReader->angleParam((*l)->getType(),(*m)->getType(),(*n)->getType()))[1];
												vector<double> param;
												if (pParReader->angleParam(param, (*l)->getType(),(*m)->getType(),(*n)->getType())) {
													icValues[3] = param[1];
												}
											}
										}
										
										if (icValues[4] == 0.0) {
											if(*n) {
												//icValues[4] = (pParReader->bondParam((*m)->getType(),(*n)->getType()))[1];
												vector<double> param;
												if (pParReader->bondParam(param, (*m)->getType(),(*n)->getType())) {
													icValues[4] = param[1];
												}
											}
										}

										_system.addIcEntry(*k, *l, *m, *n, icValues[0], icValues[1], icValues[2], icValues[3], icValues[4], improperFlag);
									}
								}
							}
						}
					}


				}
			}
		}
	}
	

	EnergySet* ESet = _system.getEnergySet();
        
	
	// Loop over each chain.
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator chItr=polymerDefi.begin(); chItr!=polymerDefi.end(); chItr++) {
		string chainId = _sequence.getChainId(chItr-polymerDefi.begin());
		
		// Loop over each position.
		for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=chItr->begin(); posItr!=chItr->end(); posItr++) {
			string resNumStr1 = _sequence.getResidueNumber(chItr-polymerDefi.begin(), posItr-chItr->begin());
		
			// Loop over each identity.
			for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin(); idItr!=posItr->end(); idItr++) {
				string resName1 = (*idItr)->getName();

				// add the bonds
				for (unsigned int l=0; l<(*idItr)->bondSize(); l++) {
					// get the atoms names
					string atom1,atom2;
					unsigned int type;
					(*idItr)->getBond(l, atom1, atom2, type);           
				
					/***********************************************************************
					 *   If multiple identities are present, and some atoms refers to previous 
					 *   (i.e. -N) or next (+N) residue, then the bonds should be added in
					 *   all combinations
					 ***********************************************************************/
				
					vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr);
					vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
					
					// Now loop over all the pAtom1 and pAtom2 atoms and add an interaction for each combination
					for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
						string type1 = (*a1)->getType();
						for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
							string type2 = (*a2)->getType();
							(*a1)->setBoundTo(*a2);
							//vector<double> params = pParReader->bondParam(type1, type2);	
							vector<double> params;
							if (pParReader->bondParam(params, type1, type2)) {
								CharmmBondInteraction *pCBI = new CharmmBondInteraction(*(*a1),*(*a2),params[0],params[1]);
								ESet->addInteraction(pCBI);
							} else {
								cerr << "WARNING: CHARMM bond parameters not found for atom types " << type1 << ", " << type2 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
							}
						//	if (params[0] != 0.0) {
						//		CharmmBondInteraction *pCBI = new CharmmBondInteraction(*(*a1),*(*a2),params[0],params[1]);
						//		ESet->addInteraction(pCBI);
						//	}
						}
					}
						
				}
				/************************* DONE SETTING UP THE BONDS **********************************/

				// add the angles if they are not autogenerated
				if (!pTopReader->getAutoGenerateAngles()) {
					for (unsigned int l=0; l<(*idItr)->angleSize(); l++) {
						// Loop over each improper in the residue.
						
						string atom1,atom2,atom3;
						// get the atoms names
						(*idItr)->getAngle(l, atom1, atom2, atom3);
					
						vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
						vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
						vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
						
						for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
							string type1 = (*a1)->getType();
							for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
								string type2 = (*a2)->getType();
								for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
									string type3 = (*a3)->getType();
									//vector<double> params = pParReader->ureyBradleyParam(type1, type2, type3);
									vector<double> params;
									if (pParReader->ureyBradleyParam(params, type1, type2, type3)) {
										CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*a1),*(*a3),params[0],params[1]);
										ESet->addInteraction(pCUI);
							//		} else {
							//			cerr << "WARNING: CHARMM urey bradley parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
									}
								//	if (params[0] != 0.0) {
								//		CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*a1),*(*a3),params[0],params[1]);
								//		ESet->addInteraction(pCUI);
								//	}
									//params = pParReader->angleParam(type1, type2, type3);
									if (pParReader->angleParam(params, type1, type2, type3)) {
										CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*a1),*(*a2),*(*a3),params[0],params[1] * M_PI / 180.0); 
										ESet->addInteraction(pCAI);
									} else {
										cerr << "WARNING: CHARMM angle parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
									}
								//	if (params[0] != 0.0) {
								//		CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*a1),*(*a2),*(*a3),params[0],params[1] * M_PI / 180.0); 
								//		ESet->addInteraction(pCAI);
								//	}
								}
							}
						}
							
					}
				}
				/************************* DONE SETTING UP THE ANGLES **********************************/

				// add the dihedrals if they are not autogenerated
				if (!pTopReader->getAutoGenerateDihedrals()) {
					//for (unsigned int l=0; l<(*idItr)->improperSize(); l++) {
					for (unsigned int l=0; l<(*idItr)->dihedralSize(); l++) {
						// Loop over each improper in the residue.
						
						string atom1,atom2,atom3,atom4;
						// get the atoms names
						(*idItr)->getDihedral(l, atom1, atom2, atom3, atom4);
					
						vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
						vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
						vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
						vector<Atom*> pAtom4 = getAtomPointers(atom4, chItr, posItr, idItr);
						
						for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
							string type1 = (*a1)->getType();
							for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
								string type2 = (*a2)->getType();
								for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
									string type3 = (*a3)->getType();
									for(vector<Atom*>::iterator a4 = pAtom4.begin() ; a4 != pAtom4.end(); a4++) {
										string type4 = (*a4)->getType();

										//vector<vector <double> > dihedralEntries = pParReader->dihedralParam(type1, type2, type3, type4);
										vector<vector <double> > dihedralEntries;
										if (pParReader->dihedralParam(dihedralEntries, type1, type2, type3, type4)) {
											// there could be multiple entries for a single dihedral
											for(int m = 0; m < dihedralEntries.size() ; m++) {
												// the delta should be expressed in radians
												dihedralEntries[m][2] = dihedralEntries[m][2]*M_PI/180.0;

											}
											CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*a1),*(*a2),*(*a3),*(*a4),dihedralEntries);
											ESet->addInteraction(pCDI);
										} else {
											cerr << "WARNING: CHARMM dihedral parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << ", " << type4 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
										}
									
									}
								}
							}
						}
							
					}
				}
				/************************* DONE SETTING UP THE DIHEDRALS **********************************/

				// add the impropers
				for (unsigned int l=0; l<(*idItr)->improperSize(); l++) {
					// Loop over each improper in the residue.
					
					string atom1,atom2,atom3,atom4;
					// get the atoms names
					(*idItr)->getImproper(l, atom1, atom2, atom3, atom4);
				
					vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
					vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
					vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
					vector<Atom*> pAtom4 = getAtomPointers(atom4, chItr, posItr, idItr);
					
					for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
						string type1 = (*a1)->getType();
						for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
							string type2 = (*a2)->getType();
							for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
								string type3 = (*a3)->getType();
								for(vector<Atom*>::iterator a4 = pAtom4.begin() ; a4 != pAtom4.end(); a4++) {
									string type4 = (*a4)->getType();
									//vector<double> improperParams = pParReader->improperParam(type1, type2, type3, type4);	
									vector<double> improperParams;
									if (pParReader->improperParam(improperParams, type1, type2, type3, type4)) {
										CharmmImproperInteraction *pCII = new CharmmImproperInteraction(*(*a1),*(*a2),*(*a3),*(*a4),improperParams[0],improperParams[1]*M_PI/180.0);
										ESet->addInteraction(pCII);
									} else {
										cerr << "WARNING: CHARMM improper parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << ", " << type4 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
									}
								}
							}
						}
					}
						
				}
				/************************* DONE SETTING UP THE IMPROPERS **********************************/
			}
		}
	}

	
	AtomPointerVector atoms = _system.getAllAtoms();

	/*********************************************************************************
	 *
	 *  ADD THE INTERACTION TERMS:
	 *
	 *  For each pair of atoms, add the appropriate interactions:
	 *   - bond if 1-2
	 *   - angle if 1-3
	 *   - dihedral if 1-4
	 *   - special vdw and elec with e14fac if 1-4 (and NOT also 1-3 or 1-2, it is possible)
	 *   - otherwise regular vdw and elec
	 **********************************************************************************/
	for(AtomPointerVector::iterator atomI = atoms.begin(); atomI < atoms.end(); atomI++) {
		//cout << "UUU A atomI " << **atomI << endl;
		string atomItype = (*atomI)->getType();
		//vector<double> vdwParamsI = pParReader->vdwParam(atomItype);

		for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
			if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
				continue;
			}
			string atomJtype = (*atomJ)->getType();
			bool special = false;
			if ((*atomI)->isBoundTo(*atomJ)) {
				special = true;
			}
			if ((*atomI)->isOneThree(*atomJ)) {
				if (pTopReader->getAutoGenerateAngles()) {
					// autogenerate the angle interactions
					vector<Atom*> middle = (*atomI)->getOneThreeMiddleAtoms(*atomJ);
						
					for (vector<Atom*>::iterator k=middle.begin(); k!=middle.end(); k++) {
						string middleType = (*k)->getType();
						//vector<double> params = pParReader->ureyBradleyParam(atomItype, middleType, atomJtype);
						vector<double> params;
						if (pParReader->ureyBradleyParam(params, atomItype, middleType, atomJtype)) {
							CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*atomI),*(*atomJ),params[0],params[1]);
							ESet->addInteraction(pCUI);
//						} else {
//							cerr << "WARNING: CHARMM urey bradley parameters not found for atom types " << atomItype << ", " << middleType << ", " << atomJtype << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
						}
					//	if (params[0] != 0.0) {
					//		CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*atomI),*(*atomJ),params[0],params[1]);
					//		ESet->addInteraction(pCUI);
					//	}
						//params = pParReader->angleParam(atomItype, middleType, atomJtype);
						if (pParReader->angleParam(params, atomItype, middleType, atomJtype)) {
							CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*atomI),*(*k),*(*atomJ),params[0],params[1] * M_PI / 180.0); 
							ESet->addInteraction(pCAI);
						} else {
							cerr << "WARNING: CHARMM angle parameters not found for atom types " << atomItype << ", " << middleType << ", " << atomJtype << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
						}
					//	if (params[0] != 0.0) {
					//		CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*atomI),*(*k),*(*atomJ),params[0],params[1] * M_PI / 180.0); 
					//		ESet->addInteraction(pCAI);
					//	}
					}
				}
				special = true;
			}
			if ((*atomI)->isOneFour(*atomJ)) {
				if (pTopReader->getAutoGenerateDihedrals()) {
					// autogenerate the dihedral interactions
					vector<vector<Atom*> > middle = (*atomI)->getOneFourMiddleAtoms(*atomJ);
					for (vector<vector<Atom*> >::iterator k=middle.begin(); k!=middle.end(); k++) {
					//	cout <<  "UUU " << "Add dihe " << k-middle.begin() << ": " << (*atomI)->getName();
						vector<Atom*> middleAtoms;
						vector<string> middleTypes;
						for (vector<Atom*>::iterator l=k->begin(); l!=k->end(); l++) {
							middleAtoms.push_back(*l);
							middleTypes.push_back((*l)->getType());
					//		cout <<  "UUU " << " " << (*l)->getName();
						}
					//	cout <<  "UUU " << " " << (*atomJ)->getName() << endl;

						//vector<vector <double> > dihedralEntries = pParReader->dihedralParam(atomItype, middleTypes[0], middleTypes[1], atomJtype);
						vector<vector <double> > dihedralEntries;
						if (pParReader->dihedralParam(dihedralEntries, atomItype, middleTypes[0], middleTypes[1], atomJtype)) {
							// there could be multiple entries for a single dihedral
							for(int m = 0; m < dihedralEntries.size() ; m++) {
								// the delta should be expressed in radians
								dihedralEntries[m][2] = dihedralEntries[m][2]*M_PI/180.0;

							}
						
							CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*atomI),*(middleAtoms[0]),*(middleAtoms[1]),*(*atomJ),dihedralEntries);
							ESet->addInteraction(pCDI);
						} else {
							cerr << "WARNING: CHARMM dihedral parameters not found for atom types " << atomItype << ", " << middleTypes[0] << ", " << middleTypes[1] << ", " << atomJtype << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
						}
					}
						
				}
				/*
				if (buildNonBondedInteractions && !special) {
					// if it is also 1-3 or 1-2 do not add the vdw and elec term
					vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
					CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[3]+vdwParamsJ[3]) * vdwRescalingFactor, sqrt(vdwParamsI[2] * vdwParamsJ[2]) );
					CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant,elec14factor, useRdielectric);
				//	pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
				//	pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
					ESet->addInteraction(pCVI);
					ESet->addInteraction(pCEI);
				}
				special = true;
				*/
			}
			/*
			// Remove any non-bonded if flag is set
			if (buildNonBondedInteractions && !special) {
				vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
				CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ),(vdwParamsI[1]+vdwParamsJ[1]) * vdwRescalingFactor, sqrt(vdwParamsI[0] * vdwParamsJ[0]) );
				CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant, 1.0, useRdielectric);
			//	pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
			//	pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
				ESet->addInteraction(pCVI);
				ESet->addInteraction(pCEI);

			}
			*/
		}
	}
	if (buildNonBondedInteractions) {
		updateNonBonded(_system, 0.0, 0.0, 0.0);
	}

	return true;
}

bool CharmmSystemBuilder::updateNonBonded(System & _system, double _ctonnb, double _ctofnb, double _cutnb) {
	/********************************************************************************
	 *  About the cutoffs:
	 *
	 *   - if _cutnb is not zero, a distance cutoff is applied to exclude interactions between far atoms
	 *   - the _ctonnb is the cutoff in which the switching function is applied to bring the energy
	 *     smoothly to zero.  
	 *   - the energy goes to zero at _ctofnb
	 *  That is:
	 *   - between 0 and _ctonnb E = full energy
	 *   - between _ctonnb and _ctofnb E = energy * switching function
	 *   - between _ctofnb and _cutnb E = 0.0
	 *   - between _cutnb and infinity, the interaction is not in the list
	 ********************************************************************************/
	/********************************************************************************
	 *  Lazaridis solvation EEF1:
	 *
	 *  Includes a single body term (reference energy, Gref, CharmmEEF1RefInteraction) 
	 *  and a two body term CharmmEEF1Interaction
	 *
	 ********************************************************************************/
	EnergySet* ESet = _system.getEnergySet();
	ESet->resetTerm("CHARMM_VDW");
	ESet->resetTerm("CHARMM_ELEC");
	ESet->resetTerm("CHARMM_EFF1");
	ESet->resetTerm("CHARMM_EF1R");
	AtomPointerVector atoms = _system.getAllAtoms();

	bool useSolvation_local = useSolvation;
	if (!pEEF1ParReader->solventExists(solvent)) {
		useSolvation_local = false;
	}

	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 **********************************************************************/
	unsigned int stamp = MslTools::getRandomInt(1000000);

	/*********************************************************************************
	 *
	 *  ADD THE NON-BONDED INTERACTION TERMS:
	 *   - special vdw and elec with e14fac if 1-4 (and NOT also 1-3 or 1-2, it is possible)
	 *   - otherwise regular vdw and elec
	 *
	 *     EEF1Param
	 *         0 = V    
	 *         1 = Gref 
	 *         2 = Gfree
	 *         3 = Href 
	 *         4 = CPref
	 *         5 = Sigw 
	 *                            0             2                3             vdw 0
	 *    void setParams(double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j);
	 **********************************************************************************/
	for(AtomPointerVector::iterator atomI = atoms.begin(); atomI < atoms.end(); atomI++) {
		if (_cutnb > 0.0 && !(*atomI)->hasCoor()) {
			// no coordinates, skip this atom
			continue;
		}
		string atomItype = (*atomI)->getType();
		bool foundVdw = true;
		//vector<double> vdwParamsI = pParReader->vdwParam(atomItype);
		vector<double> vdwParamsI;
		if (!pParReader->vdwParam(vdwParamsI, atomItype)) {
			cerr << "WARNING 49319: VDW parameters not found for type " << atomItype << " in bool CharmmSystemBuilder::updateNonBonded(System & _system, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
			foundVdw = false;
		}
		bool foundEEF1 = true; 
		vector<double> EEF1ParamsI; 
		if (useSolvation_local) {
			if (!pEEF1ParReader->EEF1Param(EEF1ParamsI, atomItype, solvent)) {
				cerr << "WARNING 49387: EEF1 parameters not found for type " << atomItype << ", solvent " << solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
				foundEEF1 = false;
			} else {
				// add the single body term
				CharmmEEF1RefInteraction *pCERI = new CharmmEEF1RefInteraction(*(*atomI),EEF1ParamsI[1]);
				ESet->addInteraction(pCERI);
			}
		}

		for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
			if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
				continue;
			}
			if (_cutnb > 0.0 && (!(*atomJ)->hasCoor() || (*atomI)->groupDistance(**atomJ, stamp) > _cutnb)) {
				// no coordinates or atom further away from cutoff, skip this atom
				continue;
			}

			string atomJtype = (*atomJ)->getType();
			bool special = false;
			if ((*atomI)->isBoundTo(*atomJ) || (*atomI)->isOneThree(*atomJ)) {
				special = true;
			}
			bool foundVdw2 = true;
			vector<double> vdwParamsJ;;
			if (!pParReader->vdwParam(vdwParamsJ, atomJtype)) {
				cerr << "WARNING 49319: VDW parameters not found for type " << atomJtype << " in bool CharmmSystemBuilder::updateNonBonded(System & _system, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
				foundVdw2 = false;
				//continue;
			}
			if ((*atomI)->isOneFour(*atomJ)) {
				if (!special) {
					// if it is also 1-3 or 1-2 do not add the vdw and elec term
					//vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
					CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant,elec14factor, useRdielectric);
					if (_cutnb > 0.0) {
						// if we are using a cutoff, set the Charmm VDW interaction with
						// the cutoffs for the switching function
						pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
					} else {
						pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
					}
					ESet->addInteraction(pCEI);
					if (foundVdw && foundVdw2) {
						CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[3]+vdwParamsJ[3]) * vdwRescalingFactor, sqrt(vdwParamsI[2] * vdwParamsJ[2]) );
						if (_cutnb > 0.0) {
							// if we are using a cutoff, set the Charmm VDW interaction with
							// the cutoffs for the switching function
							pCVI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
						} else {
							pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
						}
						ESet->addInteraction(pCVI);
						if (useSolvation_local && foundEEF1) {
							vector<double> EEF1ParamsJ;
							if (!pEEF1ParReader->EEF1Param(EEF1ParamsJ, atomJtype, solvent)) {
								cerr << "WARNING 49387: EEF1 parameters not found for type " << atomJtype << ", solvent " << solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
								continue;
							}
							//if (EEF1ParamsI[0] != 0.0 && EEF1ParamsI[2] != 0.0 && EEF3ParamsI[3] != 0.0 && EEF1ParamsJ[0] != 0.0 && EEF1ParamsJ[2] != 0.0 && EEF1ParamsJ[3] != 0.0) {
								CharmmEEF1Interaction *pCSI = new CharmmEEF1Interaction(*(*atomI),*(*atomJ), EEF1ParamsI[0], EEF1ParamsI[2], EEF1ParamsI[5], vdwParamsI[1], EEF1ParamsJ[0], EEF1ParamsJ[2], EEF1ParamsJ[5], vdwParamsJ[1]);
								if (_cutnb > 0.0) {
									// if we are using a cutoff, set the Charmm VDW interaction with
									// the cutoffs for the switching function
									pCSI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
								} else {
									pCSI->setUseNonBondCutoffs(false, 0.0, 0.0);
								}
								ESet->addInteraction(pCSI);
							//}
						}

					}
				}
				special = true;
			}
			if (!special) {
//				vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
				CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant, 1.0, useRdielectric);
				if (_cutnb > 0.0) {
					// if we are using a cutoff, set the Charmm VDW interaction with
					// the cutoffs for the switching function
					pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
				} else {
					pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
				}
				ESet->addInteraction(pCEI);
				if (foundVdw && foundVdw2) {
					CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[1]+vdwParamsJ[1]) * vdwRescalingFactor, sqrt(vdwParamsI[0] * vdwParamsJ[0]) );
					if (_cutnb > 0.0) {
						// if we are using a cutoff, set the Charmm VDW interaction with
						// the cutoffs for the switching function
						pCVI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
					} else {
						pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
					}
					ESet->addInteraction(pCVI);
					if (useSolvation_local && foundEEF1) {
						vector<double> EEF1ParamsJ;
						if (!pEEF1ParReader->EEF1Param(EEF1ParamsJ, atomJtype, solvent)) {
							cerr << "WARNING 49387: EEF1 parameters not found for type " << atomJtype << ", solvent " << solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
							continue;
						}
						//if (EEF1ParamsI[0] != 0.0 && EEF1ParamsI[2] != 0.0 && EEF1ParamsI[3] != 0.0 && EEF1ParamsJ[0] != 0.0 && EEF1ParamsJ[2] != 0.0 && EEF1ParamsJ[3] != 0.0) {
							CharmmEEF1Interaction *pCSI = new CharmmEEF1Interaction(*(*atomI),*(*atomJ), EEF1ParamsI[0], EEF1ParamsI[2], EEF1ParamsI[5], vdwParamsI[1], EEF1ParamsJ[0], EEF1ParamsJ[2], EEF1ParamsJ[5], vdwParamsJ[1]);
							if (_cutnb > 0.0) {
								// if we are using a cutoff, set the Charmm VDW interaction with
								// the cutoffs for the switching function
								pCSI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
							} else {
								pCSI->setUseNonBondCutoffs(false, 0.0, 0.0);
							}
							ESet->addInteraction(pCSI);
						//}
					}
				}
			}
		}
	}

	return true;
}


/*
bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb) {
	/ ********************************************************************************
	 *  Lazaridis solvation EEF1:
	 *
	 *  Includes a single body term (reference energy, Gref, CharmmEEF1RefInteraction) 
	 *  and a two body term CharmmEEF1Interaction
	 *
	 ******************************************************************************** /
	EnergySet* ESet = _system.getEnergySet();
	ESet->resetTerm("CHARMM_EFF1");
	ESet->resetTerm("CHARMM_EEF1REF");

	if (!pEEF1ParReader->solventExists(_solvent)) {
		return false;
	}

	AtomPointerVector atoms = _system.getAllAtoms();

	/ **********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 ********************************************************************** /
	unsigned int stamp = MslTools::getRandomInt(1000000);

	/ *********************************************************************************
	 *
	 *     EEF1Param
	 *         0 = V    
	 *         1 = Gref 
	 *         2 = Gfree
	 *         3 = Href 
	 *         4 = CPref
	 *         5 = Sigw 
	 *                            0             2                3             vdw 0
	 *    void setParams(double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j);
	 ********************************************************************************** /
	for(AtomPointerVector::iterator atomI = atoms.begin(); atomI < atoms.end(); atomI++) {
		if (_cutnb > 0.0 && !(*atomI)->hasCoor()) {
			// no coordinates, skip this atom
			continue;
		}
		string atomItype = (*atomI)->getType();
	//	vector<double> vdwParamsI = pParReader->vdwParam(atomItype);
		vector<double> vdwParamsI;
		if (!pParReader->vdwParam(vdwParamsI, atomItype)) {
			cerr << "WARNING 49319: VDW parameters not found for type " << atomItype << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
			continue;
		}
		vector<double> EEF1ParamsI; 
		if (!pEEF1ParReader->EEF1Param(EEF1ParamsI, atomItype, _solvent)) {
			cerr << "WARNING 49387: EEF1 parameters not found for type " << atomItype << ", solvent " << _solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
			continue;
		}
		
		// add the single body term
		CharmmEEF1RefInteraction *pCERI = new CharmmEEF1RefInteraction(*(*atomI),EEF1ParamsI[1]);
		ESet->addInteraction(pCERI);

		for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
			if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
				continue;
			}
			if (_cutnb > 0.0 && (!(*atomJ)->hasCoor() || (*atomI)->groupDistance(**atomJ, stamp) > _cutnb)) {
				// no coordinates or atom further away from cutoff, skip this atom
				continue;
			}

			string atomJtype = (*atomJ)->getType();
			bool special = false;
			if ((*atomI)->isBoundTo(*atomJ) || (*atomI)->isOneThree(*atomJ)) {
				special = true;
			}
			/ *  // I don't think EEF1 discriminates 1-4, what about the vdw radius.  Test against charmm
			if ((*atomI)->isOneFour(*atomJ)) {
				if (!special) {
					// if it is also 1-3 or 1-2 do not add the vdw and elec term
					vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
					CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[3]+vdwParamsJ[3]) * vdwRescalingFactor, sqrt(vdwParamsI[2] * vdwParamsJ[2]) );
					CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant,elec14factor, useRdielectric);
					if (_cutnb > 0.0) {
						// if we are using a cutoff, set the Charmm VDW interaction with
						// the cutoffs for the switching function
						pCVI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
						pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
					} else {
						pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
						pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
					}
					ESet->addInteraction(pCVI);
					ESet->addInteraction(pCEI);
				}
				special = true;
			}
			* /
			if (!special) {
			//	vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
				vector<double> vdwParamsJ;;
				if (!pParReader->vdwParam(vdwParamsJ, atomJtype)) {
					cerr << "WARNING 49319: VDW parameters not found for type " << atomItype << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
					continue;
				}
				vector<double> EEF1ParamsJ;
				if (!pEEF1ParReader->EEF1Param(EEF1ParamsJ, atomJtype, _solvent)) {
					cerr << "WARNING 49387: EEF1 parameters not found for type " << atomJtype << ", solvent " << _solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
					continue;
				}
				CharmmEEF1Interaction *pCEI = new CharmmEEF1Interaction(*(*atomI),*(*atomJ), EEF1ParamsI[0], EEF1ParamsI[2], EEF1ParamsI[3], vdwParamsI[0], EEF1ParamsJ[0], EEF1ParamsJ[2], EEF1ParamsJ[3], vdwParamsJ[0]);
				if (_cutnb > 0.0) {
					// if we are using a cutoff, set the Charmm VDW interaction with
					// the cutoffs for the switching function
					pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctonnb);
				} else {
					pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
				}
				ESet->addInteraction(pCEI);
			}
		}
	}

	return true;
}
*/

void CharmmSystemBuilder::getAtomPointersFromMulti(string _name, vector<Atom*> & _out, vector<CharmmTopologyResidue*> & _position, vector<map<string, Atom*> > & _atomMap) {
	for(vector<CharmmTopologyResidue*>::iterator idItr = _position.begin(); idItr != _position.end(); idItr++ ) {
		map <string, Atom*>::iterator found =  _atomMap[idItr-_position.begin()].find(_name);
		if (found != _atomMap[idItr-_position.begin()].end()) {	
			_out.push_back(_atomMap[idItr-_position.begin()][_name]);
		}	
	}
}

vector<Atom*> CharmmSystemBuilder::getAtomPointers(string _name, vector<vector<vector<CharmmTopologyResidue*> > >::iterator & _chItr, vector<vector<CharmmTopologyResidue*> >::iterator & _posItr, vector<CharmmTopologyResidue*>::iterator & _idItr) {
	vector<Atom*> out;
	if (_name.substr(0,1) == "-" && _posItr > _chItr->begin()) {
		// get the given atom from all identities of the previous position
		string atomName = _name.substr(1,_name.size()-1);
		vector<CharmmTopologyResidue*> & prevPos = *(_posItr-1);
		vector<map<string, Atom*> > & prevPosAtomMap = atomMap[_chItr-polymerDefi.begin()][_posItr-1-_chItr->begin()];
		getAtomPointersFromMulti(atomName, out, prevPos, prevPosAtomMap);
	} else if (_name.substr(0,1) == "+" && _posItr < (_chItr->end() - 1)) {
		// get the given atom from all identities of the next position
		string atomName = _name.substr(1,_name.size()-1);
		vector<CharmmTopologyResidue*> & nextPos = *(_posItr+1);
		vector<map<string, Atom*> > & nextPosAtomMap = atomMap[_chItr-polymerDefi.begin()][_posItr+1-_chItr->begin()];
		getAtomPointersFromMulti(atomName, out, nextPos, nextPosAtomMap);
	} else {
		// get the atom from the current identity of this position
		map <string, Atom *>::iterator found =  atomMap[_chItr-polymerDefi.begin()][_posItr-_chItr->begin()][_idItr-_posItr->begin()].find(_name);						
		if( found != atomMap[_chItr-polymerDefi.begin()][_posItr-_chItr->begin()][_idItr-_posItr->begin()].end() ) {
			out.push_back(atomMap[_chItr - polymerDefi.begin()][_posItr-_chItr->begin()][_idItr-_posItr->begin()][_name]);
		}
	}
	return out;
}
