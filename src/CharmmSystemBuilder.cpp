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

#include "CharmmSystemBuilder.h"

using namespace MSL;
using namespace std;


CharmmSystemBuilder::CharmmSystemBuilder() {
	setup();
}

CharmmSystemBuilder::CharmmSystemBuilder(System & _system, string _topologyFile, string _parameterFile, string _solvationFile) {
	setup();
	pSystem = &_system;
	bool OK = true;
	if (!readTopology(_topologyFile)) {
		OK = false;
	}
	if (!readParameters(_parameterFile)) {
		OK = false;
	}
	if (_solvationFile != "") {
		if (!readSolvation(_solvationFile)) {
			OK = false;
		}
	}
	fail_flag = !OK;
}

CharmmSystemBuilder::CharmmSystemBuilder(const CharmmSystemBuilder & _sysBuild) {
	setup();
	copy(_sysBuild);
}

CharmmSystemBuilder::~CharmmSystemBuilder() {
	delete pTopReader;
	delete pParReader;
	delete pEEF1ParReader;
	deletePointers();
}

void CharmmSystemBuilder::operator=(const CharmmSystemBuilder & _sysBuild) {
	copy(_sysBuild);
}


void CharmmSystemBuilder::setup() {
	setBuildAllTerms();
	fail_flag = false;
	pSystem = NULL;
	pTopReader = new CharmmTopologyReader;
	pParReader = new CharmmParameterReader;
	pEEF1ParReader = new CharmmEEF1ParameterReader;
	vdwRescalingFactor = 1.0;
	buildNonBondedInteractions = true;
	elec14factor = 1;
	dielectricConstant = 1;
	useRdielectric = true;
	useSolvation = false;
	useGroupCutoffs = true;
	halfThickness = 15;
	exponent = 10;
	solvent = pEEF1ParReader->getDefaultSolvent();
}

void CharmmSystemBuilder::copy(const CharmmSystemBuilder & _sysBuild) {
	pSystem = _sysBuild.pSystem;
	*pTopReader = *_sysBuild.pTopReader;
	*pParReader = *_sysBuild.pParReader;
	*pEEF1ParReader = *_sysBuild.pEEF1ParReader;
	buildNonBondedInteractions = _sysBuild.buildNonBondedInteractions;
	vdwRescalingFactor = _sysBuild.vdwRescalingFactor;
	elec14factor = _sysBuild.elec14factor;
	dielectricConstant = _sysBuild.dielectricConstant;
	useRdielectric = _sysBuild.useRdielectric;
	termsToBuild = _sysBuild.termsToBuild;
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

bool CharmmSystemBuilder::addIdentity(string _positionId, string _resName, string _bbAtoms) {
	//vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_positionId, vector<string>(1, _resName),_bbAtoms);
}

bool CharmmSystemBuilder::addIdentity(Position & _pos, string _resName, string _bbAtoms) {
	//vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_pos, vector<string>(1, _resName), _bbAtoms);
}

bool CharmmSystemBuilder::addIdentity(string _positionId, const vector<string> & _resNames, string _bbAtoms) {
	//vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	if (pSystem->positionExists(_positionId)) {
		Position & pos = pSystem->getLastFoundPosition();
		return addIdentity(pos, _resNames, _bbAtoms);
	} else {
		cerr << "WARNING 32978: position " << _positionId << " not found in System at bool CharmmSystemBuilder::addIdentity(string _positionId, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}
	//return addIdentity(_positionId, _resNames, _bbAtoms);
}

bool CharmmSystemBuilder::removeIdentity(std::string _positionId, string _resName) {
	if (pSystem->positionExists(_positionId)) {
		Position & pos = pSystem->getLastFoundPosition();
		return removeIdentity(pos, _resName);
	} else {
		cerr << "WARNING 32983: position " << _positionId << " not found in System at bool CharmmSystemBuilder::removeIdentity(std::string _positionId, string _resName, string _bbAtoms)" << endl;
		return false;
	}
}
bool CharmmSystemBuilder::removeIdentity(Position & _pos, string _resName) {

	if (!_pos.identityExists(_resName)) {
		cerr << "WARNING 32981: position " << _pos << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}
	Residue * res = &_pos.getLastFoundIdentity();
	AtomPointerVector resAtoms = res->getAtomPointers();

	/**********************************************************
	 * FIND THE INDECES OF CHAIN AND POSITION IN CHAIN
	 *
	 * TODO: once the getPositionIndex from Chain and System
	 *       are implemented convert to those for finding the
	 *       position index
	 *
	 **********************************************************/
	// first add the atoms to the position and update the atomMap
	bool posIndexFound = false;
	Chain * pParentChain = _pos.getParentChain();
	unsigned int posIndex = 0;
	for (unsigned int i=0; i<pParentChain->positionSize(); i++) {
		Position * chainPos = &(pParentChain->getPosition(i));
		if (&_pos == chainPos) {
			posIndex = i;
			posIndexFound = true;
			break;
		}
	}
	if (!posIndexFound) {
		cerr << "WARNING 32983: position " << _pos << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}
	bool chainIndexFound = false;
	System * pParentSys = _pos.getParentSystem();
	unsigned int chainIndex = 0;
	for (unsigned int i=0; i<pParentSys->chainSize(); i++) {
		Chain * sysChain = &(pParentSys->getChain(i));
		if (pParentChain == sysChain) {
			chainIndex = i;
			chainIndexFound = true;
			break;
		}
	}
	if (!chainIndexFound) {
		cerr << "WARNING 32983: chain " << *pParentChain << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}

	// find the identity in the polymerDefi and erase it
	bool found = false;
	unsigned int idIndex = 0;
	for (unsigned int i=0; i<polymerDefi[chainIndex][posIndex].size(); i++) {
		if (polymerDefi[chainIndex][posIndex][i]->getName() == _resName) {
			found = true;
			idIndex = i;
			delete polymerDefi[chainIndex][posIndex][i];
			polymerDefi[chainIndex][posIndex].erase(polymerDefi[chainIndex][posIndex].begin() + i);
			break;
		}
	}
	if (!found) {
		cerr << "WARNING 32988: identity " << _pos.getResidueName() << " not found in position in bool CharmmSystemBuilder::removeIdentity(Position & _pos, string _resName, string _bbAtoms)" << endl;
		return false;
	}
	// remove the atoms from the atom map
	atomMap[chainIndex][posIndex].erase(atomMap[chainIndex][posIndex].begin()+idIndex);

	// remove the atoms from the energy set
	EnergySet* ESet = _pos.getParentSystem()->getEnergySet();

	ESet->deleteInteractionsWithAtoms(resAtoms);
	
	// remove the identity from the position
	_pos.removeIdentity(_resName, true);

	// remove all the IcEntries that now refer to atoms that are gone
	pParentSys->purgeIcTable(); 

	return true;
}

/*
HERE!!!
bool CharmmSystemBuilder::mutate(string _positionId, string _newResname, string _bbAtoms="N CA C O HN") {

}

bool CharmmSystemBuilder::mutate(Position & _pos, , string _newResname, string _bbAtoms="N CA C O HN") {
	if (pSystem->positionExists(_positionId)) {
		Position & pos = pSystem->getLastFoundPosition();
		return addIdentity(pos, _resNames, _bbAtoms);
	} else {
		cerr << "WARNING 32978: position " << _positionId << " not found in System at bool CharmmSystemBuilder::addIdentity(string _positionId, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}
}

bool CharmmSystemBuilder::mutate(string _positionId, string _oldResName, string _newResname, string _bbAtoms="N CA C O HN") {
}

bool CharmmSystemBuilder::mutate(Position & _pos, , string _oldResName, string _newResname, string _bbAtoms="N CA C O HN") {
}
*/

bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms) {

	vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	/**********************************************************
	 * FIND THE INDECES OF CHAIN AND POSITION IN CHAIN
	 *
	 * TODO: once the getPositionIndex from Chain and System
	 *       are implemented convert to those for finding the
	 *       position index
	 *
	 **********************************************************/
	// first add the atoms to the position and update the atomMap
	bool posIndexFound = false;
	Chain * pParentChain = _pos.getParentChain();
	unsigned int posIndex = 0;
	for (unsigned int i=0; i<pParentChain->positionSize(); i++) {
		Position * chainPos = &(pParentChain->getPosition(i));
		if (&_pos == chainPos) {
			posIndex = i;
			posIndexFound = true;
			break;
		}
	}
	if (!posIndexFound) {
		cerr << "WARNING 32983: position " << _pos << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}
	bool chainIndexFound = false;
	System * pParentSys = _pos.getParentSystem();
	unsigned int chainIndex = 0;
	for (unsigned int i=0; i<pParentSys->chainSize(); i++) {
		Chain * sysChain = &(pParentSys->getChain(i));
		if (pParentChain == sysChain) {
			chainIndex = i;
			chainIndexFound = true;
			break;
		}
	}
	if (!chainIndexFound) {
		cerr << "WARNING 32983: chain " << *pParentChain << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
		return false;
	}

	/*********************************************************
	 *  These variables will be used to cycle in the position where we
	 *  will add the identities
	 *********************************************************/
	unsigned int idIndexStart = polymerDefi[chainIndex][posIndex].size();
	unsigned int idIndexEnd = idIndexStart + _resNames.size() - 1;

	Residue & currentRes = _pos.getCurrentIdentity();
	map<string, bool> bbAtomsMap;
	//for (vector<string>::iterator k=_bbAtoms.begin(); k!=_bbAtoms.end(); k++) {
	for (vector<string>::iterator k=bbAtoms.begin(); k!=bbAtoms.end(); k++) {
		bbAtomsMap[*k] = true;
	}

	/*********************************************************
	 *  Create an atom vector with all the atoms of the new
	 *  identities and add them to the position (the Position
	 *  object will take care of creating separate residues
	 *********************************************************/
	vector<string> processedResNames; // after removing a possible patch
	for (vector<string>::const_iterator k=_resNames.begin(); k!=_resNames.end(); k++) {
		// Extend the polymerDefi
		vector<string> split = MslTools::tokenize(*k, "-");
		processedResNames.push_back(split[0]);
		if (!pTopReader->residueExists(split[0])) {
			cerr << "WARNING 32988: residue " << split[0] << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
			return false;
		}
		polymerDefi[chainIndex][posIndex].push_back(new CharmmTopologyResidue(pTopReader->getLastFoundResidue()));
		atomMap[chainIndex][posIndex].push_back(map<string, Atom*>());

		// if the patches for the first or the last residue were not declared
		// get them from the residues
		if (split.size() == 1 && posIndex == 0) {
			// add default first residue patch to current identity
			split.push_back(polymerDefi[chainIndex][posIndex].back()->getFirstDefaultPatch());
		}
		if (split.size() == 1 && posIndex == polymerDefi[chainIndex].size()-1) {
			// add default last residue patch to current identity
			split.push_back(polymerDefi[chainIndex][posIndex].back()->getLastDefaultPatch());
		}
		/************************************************************
		 *  APPLY THE PATCHES (if any)
		 ************************************************************/
		for (vector<string>::iterator n=split.begin()+1; n!=split.end(); n++) {
			// patch current identity
			if (pTopReader->residueExists(*n)) {
				if (!polymerDefi[chainIndex][posIndex].back()->applyPatch(pTopReader->getLastFoundResidue())) {
					cerr << "WARNING 19134: cannot apply patch " << *n << " to residue, in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)";
					return false;
				}
			}
		}

		// Add the atoms to the position
		string chainId = _pos.getChainId();
		int resNum = _pos.getResidueNumber();
		string iCode = _pos.getResidueIcode();
		AtomPointerVector posAtoms;
		for (unsigned int l=0; l<polymerDefi[chainIndex][posIndex].back()->atomSize(); l++) {
			string name;
			string type;
			double partialCharge;
			string element;
			int group;
			polymerDefi[chainIndex][posIndex].back()->getTopolAtom(l, name, type, partialCharge, element, group);
			Atom * a = new Atom(name);
			a->setChainId(chainId);
			a->setResidueName(*k);
			a->setResidueNumber(resNum);
			a->setResidueIcode(iCode);
			a->setCharge(partialCharge);
			a->setType(type);
			a->setElement(element);
			a->setGroupNumber(group);
			a->wipeCoordinates();
			if (bbAtomsMap.find(name) != bbAtomsMap.end()) {
				// it is a backbone atom, let's copy the coordinates from the current identity
				if (currentRes.atomExists(name)) {
					a->setCoor(currentRes.getLastFoundAtom().getCoor());
				}
			}
			posAtoms.push_back(a);
		}

		_pos.addAtoms(posAtoms);
		// populate the new identity in the position
		// garbage collection
		for (AtomPointerVector::iterator k=posAtoms.begin(); k!=posAtoms.end(); k++) {
			delete *k;
		}
		posAtoms.clear();
	}


	/************************************************************
	 *  Update the atoms map
	 * 
	 *  Also, create a lookup table (newAtomsLookupMap) with all 
	 *  the atoms (by their address) that are new to the position 
	 *  (used to identify if a new IC entry or bonded term includes 
	 *  a new atom or not.
	 ************************************************************/
	map<Atom*, bool> newAtomsLookupMap;
	for (unsigned int i=idIndexStart; i<=idIndexEnd; i++) {
		// for all the new identities
		for (unsigned int j=0; j<polymerDefi[chainIndex][posIndex][i]->atomSize(); j++) {
			// find the atom created in the Position and add it to the map
			string name;
			string type;
			double partialCharge;
			string element;
			int group;
			polymerDefi[chainIndex][posIndex][i]->getTopolAtom(j, name, type, partialCharge, element, group);
			if (_pos.atomExists(processedResNames[i-idIndexStart], name)) {
				atomMap[chainIndex][posIndex][i][name] = &(_pos.getLastFoundAtom());
				newAtomsLookupMap[ &(_pos.getLastFoundAtom()) ] = true;
			} else {
				cerr << "WARNING 32998: atom " << processedResNames[i-idIndexStart] << "," << name << " not found in topology at bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
				return false;
			}
		}
	}

	/************************************************************
	 *  Add the IC table
	 * 
	 *  NOTE: for the IC table and the bonded interactions we need to cycle
	 *        with the previous and next position to get any references to atoms
	 *        in the added residue (provided as bond to +N in the previous, for example)
	 ************************************************************/
	vector<vector<vector<CharmmTopologyResidue*> > >::iterator chItr=polymerDefi.begin()+chainIndex;
	vector<vector<CharmmTopologyResidue*> >::iterator posThis=chItr->begin()+posIndex;
	vector<vector<CharmmTopologyResidue*> >::iterator posStart=chItr->begin()+posIndex;
	bool checkPrevious = false;
	if (posIndex > 0) {
		checkPrevious = true;
		posStart=chItr->begin()+posIndex-1;
	}
	vector<vector<CharmmTopologyResidue*> >::iterator posEnd=chItr->begin()+posIndex;
	bool checkNext = false;
	if (posIndex < polymerDefi[chainIndex].size()-1) {
		checkNext = true;
		posEnd=chItr->begin()+posIndex+1;
	}
	//vector<CharmmTopologyResidue*>::iterator idThig=posThis->end()-1;

	vector<string> icAtoms;
	vector<double> icValues;
	bool improperFlag;
	for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=posStart; posItr<=posEnd; posItr++) {
		// for the postion before, this position and the next one...
		unsigned int start = 0;
		if (posItr == posThis) {
			// if this is the position where we added the atoms
			// cycles only from the new identities to the end
			start = idIndexStart;
		}
		for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin()+start; idItr!=posItr->end(); idItr++) {
			// for each identity
			for (unsigned int i=0; i<(*idItr)->icSize(); i++) {
				(*idItr)->getIcLine(i, icAtoms, icValues, improperFlag);
				if (checkPrevious && posItr == posThis-1 && icAtoms[0] != "+" && icAtoms[1] != "+" && icAtoms[2] != "+" && icAtoms[3] != "+") {
					// this is the previous position but atoms have reference to the added residue
					continue;
				} else if (checkNext && posItr == posThis+1 && icAtoms[0] != "-" && icAtoms[1] != "-" && icAtoms[2] != "-" && icAtoms[3] != "-") {
					// this is the next position but atoms have reference to the added residue
					continue;
				}
				vector<vector<Atom*> > icAtomPointers(4, vector<Atom*>());
				for (vector<string>::iterator n=icAtoms.begin(); n!=icAtoms.end(); n++) {
					if (n->substr(0,1) == "-") {
						if (posItr > chItr->begin()) {
							// not the first residue of the chain
							vector<vector<CharmmTopologyResidue*> >::iterator ll=posItr-1;
							for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
								// for each identity get all the -X atoms

								map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
								if (found == atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].end()) {
									icAtomPointers[n-icAtoms.begin()].push_back(NULL);
								} else {
									icAtomPointers[n-icAtoms.begin()].push_back(found->second);
								}

							}
						} else {
							icAtomPointers[n-icAtoms.begin()].push_back(NULL);
						}
					} else if (n->substr(0,1) == "+") {
						if (posItr<chItr->end()-1) {
							// not the last residue of the chain
							vector<vector<CharmmTopologyResidue*> >::iterator ll=posItr+1;
							for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
								// for each identity get all the -X atoms

								map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
								if (found == atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].end()) {
									icAtomPointers[n-icAtoms.begin()].push_back(NULL);
								} else {
									icAtomPointers[n-icAtoms.begin()].push_back(found->second);
								}

							}
						} else {
							icAtomPointers[n-icAtoms.begin()].push_back(NULL);
						}
					} else {
						map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][posItr-chItr->begin()][idItr-posItr->begin()].find(*n);
						if (found == atomMap[chItr-polymerDefi.begin()][posItr-chItr->begin()][idItr-posItr->begin()].end()) {
							icAtomPointers[n-icAtoms.begin()].push_back(NULL);
						} else {
							icAtomPointers[n-icAtoms.begin()].push_back(found->second);
						}
					}
				}
				// add the IC terms to the system
				for (vector<Atom*>::iterator atmK=icAtomPointers[0].begin(); atmK!=icAtomPointers[0].end(); atmK++) {
					for (vector<Atom*>::iterator atmL=icAtomPointers[1].begin(); atmL!=icAtomPointers[1].end(); atmL++) {
						for (vector<Atom*>::iterator atmM=icAtomPointers[2].begin(); atmM!=icAtomPointers[2].end(); atmM++) {
							for (vector<Atom*>::iterator atmN=icAtomPointers[3].begin(); atmN!=icAtomPointers[3].end(); atmN++) {
								if ((*atmK != NULL || *atmN != NULL) && *atmL != NULL && *atmM != NULL) {
									
									bool foundIdentityAtom = false;
									if ( (*atmK != NULL && newAtomsLookupMap.find(*atmK) != newAtomsLookupMap.end()) || (*atmL != NULL && newAtomsLookupMap.find(*atmL) != newAtomsLookupMap.end()) || (*atmM != NULL && newAtomsLookupMap.find(*atmM) != newAtomsLookupMap.end()) || (*atmN != NULL && newAtomsLookupMap.find(*atmN) != newAtomsLookupMap.end())) {
										// verify that at least one atom is a new atom
										foundIdentityAtom = true;
									}
									if (!foundIdentityAtom) {
										// no atom is a new atom, do not add this IC
										continue;
									}

									// If icValues are zero (usually zero for the PatchResidues' atoms, get the minimum values from the parameter file
									if (icValues[0] == 0.0) {
										if(improperFlag) {
											if(*atmK) {
												vector<double> param;
												if (pParReader->bondParam(param, (*atmK)->getType(),(*atmM)->getType())) {
													icValues[0] = param[1];
												}
											}
										} else {
											if(*atmK) {
												vector<double> param;
												if (pParReader->bondParam(param, (*atmK)->getType(),(*atmL)->getType())) {
													icValues[0] = param[1];
												}
											}
										}
									}	

									if (icValues[1] == 0.0) {
										if(improperFlag) {
											if(*atmK) {
												vector<double> param;
												if (pParReader->angleParam(param, (*atmK)->getType(),(*atmM)->getType(),(*atmL)->getType())) {
													icValues[1] = param[1];
												}
											}
										} else {
											if(*atmK) {
												vector<double> param;
												if (pParReader->angleParam(param, (*atmK)->getType(),(*atmL)->getType(),(*atmM)->getType())) {
													icValues[1] = param[1];
												}
											}
										}
									}
									
									if (icValues[3] == 0.0) {
										if(*atmN) {
											vector<double> param;
											if (pParReader->angleParam(param, (*atmL)->getType(),(*atmM)->getType(),(*atmN)->getType())) {
												icValues[3] = param[1];
											}
										}
									}
									
									if (icValues[4] == 0.0) {
										if(*atmN) {
											vector<double> param;
											if (pParReader->bondParam(param, (*atmM)->getType(),(*atmN)->getType())) {
												icValues[4] = param[1];
											}
										}
									}

									_pos.getParentSystem()->addIcEntry(*atmK, *atmL, *atmM, *atmN, icValues[0], icValues[1], icValues[2], icValues[3], icValues[4], improperFlag);
								}
							}
						}
					}
				}


			}
		}
	}
	
	/*********************************************************************************
	 *
	 *  Attempt to build the atoms from IC
	 * 
	 **********************************************************************************/
	for (unsigned int i=0; i<_resNames.size(); i++) {
		Residue * pRes = &_pos.getResidue(_resNames[i]);
		AtomPointerVector atoms = pRes->getAtomPointers();
		for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
			(*k)->buildFromIc(false); // build only from active atoms = false
		}
	}


	/*********************************************************************************
	 *
	 *  ADD THE INTERACTION TERMS:
	 *
	 *  For each pair of atoms, add the appropriate interactions:
	 *   - bond if 1-2
	 *   - angle if 1-3 and autogenerate is on, or as declared
	 *   - dihedral if 1-4 and autogenerate is on, or as declared
	 *   - improper dihedrals, as declared
	 **********************************************************************************/
	EnergySet* ESet = _pos.getParentSystem()->getEnergySet();
        

	// Loop over each position.
	for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=posStart; posItr<=posEnd; posItr++) {
		// for the postion before, this position and the next one...
		unsigned int start = 0;
		if (posItr == posThis) {
			// if this is the position where we added the atoms
			// cycles only from the new identities to the end
			start = idIndexStart;
		}
		for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin()+start; idItr!=posItr->end(); idItr++) {
			// for each identity

			// add the bonds
			for (unsigned int l=0; l<(*idItr)->bondSize(); l++) {
				// get the atoms names
				string atom1,atom2;
				unsigned int type;
				(*idItr)->getBond(l, atom1, atom2, type);           
				if (checkPrevious && posItr == posThis-1 && atom1.substr(0,1) != "+" && atom2.substr(0,1) != "+") {
					// this is the previous position but atoms have reference to the added residue
					continue;
				} else if (checkNext && posItr == posThis+1 && atom1.substr(0,1) != "-" && atom2.substr(0,1) != "-") {
					// this is the next position but atoms have reference to the added residue
					continue;
				}
			
				/***********************************************************************
				 *   If multiple identities are present, and some atoms refers to previous 
				 *   (i.e. -N) or next (+N) residue, then the bonds should be added in
				 *   all combinations
				 ***********************************************************************/
			
				vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr);
				vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
				
				bool foundIdentityAtom = false;
				// Now loop over all the pAtom1 and pAtom2 atoms and add an interaction for each combination
				for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
					if (newAtomsLookupMap.find(*a1) != newAtomsLookupMap.end()) {
						// found a new atom in this bonded term
						foundIdentityAtom = true;
					}
					string type1 = (*a1)->getType();
					for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
						if (!foundIdentityAtom && newAtomsLookupMap.find(*a2) == newAtomsLookupMap.end()) {
							// no atom is a new atom, do not add this bonded term
							continue;
						}
						string type2 = (*a2)->getType();
						(*a1)->setBoundTo(*a2);
						vector<double> params;
						if (pParReader->bondParam(params, type1, type2)) {
							if (termsToBuild["CHARMM_BOND"]) {
								CharmmBondInteraction *pCBI = new CharmmBondInteraction(*(*a1),*(*a2),params[0],params[1]);
								ESet->addInteraction(pCBI);
							}
						} else {
							cerr << "WARNING: CHARMM bond parameters not found for atom types " << type1 << ", " << type2 << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
						}
					}
				}
					
			}
		}
	}
	/************************* DONE SETTING UP THE BONDS **********************************/

	// add the angles if they are not autogenerated
	if (!pTopReader->getAutoGenerateAngles()) {
		// Loop over each position.
		for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=posStart; posItr<=posEnd; posItr++) {
			// for the postion before, this position and the next one...
			unsigned int start = 0;
			if (posItr == posThis) {
				// if this is the position where we added the atoms
				// cycles only from the new identities to the end
				start = idIndexStart;
			}
			for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin()+start; idItr!=posItr->end(); idItr++) {

				for (unsigned int l=0; l<(*idItr)->angleSize(); l++) {
					// Loop over each improper in the residue.
					
					string atom1,atom2,atom3;
					// get the atoms names
					(*idItr)->getAngle(l, atom1, atom2, atom3);
					if (checkPrevious && posItr == posThis-1 && atom1.substr(0,1) != "+" && atom2.substr(0,1) != "+" && atom3.substr(0,1) != "+") {
						// this is the previous position but atoms have reference to the added residue
						continue;
					} else if (checkNext && posItr == posThis+1 && atom1.substr(0,1) != "-" && atom2.substr(0,1) != "-" && atom3.substr(0,1) != "-") {
						// this is the next position but atoms have reference to the added residue
						continue;
					}
				
					vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
					vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
					vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
					
					bool foundIdentityAtom = false;
					for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
						if (newAtomsLookupMap.find(*a1) != newAtomsLookupMap.end()) {
							// found a new atom in this bonded term
							foundIdentityAtom = true;
						}
						string type1 = (*a1)->getType();
						for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
							if (!foundIdentityAtom && newAtomsLookupMap.find(*a2) != newAtomsLookupMap.end()) {
								// found a new atom in this bonded term
								foundIdentityAtom = true;
							}
							string type2 = (*a2)->getType();
							for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
								if (!foundIdentityAtom && newAtomsLookupMap.find(*a3) == newAtomsLookupMap.end()) {
									// no atom is a new atom, do not add this bonded term
									continue;
								}
								string type3 = (*a3)->getType();
								//vector<double> params = pParReader->ureyBradleyParam(type1, type2, type3);
								vector<double> params;
								if (pParReader->ureyBradleyParam(params, type1, type2, type3)) {
									if (termsToBuild["CHARMM_U-BR"]) {
										CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*a1),*(*a3),params[0],params[1]);
										ESet->addInteraction(pCUI);
									}
								}
								if (pParReader->angleParam(params, type1, type2, type3)) {
									if (termsToBuild["CHARMM_ANGL"]) {
										CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*a1),*(*a2),*(*a3),params[0],params[1] * M_PI / 180.0); 
										ESet->addInteraction(pCAI);
									}
								} else {
									cerr << "WARNING: CHARMM angle parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
								}
							}
						}
					}
						
				}
			}
		}
	}
	/************************* DONE SETTING UP THE ANGLES **********************************/

	// add the dihedrals if they are not autogenerated
	if (!pTopReader->getAutoGenerateDihedrals()) {
		for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=posStart; posItr<=posEnd; posItr++) {
			// for the postion before, this position and the next one...
			unsigned int start = 0;
			if (posItr == posThis) {
				// if this is the position where we added the atoms
				// cycles only from the new identities to the end
				start = idIndexStart;
			}
			for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin()+start; idItr!=posItr->end(); idItr++) {


				for (unsigned int l=0; l<(*idItr)->dihedralSize(); l++) {
					// Loop over each dihedral in the residue.
					
					string atom1,atom2,atom3,atom4;
					// get the atoms names
					(*idItr)->getDihedral(l, atom1, atom2, atom3, atom4);
					if (checkPrevious && posItr == posThis-1 && atom1.substr(0,1) != "+" && atom2.substr(0,1) != "+" && atom3.substr(0,1) != "+" && atom4.substr(0,1) != "+") {
						// this is the previous position but atoms have reference to the added residue
						continue;
					} else if (checkNext && posItr == posThis+1 && atom1.substr(0,1) != "-" && atom2.substr(0,1) != "-" && atom3.substr(0,1) != "-" && atom4.substr(0,1) != "-") {
						// this is the next position but atoms have reference to the added residue
						continue;
					}
				
					vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
					vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
					vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
					vector<Atom*> pAtom4 = getAtomPointers(atom4, chItr, posItr, idItr);
					
					bool foundIdentityAtom = false;
					for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
						if (newAtomsLookupMap.find(*a1) != newAtomsLookupMap.end()) {
							// found a new atom in this bonded term
							foundIdentityAtom = true;
						}
						string type1 = (*a1)->getType();
						for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
							if (!foundIdentityAtom && newAtomsLookupMap.find(*a2) != newAtomsLookupMap.end()) {
								// found a new atom in this bonded term
								foundIdentityAtom = true;
							}
							string type2 = (*a2)->getType();
							for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
								if (!foundIdentityAtom && newAtomsLookupMap.find(*a3) != newAtomsLookupMap.end()) {
									// found a new atom in this bonded term
									foundIdentityAtom = true;
								}
								string type3 = (*a3)->getType();
								for(vector<Atom*>::iterator a4 = pAtom4.begin() ; a4 != pAtom4.end(); a4++) {
									if (!foundIdentityAtom && newAtomsLookupMap.find(*a4) == newAtomsLookupMap.end()) {
										// no atom is a new atom, do not add this bonded term
										continue;
									}
									string type4 = (*a4)->getType();

									vector<vector <double> > dihedralEntries;
									if (pParReader->dihedralParam(dihedralEntries, type1, type2, type3, type4)) {
										// there could be multiple entries for a single dihedral
										for(int m = 0; m < dihedralEntries.size() ; m++) {
											// the delta should be expressed in radians
											dihedralEntries[m][2] = dihedralEntries[m][2]*M_PI/180.0;

										}
										if (termsToBuild["CHARMM_DIHE"]) {
											CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*a1),*(*a2),*(*a3),*(*a4),dihedralEntries);
											ESet->addInteraction(pCDI);
										}
									} else {
										cerr << "WARNING: CHARMM dihedral parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << ", " << type4 << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
									}
								
								}
							}
						}
					}
						
				}
			}
		}
	}
	/************************* DONE SETTING UP THE DIHEDRALS **********************************/

	// add the impropers
	for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=posStart; posItr<=posEnd; posItr++) {
		// for the postion before, this position and the next one...
		unsigned int start = 0;
		if (posItr == posThis) {
			// if this is the position where we added the atoms
			// cycles only from the new identities to the end
			start = idIndexStart;
		}
		for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin()+start; idItr!=posItr->end(); idItr++) {
			for (unsigned int l=0; l<(*idItr)->improperSize(); l++) {
				// Loop over each improper in the residue.
				
				string atom1,atom2,atom3,atom4;
				// get the atoms names
				(*idItr)->getImproper(l, atom1, atom2, atom3, atom4);
				if (checkPrevious && posItr == posThis-1 && atom1.substr(0,1) != "+" && atom2.substr(0,1) != "+" && atom3.substr(0,1) != "+" && atom4.substr(0,1) != "+") {
					// this is the previous position but atoms have reference to the added residue
					continue;
				} else if (checkNext && posItr == posThis+1 && atom1.substr(0,1) != "-" && atom2.substr(0,1) != "-" && atom3.substr(0,1) != "-" && atom4.substr(0,1) != "-") {
					// this is the next position but atoms have reference to the added residue
					continue;
				}
			
				vector<Atom*> pAtom1 = getAtomPointers(atom1, chItr, posItr, idItr); 
				vector<Atom*> pAtom2 = getAtomPointers(atom2, chItr, posItr, idItr);
				vector<Atom*> pAtom3 = getAtomPointers(atom3, chItr, posItr, idItr);
				vector<Atom*> pAtom4 = getAtomPointers(atom4, chItr, posItr, idItr);
				
				bool foundIdentityAtom = false;
				for(vector<Atom*>::iterator a1 = pAtom1.begin() ; a1 != pAtom1.end(); a1++) {
					if (newAtomsLookupMap.find(*a1) != newAtomsLookupMap.end()) {
						// found a new atom in this bonded term
						foundIdentityAtom = true;
					}
					string type1 = (*a1)->getType();
					for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
						if (!foundIdentityAtom && newAtomsLookupMap.find(*a2) != newAtomsLookupMap.end()) {
							// found a new atom in this bonded term
							foundIdentityAtom = true;
						}
						string type2 = (*a2)->getType();
						for(vector<Atom*>::iterator a3 = pAtom3.begin() ; a3 != pAtom3.end(); a3++) {
							if (!foundIdentityAtom && newAtomsLookupMap.find(*a3) != newAtomsLookupMap.end()) {
								// found a new atom in this bonded term
								foundIdentityAtom = true;
							}
							string type3 = (*a3)->getType();
							for(vector<Atom*>::iterator a4 = pAtom4.begin() ; a4 != pAtom4.end(); a4++) {
								if (!foundIdentityAtom && newAtomsLookupMap.find(*a4) == newAtomsLookupMap.end()) {
									// no atom is a new atom, do not add this bonded term
									continue;
								}
								string type4 = (*a4)->getType();
								//vector<double> improperParams = pParReader->improperParam(type1, type2, type3, type4);	
								vector<double> improperParams;
								if (pParReader->improperParam(improperParams, type1, type2, type3, type4)) {
									if (termsToBuild["CHARMM_IMPR"]) {
										CharmmImproperInteraction *pCII = new CharmmImproperInteraction(*(*a1),*(*a2),*(*a3),*(*a4),improperParams[0],improperParams[1]*M_PI/180.0);
										ESet->addInteraction(pCII);
									}
								} else {
									cerr << "WARNING: CHARMM improper parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << ", " << type4 << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
								}
							}
						}
					}
				}
			}
		}
	}
	/************************* DONE SETTING UP THE IMPROPERS **********************************/


	// add autogenerated angles and dihedral if nececssary (NOTE, this could be coded better:
	// essentially it is cycling around all atoms and adding interactions if at least one atom
	// belongs to the newly added residue
	if (pTopReader->getAutoGenerateAngles()	|| pTopReader->getAutoGenerateDihedrals()) {
		AtomPointerVector atoms = _pos.getParentSystem()->getAllAtomPointers();

		for(AtomPointerVector::iterator atomI = atoms.begin(); atomI < atoms.end(); atomI++) {
			bool foundIdentityAtom = false;
			if (newAtomsLookupMap.find(*atomI) != newAtomsLookupMap.end()) {
				foundIdentityAtom = true;
			}

			string atomItype = (*atomI)->getType();

			for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
				if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
					continue;
				}
				if (!foundIdentityAtom && newAtomsLookupMap.find(*atomJ) != newAtomsLookupMap.end()) {
					foundIdentityAtom = true;
				}
				string atomJtype = (*atomJ)->getType();
				if (pTopReader->getAutoGenerateAngles() && (*atomI)->isOneThree(*atomJ)) {
					/*******************************************************
					 *  Add the bond term only if at least one of the 3 atoms
					 *  belongs to the newly added residue
					 *******************************************************/

					// autogenerate the angle interactions
					vector<Atom*> middle = (*atomI)->getOneThreeMiddleAtoms(*atomJ);
					for (vector<Atom*>::iterator k=middle.begin(); k!=middle.end(); k++) {
						if (foundIdentityAtom || newAtomsLookupMap.find(*k) != newAtomsLookupMap.end()) {
							// this angle includes at least one atom in the new residue
							string middleType = (*k)->getType();
							vector<double> params;
							if (pParReader->ureyBradleyParam(params, atomItype, middleType, atomJtype)) {
								if (termsToBuild["CHARMM_U-BR"]) {
									CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*atomI),*(*atomJ),params[0],params[1]);
									ESet->addInteraction(pCUI);
								}
							}
							if (pParReader->angleParam(params, atomItype, middleType, atomJtype)) {
								if (termsToBuild["CHARMM_ANGL"]) {
									CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*atomI),*(*k),*(*atomJ),params[0],params[1] * M_PI / 180.0); 
									ESet->addInteraction(pCAI);
								}
							} else {
								cerr << "WARNING: CHARMM angle parameters not found for atom types " << atomItype << ", " << middleType << ", " << atomJtype << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
							}
						}
					}
				}
				if (pTopReader->getAutoGenerateDihedrals() && (*atomI)->isOneFour(*atomJ)) {
					// autogenerate the dihedral interactions
					vector<vector<Atom*> > middle = (*atomI)->getOneFourMiddleAtoms(*atomJ);
					for (vector<vector<Atom*> >::iterator k=middle.begin(); k!=middle.end(); k++) {
						vector<Atom*> middleAtoms;
						vector<string> middleTypes;
						for (vector<Atom*>::iterator l=k->begin(); l!=k->end(); l++) {
							middleAtoms.push_back(*l);
							middleTypes.push_back((*l)->getType());
						}
						if (foundIdentityAtom || newAtomsLookupMap.find(middleAtoms[0]) != newAtomsLookupMap.end() || newAtomsLookupMap.find(middleAtoms[1]) != newAtomsLookupMap.end()) {
							// this dihedral includes at least one atom in the new residue
							vector<vector <double> > dihedralEntries;
							if (pParReader->dihedralParam(dihedralEntries, atomItype, middleTypes[0], middleTypes[1], atomJtype)) {
								// there could be multiple entries for a single dihedral
								for(int m = 0; m < dihedralEntries.size() ; m++) {
									// the delta should be expressed in radians
									dihedralEntries[m][2] = dihedralEntries[m][2]*M_PI/180.0;

								}
							
								if (termsToBuild["CHARMM_DIHE"]) {
									CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*atomI),*(middleAtoms[0]),*(middleAtoms[1]),*(*atomJ),dihedralEntries);
									ESet->addInteraction(pCDI);
								}
							} else {
								cerr << "WARNING: CHARMM dihedral parameters not found for atom types " << atomItype << ", " << middleTypes[0] << ", " << middleTypes[1] << ", " << atomJtype << " in bool CharmmSystemBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms)" << endl;
							}
						}
					}
				}
			}
		}
	}

	return true;
}

bool CharmmSystemBuilder::buildSystem(const PolymerSequence & _sequence) {

	if (pSystem == NULL) {
		cerr << "WARNING, undefined System in bool CharmmSystemBuilder::buildSystem(const PolymerSequence & _sequence), set System first" << endl;
		return false;
	}

	deletePointers();
	pSystem->reset();

	vector<vector<vector<string> > > seq = _sequence.getSequence();

	
	/************************************************************
	 *  Populate polymerDefi according to the sequence.  polymerDefi
	 *  is a 3D vector vector<vector<vector<CharmmTopologyResidue*> > >
	 *
	 *  The first level corresponds to the chain, the second to
	 *  the position, the third to the identity
	 ************************************************************/
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

				/************************************************************
				 *  ADD THE TOPOLOGY RESIDUE
				 ************************************************************/
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
				int splitSize = split.size();
				// get them from the residues
				if (splitSize == 1 && l==k->begin()) {
					// add default first residue patch to current identity
					split.push_back(polymerDefi.back().back().back()->getFirstDefaultPatch());
				}
				if (splitSize == 1 && l==k->end() - 1) {
					// add default last residue patch to current identity
					split.push_back(polymerDefi.back().back().back()->getLastDefaultPatch());
				}
				/************************************************************
				 *  APPLY THE PATCHES (if any)
				 ************************************************************/
				for (vector<string>::iterator n=split.begin()+1; n!=split.end(); n++) {
					// patch current identity
					if (pTopReader->residueExists(*n)) {
						if (!polymerDefi.back().back().back()->applyPatch(pTopReader->getLastFoundResidue())) {
							cerr << "WARNING 19134: cannot apply patch " << *n << " to residue, in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
							return false;
						}
					} else {
						cerr << "WARNING 19139: cannot apply unexistent patch " << *n << " to residue, in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
						return false;
					}
				}

			}
		}
	}

	/************************************************************
	 *  From the polymerDefi create all the atoms (sysAtoms).
	 *  They will be added to the system that will sort them
	 *  into the appropriate chains, positions and identities
	 ************************************************************/
	AtomPointerVector sysAtoms;
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator i=polymerDefi.begin(); i!=polymerDefi.end(); i++) {
		string chainId = _sequence.getChainId(i-polymerDefi.begin());
		for (vector<vector<CharmmTopologyResidue*> >::iterator j=i->begin(); j!=i->end(); j++) {
			string resNumStr = _sequence.getResidueNumber(i-polymerDefi.begin(), j-i->begin());
			int resNum;
			string iCode;
			MslTools::splitIntAndString(resNumStr, resNum, iCode);
			for (vector<CharmmTopologyResidue*>::iterator k=j->begin(); k!=j->end(); k++) {
				string resName = (*k)->getName();
				for (unsigned int l=0; l<(*k)->atomSize(); l++) {
					string name;
					string type;
					double partialCharge;
					string element;
					int group;
					(*k)->getTopolAtom(l, name, type, partialCharge, element, group);
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
				};
			}
		}
	}
	// populate the system
	pSystem->addAtoms(sysAtoms);

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
				string resName = (*m)->getName();
				for (unsigned int i=0; i<(*m)->atomSize(); i++) {
					string name;
					string type;
					double partialCharge;
					string element;
					int group;
					(*m)->getTopolAtom(i, name, type, partialCharge, element, group);
					if (pSystem->atomExists(chainId, resNum, iCode, resName, name)) {
						atomMap[k-polymerDefi.begin()][l-k->begin()][m-l->begin()][name] = &(pSystem->getLastFoundAtom());
					} else {
						return false;
					}
				}
			}
		}
	}

	/************************************************************
	 *  Add the IC table
	 ************************************************************/
	vector<string> icAtoms;
	vector<double> icValues;
	bool improperFlag;
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator chItr=polymerDefi.begin(); chItr!=polymerDefi.end(); chItr++) {
		// for each chain
		for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=chItr->begin(); posItr!=chItr->end(); posItr++) {
			// for each postion
			for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin(); idItr!=posItr->end(); idItr++) {
				// for each identity
				for (unsigned int i=0; i<(*idItr)->icSize(); i++) {
					(*idItr)->getIcLine(i, icAtoms, icValues, improperFlag);
					vector<vector<Atom*> > icAtomPointers(4, vector<Atom*>());
					for (vector<string>::iterator n=icAtoms.begin(); n!=icAtoms.end(); n++) {
						if (n->substr(0,1) == "-") {
							if (posItr > chItr->begin()) {
								// not the first residue of the chain
								vector<vector<CharmmTopologyResidue*> >::iterator ll=posItr-1;
								for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
									// for each identity get all the -X atoms

									map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
									if (found == atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].end()) {
										icAtomPointers[n-icAtoms.begin()].push_back(NULL);
									} else {
										icAtomPointers[n-icAtoms.begin()].push_back(found->second);
									}

								}
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
							}
						} else if (n->substr(0,1) == "+") {
							if (posItr<chItr->end()-1) {
								// not the last residue of the chain
								vector<vector<CharmmTopologyResidue*> >::iterator ll=posItr+1;
								for (vector<CharmmTopologyResidue*>::iterator o=ll->begin(); o!=ll->end(); o++) {
									// for each identity get all the -X atoms

									map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].find(n->substr(1, n->size()-1));
									if (found == atomMap[chItr-polymerDefi.begin()][ll-chItr->begin()][o-ll->begin()].end()) {
										icAtomPointers[n-icAtoms.begin()].push_back(NULL);
									} else {
										icAtomPointers[n-icAtoms.begin()].push_back(found->second);
									}

								}
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
							}
						} else {
							map<string, Atom*>::iterator found = atomMap[chItr-polymerDefi.begin()][posItr-chItr->begin()][idItr-posItr->begin()].find(*n);
							if (found == atomMap[chItr-polymerDefi.begin()][posItr-chItr->begin()][idItr-posItr->begin()].end()) {
								icAtomPointers[n-icAtoms.begin()].push_back(NULL);
							} else {
								icAtomPointers[n-icAtoms.begin()].push_back(found->second);
							}
						}
					}
					// add the IC terms to the system
					//vector<vector<Atom*> > icAtomPointers(4, vector<Atom*>());
					for (vector<Atom*>::iterator atmK=icAtomPointers[0].begin(); atmK!=icAtomPointers[0].end(); atmK++) {
						for (vector<Atom*>::iterator atmL=icAtomPointers[1].begin(); atmL!=icAtomPointers[1].end(); atmL++) {
							for (vector<Atom*>::iterator atmM=icAtomPointers[2].begin(); atmM!=icAtomPointers[2].end(); atmM++) {
								for (vector<Atom*>::iterator atmN=icAtomPointers[3].begin(); atmN!=icAtomPointers[3].end(); atmN++) {
									if ((*atmK != NULL || *atmN != NULL) && *atmL != NULL && *atmM != NULL) {
										
										// If icValues are zero (usually zero for the PatchResidues' atoms, get the minimum values from the parameter file
										if (icValues[0] == 0.0) {
											if(improperFlag) {
												if(*atmK) {
													vector<double> param;
													if (pParReader->bondParam(param, (*atmK)->getType(),(*atmM)->getType())) {
														icValues[0] = param[1];
													}
												}
											} else {
												if(*atmK) {
													vector<double> param;
													if (pParReader->bondParam(param, (*atmK)->getType(),(*atmL)->getType())) {
														icValues[0] = param[1];
													}
												}
											}
										}	

										if (icValues[1] == 0.0) {
											if(improperFlag) {
												if(*atmK) {
													vector<double> param;
													if (pParReader->angleParam(param, (*atmK)->getType(),(*atmM)->getType(),(*atmL)->getType())) {
														icValues[1] = param[1];
													}
												}
											} else {
												if(*atmK) {
													vector<double> param;
													if (pParReader->angleParam(param, (*atmK)->getType(),(*atmL)->getType(),(*atmM)->getType())) {
														icValues[1] = param[1];
													}
												}
											}
										}
										
										if (icValues[3] == 0.0) {
											if(*atmN) {
												vector<double> param;
												if (pParReader->angleParam(param, (*atmL)->getType(),(*atmM)->getType(),(*atmN)->getType())) {
													icValues[3] = param[1];
												}
											}
										}
										
										if (icValues[4] == 0.0) {
											if(*atmN) {
												vector<double> param;
												if (pParReader->bondParam(param, (*atmM)->getType(),(*atmN)->getType())) {
													icValues[4] = param[1];
												}
											}
										}

										pSystem->addIcEntry(*atmK, *atmL, *atmM, *atmN, icValues[0], icValues[1], icValues[2], icValues[3], icValues[4], improperFlag);
									}
								}
							}
						}
					}


				}
			}
		}
	}
	

	/*********************************************************************************
	 *
	 *  ADD THE INTERACTION TERMS:
	 *
	 *  For each pair of atoms, add the appropriate interactions:
	 *   - bond if 1-2
	 *   - angle if 1-3 and autogenerate is on, or as declared
	 *   - dihedral if 1-4 and autogenerate is on, or as declared
	 *   - improper dihedrals, as declared
	 **********************************************************************************/
	EnergySet* ESet = pSystem->getEnergySet();
        
	
	// Loop over each chain.
	for (vector<vector<vector<CharmmTopologyResidue*> > >::iterator chItr=polymerDefi.begin(); chItr!=polymerDefi.end(); chItr++) {
		
		// Loop over each position.
		for (vector<vector<CharmmTopologyResidue*> >::iterator posItr=chItr->begin(); posItr!=chItr->end(); posItr++) {
		
			// Loop over each identity.
			for (vector<CharmmTopologyResidue*>::iterator idItr=posItr->begin(); idItr!=posItr->end(); idItr++) {

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
							vector<double> params;
							if (pParReader->bondParam(params, type1, type2)) {
								if (termsToBuild["CHARMM_BOND"]) {
									CharmmBondInteraction *pCBI = new CharmmBondInteraction(*(*a1),*(*a2),params[0],params[1]);
									ESet->addInteraction(pCBI);
								}
							} else {
								cerr << "WARNING: CHARMM bond parameters not found for atom types " << type1 << ", " << type2 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
							}
						}
					}
						
				}
				/************************* DONE SETTING UP THE BONDS **********************************/

				// add the angles if they are NOT autogenerated
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
										if (termsToBuild["CHARMM_U-BR"]) {
											CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*a1),*(*a3),params[0],params[1]);
											ESet->addInteraction(pCUI);
										}
									}
									if (pParReader->angleParam(params, type1, type2, type3)) {
										if (termsToBuild["CHARMM_ANGL"]) {
											CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*a1),*(*a2),*(*a3),params[0],params[1] * M_PI / 180.0); 
											ESet->addInteraction(pCAI);
										}
									} else {
										cerr << "WARNING: CHARMM angle parameters not found for atom types " << type1 << ", " << type2 << ", " << type3 << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
									}
								}
							}
						}
							
					}
				}
				/************************* DONE SETTING UP THE ANGLES **********************************/

				// add the dihedrals if they are NOT autogenerated
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
											if (termsToBuild["CHARMM_DIHE"]) {
												CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*a1),*(*a2),*(*a3),*(*a4),dihedralEntries);
												ESet->addInteraction(pCDI);
											}
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
										if (termsToBuild["CHARMM_IMPR"]) {
											CharmmImproperInteraction *pCII = new CharmmImproperInteraction(*(*a1),*(*a2),*(*a3),*(*a4),improperParams[0],improperParams[1]*M_PI/180.0);
											ESet->addInteraction(pCII);
										}
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

	// add autogenerated angles and dihedral if nececssary
	if (pTopReader->getAutoGenerateAngles()	|| pTopReader->getAutoGenerateDihedrals()) {
		AtomPointerVector atoms = pSystem->getAllAtomPointers();

		for(AtomPointerVector::iterator atomI = atoms.begin(); atomI < atoms.end(); atomI++) {
			string atomItype = (*atomI)->getType();

			for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
				if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
					continue;
				}
				string atomJtype = (*atomJ)->getType();
				if (pTopReader->getAutoGenerateAngles() && (*atomI)->isOneThree(*atomJ)) {
					// autogenerate the angle interactions
					vector<Atom*> middle = (*atomI)->getOneThreeMiddleAtoms(*atomJ);
						
					for (vector<Atom*>::iterator k=middle.begin(); k!=middle.end(); k++) {
						string middleType = (*k)->getType();
						vector<double> params;
						if (pParReader->ureyBradleyParam(params, atomItype, middleType, atomJtype)) {
							if (termsToBuild["CHARMM_U-BR"]) {
								CharmmUreyBradleyInteraction *pCUI = new CharmmUreyBradleyInteraction(*(*atomI),*(*atomJ),params[0],params[1]);
								ESet->addInteraction(pCUI);
							}
						}
						if (pParReader->angleParam(params, atomItype, middleType, atomJtype)) {
							if (termsToBuild["CHARMM_ANGL"]) {
								CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(*(*atomI),*(*k),*(*atomJ),params[0],params[1] * M_PI / 180.0); 
								ESet->addInteraction(pCAI);
							}
						} else {
							cerr << "WARNING: CHARMM angle parameters not found for atom types " << atomItype << ", " << middleType << ", " << atomJtype << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
						}
					}
				}
				if (pTopReader->getAutoGenerateDihedrals() && (*atomI)->isOneFour(*atomJ)) {
					// autogenerate the dihedral interactions
					vector<vector<Atom*> > middle = (*atomI)->getOneFourMiddleAtoms(*atomJ);
					for (vector<vector<Atom*> >::iterator k=middle.begin(); k!=middle.end(); k++) {
						vector<Atom*> middleAtoms;
						vector<string> middleTypes;
						for (vector<Atom*>::iterator l=k->begin(); l!=k->end(); l++) {
							middleAtoms.push_back(*l);
							middleTypes.push_back((*l)->getType());
						}

						vector<vector <double> > dihedralEntries;
						if (pParReader->dihedralParam(dihedralEntries, atomItype, middleTypes[0], middleTypes[1], atomJtype)) {
							// there could be multiple entries for a single dihedral
							for(int m = 0; m < dihedralEntries.size() ; m++) {
								// the delta should be expressed in radians
								dihedralEntries[m][2] = dihedralEntries[m][2]*M_PI/180.0;

							}
						
							if (termsToBuild["CHARMM_DIHE"]) {
								CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(*(*atomI),*(middleAtoms[0]),*(middleAtoms[1]),*(*atomJ),dihedralEntries);
								ESet->addInteraction(pCDI);
							}
						} else {
							cerr << "WARNING: CHARMM dihedral parameters not found for atom types " << atomItype << ", " << middleTypes[0] << ", " << middleTypes[1] << ", " << atomJtype << " in bool CharmmSystemBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)" << endl;
						}
					}
				}
			}
		}
	}

	/*********************************************************************************
	 *  ADD THE NON-BONDED INTERACTION TERMS:
	 **********************************************************************************/
	if (buildNonBondedInteractions) {
		updateNonBonded(0.0, 0.0, 0.0);
	}

	return true;
}

bool CharmmSystemBuilder::buildSystemFromPDB(string _fileName) {

	/************************************************
	 *  Read the sequence from a PDB, build the System
	 *  from topology, then assign coordinate to the
	 *  system
	 *
	 *  NOTE:  the PDB needs to be in CHARMM name
	 *         format and all atoms need to be present
	 *         or no coordinates will be assigned
	 *************************************************/

	PDBReader PR;
	if(!PR.open(_fileName)|| !PR.read()) {
		cerr << "WARNING 48293: unable to open pdb file " << _fileName << " in bool CharmmSystemBuilder::buildSystemFromPDB(string _fileName)" << endl;
		return false;
	}
	return buildSystemFromPDB(PR.getAtomPointers());
}

bool CharmmSystemBuilder::buildSystemFromPDB(const AtomPointerVector & _atoms) {

	if (pSystem == NULL) {
		cerr << "WARNING 48288: uninnitialized System inbool CharmmSystemBuilder::buildSystemFromPDB(string _fileName bool CharmmSystemBuilder::updateNonBonded(double _ctonnb, double _ctofnb, double _cutnb))" << endl;
		return false;
	}

//	PDBReader PR;
//	if(!PR.open(_fileName)|| !PR.read()) {
//		cerr << "WARNING 48293: unable to open pdb file " << _fileName << " in bool CharmmSystemBuilder::buildSystemFromPDB(string _fileName)" << endl;
//		return false;
//	}

	//PolymerSequence seq (PR.getAtomPointers());
	PolymerSequence seq (_atoms);
	buildSystem(seq);

	//pSystem->assignCoordinates(PR.getAtomPointers());
	pSystem->assignCoordinates(_atoms);
	pSystem->fillIcFromCoor();

	return true;

}

bool CharmmSystemBuilder::updateNonBonded(double _ctonnb, double _ctofnb, double _cutnb, bool _ignoreNonVariable) {
	if (pSystem == NULL) {
		cerr << "WARNING: uninnitialized System in bool CharmmSystemBuilder::updateNonBonded(double _ctonnb, double _ctofnb, double _cutnb)" << endl;
		return false;
	}
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
	EnergySet* ESet = pSystem->getEnergySet();
	ESet->eraseTerm("CHARMM_VDW");
	ESet->eraseTerm("CHARMM_ELEC");
	ESet->eraseTerm("CHARMM_EEF1");
	ESet->eraseTerm("CHARMM_EEF1REF");
	ESet->eraseTerm("CHARMM_IMM1");
	ESet->eraseTerm("CHARMM_IMM1REF");
	AtomPointerVector atoms = pSystem->getAllAtomPointers();

	bool useSolvation_local = useSolvation;
	if (!pEEF1ParReader->solventExists(solvent)) {
		useSolvation_local = false;
	}
	bool  useIMM1 = false;
	if(solvent == "MEMBRANE") {
		useIMM1 = true;
		useSolvation_local = useSolvation;
	}

	// Find which atoms are fixed (i.e. belong to non-variable positions). Interactions
	// between fixed atoms can often be safely ignored in design.
	vector<bool> fixed;
	if (_ignoreNonVariable) {
		fixed.resize(atoms.size(), true);
		for (int i = 0; i < atoms.size(); i++) {
			if (atoms[i]->getParentPosition()->getTotalNumberOfRotamers() > 1) fixed[i] = false;
		}
	}

	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 **********************************************************************/
	//unsigned int stamp = MslTools::getRandomInt(1000000);
	RandomNumberGenerator rng;
	unsigned int stamp = rng.getRandomInt(1000000);

	// To implement "smart" cutoffs, first construct a minimally bounding box around each atom
	// that encompases all of its alternative conformations
	vector<double> xmin(atoms.size()), xmax(atoms.size()), ymin(atoms.size()), ymax(atoms.size()), zmin(atoms.size()), zmax(atoms.size());
	if (_cutnb > 0) {
		for (int i = 0; i < atoms.size(); i++) {
			int ac = atoms[i]->getActiveConformation();
			atoms[i]->setActiveConformation(0);
			CartesianPoint& a = (getUseGroupCutoffs() ? atoms[i]->getGroupGeometricCenter(stamp) : atoms[i]->getCoor());
			xmin[i] = xmax[i] = a.getX();
			ymin[i] = ymax[i] = a.getY();
			zmin[i] = zmax[i] = a.getZ();
			for (int ci = 1; ci < atoms[i]->getNumberOfAltConformations(); ci++) {
				atoms[i]->setActiveConformation(ci);
				CartesianPoint& a = (getUseGroupCutoffs() ? atoms[i]->getGroupGeometricCenter(stamp) : atoms[i]->getCoor());
				if (xmin[i] > a.getX()) xmin[i] = a.getX();
				if (xmax[i] < a.getX()) xmax[i] = a.getX();
				if (ymin[i] > a.getY()) ymin[i] = a.getY();
				if (ymax[i] < a.getY()) ymax[i] = a.getY();
				if (zmin[i] > a.getZ()) zmin[i] = a.getZ();
				if (zmax[i] < a.getZ()) zmax[i] = a.getZ();
			}
			// pad the bounding box by half cutnb on all sides - this way, two boxes overlapping can be used
			// as a nececssary condition for the two corresponding atoms (sometimes) being in interaction range
			xmin[i] -= _cutnb/2; ymin[i] -= _cutnb/2; zmin[i] -= _cutnb/2;
			xmax[i] += _cutnb/2; ymax[i] += _cutnb/2; zmax[i] += _cutnb/2;
			atoms[i]->setActiveConformation(ac);
		}
	}


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
		int ai = atomI - atoms.begin();
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
		vector<double> EEF1ParamsII;  // basically the params for the other solvent
		if (useSolvation_local && !useIMM1) {
			if (!pEEF1ParReader->EEF1Param(EEF1ParamsI, atomItype, solvent)) {
				//cerr << "WARNING 49387: EEF1 parameters not found for type " << atomItype << ", solvent " << solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
				foundEEF1 = false;
			} else {
				// add the single body term
				if (termsToBuild["CHARMM_EEF1REF"]) {
					CharmmEEF1RefInteraction *pCERI = new CharmmEEF1RefInteraction(*(*atomI),EEF1ParamsI[1]);
					ESet->addInteraction(pCERI);
				}
			}
		} else if (useSolvation_local && useIMM1 && termsToBuild["CHARMM_IMM1REF"]) {
			if (pEEF1ParReader->EEF1Param(EEF1ParamsI, atomItype, "WATER") &&  pEEF1ParReader->EEF1Param(EEF1ParamsII, atomItype, "CHEX")) {
				CharmmIMM1RefInteraction *pIMM1R = new CharmmIMM1RefInteraction(*(*atomI),EEF1ParamsI[1],EEF1ParamsII[1],halfThickness,exponent);
				ESet->addInteraction(pIMM1R);
			} else {
				foundEEF1 = false;
			}
		}
		
		for(AtomPointerVector::iterator atomJ = atomI+1; atomJ < atoms.end() ; atomJ++) {
			int aj = atomJ - atoms.begin();
			if ((*atomI)->isInAlternativeIdentity(*atomJ)) {
				continue;
			}
			if (_ignoreNonVariable && fixed[ai] && fixed[aj]) continue;
			if (_cutnb > 0.0 && (!(*atomJ)->hasCoor() ||
				// sufficient condition for the two atoms never being in non-bond cutoff range range
				((xmax[ai] < xmin[aj]) || (xmax[aj] < xmin[ai]) || (ymax[ai] < ymin[aj]) || (ymax[aj] < ymin[ai]) || (zmax[ai] < zmin[aj]) || (zmax[aj] < zmin[ai])) )) {
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
			if (useSolvation_local && foundEEF1 && !special) {
				vector<double> EEF1ParamsJ;
				vector<double> EEF1ParamsJJ;
				if(!useIMM1) {
					if (pEEF1ParReader->EEF1Param(EEF1ParamsJ, atomJtype, solvent)) {
						if (termsToBuild["CHARMM_EEF1"]) {
							CharmmEEF1Interaction *pCSI = new CharmmEEF1Interaction(*(*atomI),*(*atomJ), EEF1ParamsI[0], EEF1ParamsI[2], EEF1ParamsI[5], vdwParamsI[1], EEF1ParamsJ[0], EEF1ParamsJ[2], EEF1ParamsJ[5], vdwParamsJ[1]);
							if (_cutnb > 0.0) {
								// if we are using a cutoff, set the Charmm VDW interaction with
								// the cutoffs for the switching function
								pCSI->setUseNonBondCutoffs(true, _ctonnb, _ctofnb);
							} else {
								pCSI->setUseNonBondCutoffs(false, 0.0, 0.0);
							}
							ESet->addInteraction(pCSI);
						}
					} else {
						cerr << "WARNING 49387: EEF1 parameters not found for type " << atomJtype << ", solvent " << solvent << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
					}
				} else {
					// two double body terms
					if(termsToBuild["CHARMM_IMM1"]) {
						if (pEEF1ParReader->EEF1Param(EEF1ParamsJ, atomItype, "WATER") && pEEF1ParReader->EEF1Param(EEF1ParamsJJ, atomItype, "CHEX")) {
							vector<double> IMM1ParamsW(8,0.0);
							IMM1ParamsW[0] = EEF1ParamsI[0];
							IMM1ParamsW[1] = EEF1ParamsI[2];
							IMM1ParamsW[2] = EEF1ParamsI[5];
							IMM1ParamsW[3] = vdwParamsI[1];
							IMM1ParamsW[4] = EEF1ParamsJ[0];
							IMM1ParamsW[5] = EEF1ParamsJ[2];
							IMM1ParamsW[6] = EEF1ParamsJ[5];
							IMM1ParamsW[7] = vdwParamsJ[1];

							vector<double> IMM1ParamsC(8,0.0);
							IMM1ParamsW[0] = EEF1ParamsII[0];
							IMM1ParamsW[1] = EEF1ParamsII[2];
							IMM1ParamsW[2] = EEF1ParamsII[5];
							IMM1ParamsW[3] = vdwParamsI[1];
							IMM1ParamsW[4] = EEF1ParamsJJ[0];
							IMM1ParamsW[5] = EEF1ParamsJJ[2];
							IMM1ParamsW[6] = EEF1ParamsJJ[5];
							IMM1ParamsW[7] = vdwParamsJ[1];

							CharmmIMM1Interaction *pCIMM1 = new CharmmIMM1Interaction(*(*atomI),*(*atomJ),IMM1ParamsW,IMM1ParamsC,halfThickness,exponent);
							if(_cutnb > 0.0) {
								pCIMM1->setUseNonBondCutoffs(true,_ctonnb, _ctofnb);
							} else {
								pCIMM1->setUseNonBondCutoffs(false, 0.0, 0.0);
							}
							ESet->addInteraction(pCIMM1);
						}
					}
				}
			}
			if ((*atomI)->isOneFour(*atomJ)) {
				if (!special) {
					// if it is also 1-3 or 1-2 do not add the vdw and elec term
					//vector<double> vdwParamsJ = pParReader->vdwParam(atomJtype);
					if (termsToBuild["CHARMM_ELEC"]) {
						CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant,elec14factor, useRdielectric);

						if (_cutnb > 0.0) {
							// if we are using a cutoff, set the Charmm VDW interaction with
							// the cutoffs for the switching function
							pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctofnb);
						} else {
							pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
						}
						ESet->addInteraction(pCEI);
					}
					if (foundVdw && foundVdw2) {
						if (termsToBuild["CHARMM_VDW"]) {
							CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[3]+vdwParamsJ[3]) * vdwRescalingFactor, sqrt(vdwParamsI[2] * vdwParamsJ[2]) );
							if (_cutnb > 0.0) {
								// if we are using a cutoff, set the Charmm VDW interaction with
								// the cutoffs for the switching function
								pCVI->setUseNonBondCutoffs(true, _ctonnb, _ctofnb);
							} else {
								pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
							}
							ESet->addInteraction(pCVI);
						}
						

					}
				}
				special = true;
			}
			if (!special) {
				if (termsToBuild["CHARMM_ELEC"]) {
					CharmmElectrostaticInteraction *pCEI = new CharmmElectrostaticInteraction(*(*atomI),*(*atomJ),dielectricConstant, 1.0, useRdielectric);
					if (_cutnb > 0.0) {
						// if we are using a cutoff, set the Charmm VDW interaction with
						// the cutoffs for the switching function
						pCEI->setUseNonBondCutoffs(true, _ctonnb, _ctofnb);
					} else {
						pCEI->setUseNonBondCutoffs(false, 0.0, 0.0);
					}
					ESet->addInteraction(pCEI);
				}
				if (foundVdw && foundVdw2) {
					if (termsToBuild["CHARMM_VDW"]) {
						CharmmVdwInteraction *pCVI = new CharmmVdwInteraction(*(*atomI),*(*atomJ), (vdwParamsI[1]+vdwParamsJ[1]) * vdwRescalingFactor, sqrt(vdwParamsI[0] * vdwParamsJ[0]) );
						if (_cutnb > 0.0) {
							// if we are using a cutoff, set the Charmm VDW interaction with
							// the cutoffs for the switching function
							pCVI->setUseNonBondCutoffs(true, _ctonnb, _ctofnb);
						} else {
							pCVI->setUseNonBondCutoffs(false, 0.0, 0.0);
						}
						ESet->addInteraction(pCVI);
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
	ESet->eraseTerm("CHARMM_EEF1");
	ESet->eraseTerm("CHARMM_EEF1REF");

	if (!pEEF1ParReader->solventExists(_solvent)) {
		return false;
	}

	AtomPointerVector atoms = _system.getAllAtomPointers();

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
					cerr << "WARNING 49319: VDW parameters not found for type " << atomJtype << " in bool CharmmSystemBuilder::updateSolvation(System & _system, string _solvent, double _ctonnb, double _ctofnb, double _cutnb)" << endl;
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
