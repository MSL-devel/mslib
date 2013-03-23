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

#include "PDBTopologyBuilder.h"

using namespace MSL;
using namespace std;


PDBTopologyBuilder::PDBTopologyBuilder() {
	setup();
}

PDBTopologyBuilder::PDBTopologyBuilder(System & _system, string _topologyFile) {
	setup();
	pSystem = &_system;
	bool OK = true;
	if (!readTopology(_topologyFile)) {
		OK = false;
	}
	fail_flag = !OK;
}

PDBTopologyBuilder::PDBTopologyBuilder(const PDBTopologyBuilder & _sysBuild) {
	setup();
	copy(_sysBuild);
}

PDBTopologyBuilder::~PDBTopologyBuilder() {
	delete pTopReader;
	deletePointers();
}

void PDBTopologyBuilder::operator=(const PDBTopologyBuilder & _sysBuild) {
	copy(_sysBuild);
}


void PDBTopologyBuilder::setup() {
	fail_flag = false;
	wasBuilt_flag = false;
	pSystem = NULL;
	pTopReader = new CharmmTopologyReader;
}

void PDBTopologyBuilder::copy(const PDBTopologyBuilder & _sysBuild) {
	pSystem = _sysBuild.pSystem;
	*pTopReader = *_sysBuild.pTopReader;
	fail_flag = _sysBuild.fail_flag;
	wasBuilt_flag = _sysBuild.wasBuilt_flag;
}

void PDBTopologyBuilder::deletePointers() {
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

bool PDBTopologyBuilder::addIdentity(string _positionId, string _resName, string _bbAtoms) {
	vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_positionId, vector<string>(1, _resName), bbAtoms);
}

bool PDBTopologyBuilder::addIdentity(Position & _pos, string _resName, string _bbAtoms) {
	vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_pos, vector<string>(1, _resName), bbAtoms);
}

bool PDBTopologyBuilder::addIdentity(string _positionId, const vector<string> & _resNames, string _bbAtoms) {
	vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_positionId, _resNames, bbAtoms);
}

bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, string _bbAtoms) {
	vector<string> bbAtoms = MslTools::tokenize(_bbAtoms);
	return addIdentity(_pos, _resNames, bbAtoms);;
}

bool PDBTopologyBuilder::addIdentity(string _positionId, string _resName, vector<string> _bbAtoms) {
	return addIdentity(_positionId, vector<string>(1, _resName), _bbAtoms);
}

bool PDBTopologyBuilder::addIdentity(Position & _pos, string _resName, vector<string> _bbAtoms) {
	return addIdentity(_pos, vector<string>(1, _resName), _bbAtoms);
}

bool PDBTopologyBuilder::addIdentity(string _positionId, const vector<string> & _resNames, vector<string> _bbAtoms) {
	if (pSystem->positionExists(_positionId)) {
		Position & pos = pSystem->getLastFoundPosition();
		return addIdentity(pos, _resNames, _bbAtoms);
	} else {
		cerr << "WARNING 32978: position " << _positionId << " not found in System at bool PDBTopologyBuilder::addIdentity(string _positionId, const vector<string> & _resNames, vector<string> _bbAtoms)" << endl;
		return false;
	}
}

bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms) {
	/**********************************************************
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
		cerr << "WARNING 32983: position " << _pos << " not found in topology at bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms)" << endl;
		return false;
	}
	bool chainIndexFound = false;
	System * pParentSys = _pos.getParentSystem();
	unsigned int chainIndex = 0;
	unsigned int chainSize = 0;
	for (unsigned int i=0; i<pParentSys->chainSize(); i++) {
		Chain * sysChain = &(pParentSys->getChain(i));
		if (pParentChain == sysChain) {
			chainIndex = i;
			chainSize = pParentChain->positionSize();
			chainIndexFound = true;
			break;
		}
	}
	if (!chainIndexFound) {
		cerr << "WARNING 32983: chain " << *pParentChain << " not found in topology at bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms)" << endl;
		return false;
	}

	Residue & currentRes = _pos.getCurrentIdentity();
	map<string, bool> bbAtomsMap;
	for (vector<string>::iterator k=_bbAtoms.begin(); k!=_bbAtoms.end(); k++) {
		bbAtomsMap[*k] = true;
	}

	/*********************************************************
	 *  Create an atom vector with all the atoms of the new
	 *  identities and add them to the position (the Position
	 *  object will take care of creating separate residues
	 *********************************************************/
	vector<string> processedResNames; // after removing a possible patch
	for (unsigned int i=0; i<_resNames.size(); i++) {
		// Extend the polymerDefi
		vector<string> split = MslTools::tokenize(_resNames[i], "-");
		processedResNames.push_back(split[0]);
		if (!pTopReader->residueExists(split[0])) {
			cerr << "WARNING 32988: residue " << split[0] << " not found in topology at bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms)" << endl;
			return false;
		}

		CharmmTopologyResidue * pTopRes = new CharmmTopologyResidue(pTopReader->getLastFoundResidue());
		if (wasBuilt_flag) {
			polymerDefi[chainIndex][posIndex].push_back(pTopRes);
			atomMap[chainIndex][posIndex].push_back(map<string, Atom*>());
		}

		// if the patches for the first or the last residue were not declared
		// get them from the residues
		if (split.size() == 1 && posIndex == 0) {
			// add default first residue patch to current identity
			//split.push_back(polymerDefi[chainIndex][posIndex].back()->getFirstDefaultPatch());
			split.push_back(pTopRes->getFirstDefaultPatch());
		}
		//if (split.size() == 1 && posIndex == polymerDefi[chainIndex].size()-1) {
		if (split.size() == 1 && posIndex == chainSize-1) {
			// add default last residue patch to current identity
			//split.push_back(polymerDefi[chainIndex][posIndex].back()->getLastDefaultPatch());
			split.push_back(pTopRes->getLastDefaultPatch());
		}
		/************************************************************
		 *  APPLY THE PATCHES (if any)
		 ************************************************************/
		for (vector<string>::iterator n=split.begin()+1; n!=split.end(); n++) {
			// patch current identity
			if (pTopReader->residueExists(*n)) {
				//if (!polymerDefi[chainIndex][posIndex].back()->applyPatch(pTopReader->getLastFoundResidue())) {
				if (!pTopRes->applyPatch(pTopReader->getLastFoundResidue())) {
					cerr << "WARNING 19134: cannot apply patch " << *n << " to residue, in bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms)";
					return false;
				}
			}
		}

		// Add the atoms to the position
		string chainId = _pos.getChainId();
		int resNum = _pos.getResidueNumber();
		string iCode = _pos.getResidueIcode();
		AtomPointerVector posAtoms;
		//for (unsigned int l=0; l<polymerDefi[chainIndex][posIndex].back()->atomSize(); l++) {
		for (unsigned int l=0; l<pTopRes->atomSize(); l++) {
			string name;
			string type;
			double partialCharge;
			string element;
			int group;
			//polymerDefi[chainIndex][posIndex].back()->getTopolAtom(l, name, type, partialCharge, element, group);
			pTopRes->getTopolAtom(l, name, type, partialCharge, element, group);
			Atom * a = new Atom(name);
			a->setChainId(chainId);
			a->setResidueName(_resNames[i]);
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

		if (!wasBuilt_flag) {
			// add an ic table and try to build the new atoms (only if the backbone atoms were given)
			if (_bbAtoms.size() != 0) {
				Residue * pRes = &_pos.getResidue(_resNames[i]);
				AtomPointerVector atoms = pRes->getAtomPointers();
				vector<string> icAtoms;
				vector<double> icValues;
				bool improperFlag;
				for (unsigned int i=0; i<pTopRes->icSize(); i++) {
					pTopRes->getIcLine(i, icAtoms, icValues, improperFlag);
					vector<Atom*> icAtomPointers(4, (Atom*)NULL);
					for (unsigned int j=0; j<icAtoms.size(); j++) {
						if (icAtoms[j].substr(0,1) == "-" || icAtoms[j].substr(0,1) == "+") {
							icAtomPointers[j] = NULL;
						} else {
							for (unsigned int k=0; k<atoms.size(); k++) {
								if (atoms[k]->getName() == icAtoms[j]) {
									icAtomPointers[j] = atoms[k];
									break;
								}
							}
						}
					}
					// add the IC terms to the system
					if ((icAtomPointers[0] != NULL || icAtomPointers[3] != NULL) && icAtomPointers[1] != NULL && icAtomPointers[2] != NULL) {
						_pos.getParentSystem()->addIcEntry(icAtomPointers[0], icAtomPointers[1], icAtomPointers[2], icAtomPointers[3], icValues[0], icValues[1], icValues[2], icValues[3], icValues[4], improperFlag);
					}


				}
				for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
					(*k)->buildFromIc(false); // build only from active atoms = false
				}
			}
			delete pTopRes;
			pTopRes = NULL;
		}

	}

	if (wasBuilt_flag) {

		/*********************************************************
		 *  These variables will be used to cycle in the position where we
		 *  will add the identities
		 *********************************************************/
		unsigned int idIndexStart = polymerDefi[chainIndex][posIndex].size();
		unsigned int idIndexEnd = idIndexStart + _resNames.size() - 1;

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
					cerr << "WARNING 32998: atom " << processedResNames[i-idIndexStart] << "," << name << " not found in topology at bool PDBTopologyBuilder::addIdentity(Position & _pos, const vector<string> & _resNames, vector<string> _bbAtoms)" << endl;
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
		//vector<CharmmTopologyResidue*>::iterator idThis=posThis->end()-1;

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
										_pos.getParentSystem()->addIcEntry(*atmK, *atmL, *atmM, *atmN, icValues[0], icValues[1], icValues[2], icValues[3], icValues[4], improperFlag);
									}
								}
							}
						}
					}


				}
			}
		}
	}
	


	return true;
}


bool PDBTopologyBuilder::buildSystem(const PolymerSequence & _sequence) {
	//pSystem = &_system;

	if (pSystem == NULL) {
		cerr << "WARNING, undefined System in bool PDBTopologyBuilder::buildSystem(const PolymerSequence & _sequence), set System first" << endl;
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
					cerr << "WARNING 19129: residue " << split[0] << " does not exist in topology, in bool PDBTopologyBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
					return false;
				}
				// if the patches for the first or the last residue were not declared
				// get them from the residues
				if (split.size() == 1 && l==k->begin()) {
					// add default first residue patch to current identity
					split.push_back(polymerDefi.back().back().back()->getFirstDefaultPatch());
				}
				if (split.size() == 1 && l==k->end() - 1) {
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
							cerr << "WARNING 19134: cannot apply patch " << *n << " to residue, in bool PDBTopologyBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
							return false;
						}
					} else {
						cerr << "WARNING 19139: cannot apply unexistent patch " << *n << " to residue, in bool PDBTopologyBuilder::buildSystem(System & _system, const PolymerSequence & _sequence)";
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
										
										/*
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
										*/

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
	 *  ADD THE BONDS:
	 *
	 **********************************************************************************/
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
						for(vector<Atom*>::iterator a2 = pAtom2.begin() ; a2 != pAtom2.end(); a2++) {
							(*a1)->setBoundTo(*a2);
						}
					}
						
				}
			}
		}
	}

	wasBuilt_flag = true;
	
	return true;
}

bool PDBTopologyBuilder::buildSystemFromPDB(string _fileName) {

	/************************************************
	 *  Read the sequence from a PDB, build the System
	 *  from topology, then assign coordinate to the
	 *  system
	 *
	 *  NOTE:  the PDB needs to be in CHARMM name
	 *         format and all atoms need to be present
	 *         or no coordinates will be assigned
	 *************************************************/

	if (pSystem == NULL) {
		cerr << "WARNING 48288: uninnitialized System inbool PDBTopologyBuilder::buildSystemFromPDB(string _fileName bool PDBTopologyBuilder::updateNonBonded(double _ctonnb, double _ctofnb, double _cutnb))" << endl;
		return false;
	}

	PDBReader PR;
	if(!PR.open(_fileName)|| !PR.read()) {
		cerr << "WARNING 48293: unable to open pdb file " << _fileName << " in bool PDBTopologyBuilder::buildSystemFromPDB(string _fileName)" << endl;
		return false;
	}

	PolymerSequence seq(PR.getAtomPointers());
	if (!buildSystem(seq)) {
		return false;
	}

	pSystem->assignCoordinates(PR.getAtomPointers());
	pSystem->fillIcFromCoor();

	return true;

}

void PDBTopologyBuilder::getAtomPointersFromMulti(string _name, vector<Atom*> & _out, vector<CharmmTopologyResidue*> & _position, vector<map<string, Atom*> > & _atomMap) {
	for(vector<CharmmTopologyResidue*>::iterator idItr = _position.begin(); idItr != _position.end(); idItr++ ) {
		map <string, Atom*>::iterator found =  _atomMap[idItr-_position.begin()].find(_name);
		if (found != _atomMap[idItr-_position.begin()].end()) {	
			_out.push_back(_atomMap[idItr-_position.begin()][_name]);
		}	
	}
}

vector<Atom*> PDBTopologyBuilder::getAtomPointers(string _name, vector<vector<vector<CharmmTopologyResidue*> > >::iterator & _chItr, vector<vector<CharmmTopologyResidue*> >::iterator & _posItr, vector<CharmmTopologyResidue*>::iterator & _idItr) {
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
