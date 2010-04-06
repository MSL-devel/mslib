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


#include "System.h"
#include "PolymerSequence.h"

using namespace MSL;
using namespace std;


System::System() {
	setup();
}

System::System(const Chain & _chain) {
	setup();
	addChain(_chain);
}

System::System(const AtomPointerVector & _atoms) {
	setup();
	addAtoms(_atoms);
}

System::System(const System & _system) {
	setup();
	copy(_system);
}

System::~System() {
	deletePointers();
}

void System::operator=(const System & _system) {
	copy(_system);
}


void System::setup() {
	pdbReader = new PDBReader;
	pdbWriter = new PDBWriter;
	ESet = new EnergySet; //TODO update System Copy
	noUpdateIndex_flag = false;
	foundChain = chainMap.end();
	nameSpace = "";
}

void System::copy(const System & _system) {
	/**************************************************************
	 *  TODO:  add copy of bonds
	 **************************************************************/
	reset();
	cerr << "System::copy() doesnt handle copying the EnergySet" << endl;

	for (vector<Chain*>::const_iterator k=_system.chains.begin(); k!=_system.chains.end(); k++) {
		chains.push_back(new Chain(**k));
		chains.back()->setParentSystem(this);
		chainMap[chains.back()->getChainId()] = chains.back();
	}
	foundChain = chainMap.end();
	updateIndexing();
	updateAllAtomIndexing();
	
	/************************************************
	 *  Copy IC table 
	 ************************************************/
	for (IcTable::iterator k=icTable.begin(); k!=icTable.end(); k++) {
		delete *k;
	}
	icTable.clear();
	for (IcTable::const_iterator k=_system.icTable.begin(); k!=_system.icTable.end(); k++) {
		//char c [1000];
		Atom * pAtom1 = NULL;
		Atom * pAtom2 = NULL;
		Atom * pAtom3 = NULL;
		Atom * pAtom4 = NULL;
		Atom * sys_pAtom1 = (*k)->getAtom1();
		Atom * sys_pAtom2 = (*k)->getAtom2();
		Atom * sys_pAtom3 = (*k)->getAtom3();
		Atom * sys_pAtom4 = (*k)->getAtom4();
		if (sys_pAtom1 != NULL) {
			if (atomExists(sys_pAtom1->getChainId(), sys_pAtom1->getResidueNumber(), sys_pAtom1->getResidueIcode(), sys_pAtom1->getResidueName(), sys_pAtom1->getName())) {
				pAtom1 = &getLastFoundAtom();
			}
		}
		if (sys_pAtom2 != NULL) {
			if (atomExists(sys_pAtom2->getChainId(), sys_pAtom2->getResidueNumber(), sys_pAtom2->getResidueIcode(), sys_pAtom2->getResidueName(), sys_pAtom2->getName())) {
				pAtom2 = &getLastFoundAtom();
			}
		}
		if (sys_pAtom3 != NULL) {
			if (atomExists(sys_pAtom3->getChainId(), sys_pAtom3->getResidueNumber(), sys_pAtom3->getResidueIcode(), sys_pAtom3->getResidueName(), sys_pAtom3->getName())) {
				pAtom3 = &getLastFoundAtom();
			}
		}
		if (sys_pAtom4 != NULL) {
			if (atomExists(sys_pAtom4->getChainId(), sys_pAtom4->getResidueNumber(), sys_pAtom4->getResidueIcode(), sys_pAtom4->getResidueName(), sys_pAtom4->getName())) {
				pAtom4 = &getLastFoundAtom();
			}
		}


		//if (findIcAtoms(pAtom1, pAtom2, pAtom3, pAtom4, _1_chain, _1_resNumIcode, _1_name, _2_chain, _2_resNumIcode, _2_name, _3_chain, _3_resNumIcode, _3_name, _4_chain, _4_resNumIcode, _4_name)) {
		if (pAtom2 != NULL && pAtom3 != NULL && (pAtom1 != NULL || pAtom4 != NULL)) {
			// 2 and 3 must exists and at least one between 1 and 4
			// add the ic entry
			icTable.push_back(new IcEntry(*pAtom1, *pAtom2, *pAtom3, *pAtom4, (*k)->getDistance1(), (*k)->getAngle1(), (*k)->getDihedral(), (*k)->getAngle2(), (*k)->getDistance2(), (*k)->isImproper()));
			// add any saved internal coordinate buffers
			icTable.back()->setStoredValues((*k)->getStoredValues());
		}

	}
}

void System::reset() {
	deletePointers();
	setup();
}

void System::deletePointers() {
	activeAtoms.clear();
	activeAndInactiveAtoms.clear();
	positions.clear();
	// avoid calling updates for deleted chains
	noUpdateIndex_flag = true;
	for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {
		delete *k;
		*k = NULL;
	}
	chains.clear();
	chainMap.clear();

	for (IcTable::iterator k=icTable.begin(); k!=icTable.end(); k++) {
		delete *k;
		*k = NULL;
	}
	icTable.clear();

	delete pdbReader;
	pdbReader = NULL;
	delete pdbWriter;
	pdbWriter = NULL;
	delete ESet;
	ESet = NULL;
	// reset the flag and run the updates
	noUpdateIndex_flag = false;
	updateIndexing();
	updateAllAtomIndexing();
	//delete polSeq;
}

void System::addChain(const Chain & _chain, string _chainId) {
	// how to handle the error if the chain id is already used?
	chains.push_back(new Chain(_chain));
	if (_chainId == "") {
		string chainIdSequence ="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		_chainId = _chain.getChainId();
		for (unsigned int i=0; i<chainIdSequence.size(); i++) {
			bool found = false;
			for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {
				if ((*k)->getChainId() == chainIdSequence.substr(i,1)) {
					found = true;
					break;
				}
			}
			if (!found) {
				_chainId = chainIdSequence.substr(i,1);
				break;
			}
		}

	}
	chains.back()->setChainId(_chainId);
	chains.back()->setParentSystem(this);
	chainMap[chains.back()->getChainId()] = chains.back();
	vector<Position*> positions = chains.back()->getPositions();
	updateIndexing();
	updateAllAtomIndexing();
	/*
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		AtomPointerVector resAtoms = (*k)->getAtomPointers();
		ResidueAtoms tmp;
		tmp.start = activeAtoms.size();
		tmp.size = resAtoms.size();
		residueLookupMap[*k] = tmp;
		activeAtoms.insert(activeAtoms.end(), resAtoms.begin(), resAtoms.end());
	}
	*/

}

bool System::removeChain(string _chainId) {


	foundChain=chainMap.find(_chainId);

	if (foundChain!=chainMap.end()) {
		for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {
			if ((*foundChain).second == *k) {
				noUpdateIndex_flag = true;

				// remove from the active atoms and the residue lookup map
				/*
				vector<Position*> positions = (*k)->getPositions();
				for (vector<Position*>::iterator l=positions.begin(); l!=positions.end(); l++) {
					map<Position*, ResidueAtoms>::iterator found=residueLookupMap.find(*l);
					if (found!=residueLookupMap.end()) {
						if (found->second.size > 0) {
							activeAtoms.erase(activeAtoms.begin() + found->second.start, activeAtoms.begin() + found->second.start + found->second.size - 1);
						}
						residueLookupMap.erase(found);
					}
				}
				*/

				// deallocate from memory and remove from list
				delete *k;
				chains.erase(k);
				chainMap.erase(foundChain);
				foundChain = chainMap.end();
				return true;
				noUpdateIndex_flag = false;
			}
		}
		updateIndexing();
		updateAllAtomIndexing();
	}
	return false;
}

void System::removeAllChains() {
	deletePointers();
}

bool System::duplicateChain(string _chainId, string _newChainId) {
	//if (!exists(_chainId)) {
	if (!chainExists(_chainId)) {
		return false;
	}
	addChain(*(foundChain->second), _newChainId);
	return true;
}

bool System::duplicateChain(unsigned int _n, string _newChainId) {
	if (_n >= chains.size()) {
		return false;
	}
	addChain(*(chains[_n]), _newChainId);
	return true;
}

void System::addAtoms(const AtomPointerVector & _atoms) {
	/***********************************************
	 *  This function splits the AtomPointerVector in a number
	 *  of AtomPointerVector objects by chain ID and then calls
	 *  the addAtoms(const AtomPointerVector & _atoms) function
	 *  of the chains to take care of rhe rest
	 ***********************************************/
	
	noUpdateIndex_flag = true;

	map<string, AtomPointerVector> dividedInChains;

	// store the order of the chains so that it will be preserved
	vector<string> chainOrder;
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		if (dividedInChains.find((*k)->getChainId()) == dividedInChains.end()) {
			chainOrder.push_back((*k)->getChainId());
		}
		dividedInChains[(*k)->getChainId()].push_back(*k);
	}

	int counter = 0;
	/*
	for (map<string, AtomPointerVector>::iterator k=dividedInChains.begin(); k!=dividedInChains.end(); k++) {
		counter++;
		//cout << "UUU chain " << counter << " of " << dividedInChains.size() << endl;
		foundChain=chainMap.find(k->first);
		if (foundChain!=chainMap.end()) {
			/ ***********************************************
			 *  A chain with the chainId already EXISTS, add
			 *  the atoms to it
			 *********************************************** /
			(*foundChain).second->addAtoms(k->second);
		} else {
			/ ***********************************************
			 *  A chain with the chainId DOES NOT EXIST, 
			 *  create a new chain first and add
			 *  the atoms to it
			 *********************************************** /
			chains.push_back(new Chain(k->first));
			chains.back()->setParentSystem(this);
			chainMap[chains.back()->getChainId()] = chains.back();
			chains.back()->addAtoms(k->second);
		}
	}
	*/
	// changed code to keep the original order of the chains (the map was loosing it)
	for (unsigned int i=0; i<chainOrder.size(); i++) {
		counter++;
		//cout << "UUU chain " << counter << " of " << dividedInChains.size() << endl;
		foundChain=chainMap.find(chainOrder[i]);
		if (foundChain!=chainMap.end()) {
			/***********************************************
			 *  A chain with the chainId already EXISTS, add
			 *  the atoms to it
			 ***********************************************/
			(*foundChain).second->addAtoms(dividedInChains[chainOrder[i]]);
		} else {
			/***********************************************
			 *  A chain with the chainId DOES NOT EXIST, 
			 *  create a new chain first and add
			 *  the atoms to it
			 ***********************************************/
			chains.push_back(new Chain(chainOrder[i]));
			chains.back()->setParentSystem(this);
			chainMap[chains.back()->getChainId()] = chains.back();
			chains.back()->addAtoms(dividedInChains[chainOrder[i]]);
		}
	}


	noUpdateIndex_flag = false;
	updateIndexing();
	updateAllAtomIndexing();

}

void System::updateChainMap(Chain * _chain) {
	// if a chain changes its name it calls this function
	// to update the map
	for (map<string, Chain*>::iterator k=chainMap.begin(); k!=chainMap.end(); k++) {
		if (k->second == _chain) {
			if (k->first != _chain->getChainId()) {
				chainMap.erase(k);
				chainMap[_chain->getChainId()] = _chain;
			}
		}
	}
}

/*
void System::swapInActiveList(Position * _position, AtomPointerVector & _atoms) {
	map<Position *, ResidueAtoms>::iterator found = residueLookupMap.find(_position);
	if (found != residueLookupMap.end()) {
		if (found->second.size > 0) {
			activeAtoms.erase(activeAtoms.begin() + found->second.start, activeAtoms.begin() + found->second.start + found->second.size - 1);
		}
		if (_atoms.size() > 0) {
			activeAtoms.insert(activeAtoms.begin() + found->second.start, _atoms.begin(), _atoms.end());
		}
	}
	
}
*/

void System::updateIndexing() {
	if (noUpdateIndex_flag) {
		return;
	}
	activeAtoms.clear();
	for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {
		activeAtoms.insert(activeAtoms.end(), (*k)->getAtomPointers().begin(), (*k)->getAtomPointers().end());

	}
}


void System::updateAllAtomIndexing() {
	if (noUpdateIndex_flag) {
		return;
	}
	activeAndInactiveAtoms.clear();
	positions.clear();
	for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {
		activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), (*k)->getAllAtomPointers().begin(), (*k)->getAllAtomPointers().end());
		positions.insert(positions.end(), (*k)->getPositions().begin(), (*k)->getPositions().end());

	}
}


bool System::addIcEntry(string _atomId1, string _atomId2, string _atomId3, string _atomId4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {
	/*
	vector<string> tokens = MslTools::tokenize(_1_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain1 = tokens[0];
	string resNumIcode1 = tokens[1];
	string name1 = tokens[2];

	tokens = MslTools::tokenize(_2_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain2 = tokens[0];
	string resNumIcode2 = tokens[1];
	string name2 = tokens[2];

	tokens = MslTools::tokenize(_3_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain3 = tokens[0];
	string resNumIcode3 = tokens[1];
	string name3 = tokens[2];

	tokens = MslTools::tokenize(_4_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain4 = tokens[0];
	string resNumIcode4 = tokens[1];
	string name4 = tokens[2];

	return addIcEntry(chain1, resNumIcode1, name1, chain2, resNumIcode2, name2, chain3, resNumIcode3, name3, chain4, resNumIcode4, name4, _d1, _a1, _dihe, _a2, _d2); 
	*/

	Atom * pAtom1 = NULL;
	if (atomExists(_atomId1)) {
		pAtom1 = &getLastFoundAtom();
	}

	Atom * pAtom2 = NULL;
	if (atomExists(_atomId2)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3 = NULL;
	if (atomExists(_atomId3)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom4 = NULL;
	if (atomExists(_atomId4)) {
		pAtom4 = &getLastFoundAtom();
	}
	if (pAtom1 == NULL && pAtom4 == NULL) {
		// at least one must exist
		return false;
	}
	return addIcEntry(pAtom1, pAtom2, pAtom3, pAtom4, _d1, _a1, _dihe, _a2, _d2); 

}

bool System::addIcEntry(string _chain1, unsigned int _resnum1, string _icode1, string _atomName1, string _chain2, unsigned int _resnum2, string _icode2, string _atomName2, string _chain3, unsigned int _resnum3, string _icode3, string _atomName3, string _chain4, unsigned int _resnum4, string _icode4, string _atomName4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {
	Atom * pAtom1 = NULL;
	if (atomExists(_chain1, _resnum1, _icode1, _atomName1)) {
		pAtom1 = &getLastFoundAtom();
	}

	Atom * pAtom2 = NULL;
	if (atomExists(_chain2, _resnum2, _icode2, _atomName2)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3 = NULL;
	if (atomExists(_chain3, _resnum3, _icode3, _atomName3)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom4 = NULL;
	if (atomExists(_chain4, _resnum4, _icode4, _atomName4)) {
		pAtom4 = &getLastFoundAtom();
	}
	if (pAtom1 == NULL && pAtom4 == NULL) {
		// at least one must exist
		return false;
	}
	return addIcEntry(pAtom1, pAtom2, pAtom3, pAtom4, _d1, _a1, _dihe, _a2, _d2); 
}

bool System::addIcEntry(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {

	cerr << "DEPRECATED FUNCTION bool System::addIcEntry(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag)" << endl;
	Atom * pAtom1 = NULL;
	Atom * pAtom2 = NULL;
	Atom * pAtom3 = NULL;
	Atom * pAtom4 = NULL;

	if (findIcAtoms(pAtom1, pAtom2, pAtom3, pAtom4, _1_chain, _1_resNumIcode, _1_name, _2_chain, _2_resNumIcode, _2_name, _3_chain, _3_resNumIcode, _3_name, _4_chain, _4_resNumIcode, _4_name)) {
		addIcEntry(pAtom1, pAtom2, pAtom3, pAtom4, _d1, _a1, _dihe, _a2, _d2, _improperFlag);
		return true;
	} else {
		return false;
	}
}

bool System::addIcEntry(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {
	icTable.push_back(new IcEntry(*_pAtom1, *_pAtom2, *_pAtom3, *_pAtom4, _d1, _a1, _dihe, _a2, _d2, _improperFlag));
	return true;
}

bool System::seed(string _atomId1, string _atomId2, string _atomId3) {
	Atom * pAtom1 = NULL;
	if (atomExists(_atomId1)) {
		pAtom1 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom2 = NULL;
	if (atomExists(_atomId2)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3 = NULL;
	if (atomExists(_atomId3)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}
	return seed(pAtom1, pAtom2, pAtom3);
	/*
	vector<string> tokens = MslTools::tokenize(_1_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain1 = tokens[0];
	string resNumIcode1 = tokens[1];
	string name1 = tokens[2];

	tokens = MslTools::tokenize(_2_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain2 = tokens[0];
	string resNumIcode2 = tokens[1];
	string name2 = tokens[2];

	tokens = MslTools::tokenize(_3_chain_resNumIcode_name);
	if (tokens.size() != 3) {
		return false;
	}
	string chain3 = tokens[0];
	string resNumIcode3 = tokens[1];
	string name3 = tokens[2];

	return seed(chain1, resNumIcode1, name1, chain2, resNumIcode2, name2, chain3, resNumIcode3, name3); 
	*/
}

bool System::seed(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name) {

	cerr << "DEPRECATED FUNCTION bool System::seed(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name)" << endl;
	Atom * pAtom1 = NULL;
	if (exists(_1_chain, _1_resNumIcode, _1_name)) {
		pAtom1 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom2 = NULL;
	if (exists(_2_chain, _2_resNumIcode, _2_name)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3 = NULL;
	if (exists(_3_chain, _3_resNumIcode, _3_name)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}
	return seed(pAtom1, pAtom2, pAtom3);
}

bool System::seed(string _chain1, unsigned int _resnum1, string _icode1, string _atomName1, string _chain2, unsigned int _resnum2, string _icode2, string _atomName2, string _chain3, unsigned int _resnum3, string _icode3, string _atomName3) {
	Atom * pAtom1 = NULL;
	if (atomExists(_chain1, _resnum1, _icode1, _atomName1)) {
		pAtom1 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom2 = NULL;
	if (atomExists(_chain2, _resnum2, _icode2, _atomName2)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3 = NULL;
	if (atomExists(_chain3, _resnum3, _icode3, _atomName3)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}
	return seed(pAtom1, pAtom2, pAtom3);
}

bool System::seed(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) {
	// removes all coordinates and finds 3 atoms to seed in cartesian space

	/*
	Atom * pAtom1;
	if (exists(_1_chain, _1_resNumIcode, _1_name)) {
		pAtom1 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom2;
	if (exists(_2_chain, _2_resNumIcode, _2_name)) {
		pAtom2 = &getLastFoundAtom();
	} else {
		return false;
	}

	Atom * pAtom3;
	if (exists(_3_chain, _3_resNumIcode, _3_name)) {
		pAtom3 = &getLastFoundAtom();
	} else {
		return false;
	}
	*/

	if (_pAtom1 == NULL || _pAtom2 == NULL || _pAtom3 == NULL) {
		return false;
	}

	double d1 = 0.0;
	double d2 = 0.0;
	double a = 0.0;
	bool foundD1 = false;
	bool foundD2 = false;
	bool foundA = false;
	for (IcTable::iterator k=icTable.begin(); k!=icTable.end(); k++) {
		Atom * IC1 = (*k)->getAtom1();
		Atom * IC2 = (*k)->getAtom2();
		Atom * IC3 = (*k)->getAtom3();
		Atom * IC4 = (*k)->getAtom4();
		if (IC4 == _pAtom1 && IC3 == _pAtom2) {
			d1 = (*k)->getDistance2();
			foundD1 = true;
			if (IC2 == _pAtom3) {
				a  = (*k)->getAngle2Radians();
				foundA = true;
			}
	//		cout << "UUU F1 as rev-dir" << endl;
			if (foundD2) {
				break;
			}
		} else if (IC4 == _pAtom3 && IC3 == _pAtom2) {
			d2 = (*k)->getDistance2();
			foundD2 = true;
			if (IC2 == _pAtom1) {
				a  = (*k)->getAngle2Radians();
				foundA = true;
			}
	//		cout << "UUU F2 as rev-inv" << endl;
			if (foundD1) {
				break;
			}
		} else if (!(*k)->isImproper() && IC1 == _pAtom1 && IC2 == _pAtom2) {
			d1 = (*k)->getDistance1();
			foundD1 = true;
			if (IC3 == _pAtom3) {
				a  = (*k)->getAngle1Radians();
				foundA = true;
			}
	//		cout << "UUU F1 as regfor-dir" << endl;
			if (foundD2) {
				break;
			}
		} else if (!(*k)->isImproper() && IC1 == _pAtom3 && IC2 == _pAtom2) {
			d2 = (*k)->getDistance1();
			foundD2 = true;
			if (IC3 == _pAtom1) {
				a  = (*k)->getAngle1Radians();
				foundA = true;
			}
	//		cout << "UUU F2 as regfor-inv" << endl;
			if (foundD1) {
				break;
			}
		} else if ((*k)->isImproper() && IC1 == _pAtom1 && IC3 == _pAtom2) {
			d1 = (*k)->getDistance1();
			foundD1 = true;
			if (IC2 == _pAtom3) {
				a  = (*k)->getAngle1Radians();
				foundA = true;
			}
	//		cout << "UUU F1 as impfor-dir" << endl;
			if (foundD2) {
				break;
			}
		} else if ((*k)->isImproper() && IC1 == _pAtom3 && IC3 == _pAtom2) {
			d2 = (*k)->getDistance1();
			foundD2 = true;
			if (IC2 == _pAtom1) {
				a  = (*k)->getAngle1Radians();
				foundA = true;
			}
	//		cout << "UUU F2 as impfor-inv" << endl;
			if (foundD1) {
				break;
			}
		} else if (IC4 == _pAtom2 && IC3 == _pAtom3) {
			d2 = (*k)->getDistance2();
			foundD2 = true;
	//		cout << "UUU F2 as rev-dir" << endl;
			if (foundD1 && foundA) {
				break;
			}
		} else if (IC4 == _pAtom2 && IC3 == _pAtom1) {
			d1 = (*k)->getDistance2();
			foundD1 = true;
	//		cout << "UUU F1 as rev-inv" << endl;
			if (foundD2 && foundA) {
				break;
			}
		} else if (!(*k)->isImproper() && IC1 == _pAtom2 && IC2 == _pAtom3) {
			d2 = (*k)->getDistance1();
			foundD2 = true;
	//		cout << "UUU F2 as regfor-dir" << endl;
			if (foundD1 && foundA) {
				break;
			}
		} else if (!(*k)->isImproper() && IC1 == _pAtom2 && IC2 == _pAtom1) {
			d1 = (*k)->getDistance1();
			foundD1 = true;
	//		cout << "UUU F1 as regfor-inv" << endl;
			if (foundD2 && foundA) {
				break;
			}
		} else if ((*k)->isImproper() && IC1 == _pAtom2 && IC3 == _pAtom3) {
			d2 = (*k)->getDistance1();
			foundD2 = true;
	//		cout << "UUU F2 as impfor-dir" << endl;
			if (foundD1 && foundA) {
				break;
			}
		} else if ((*k)->isImproper() && IC1 == _pAtom2 && IC3 == _pAtom1) {
			d1 = (*k)->getDistance1();
			foundD1 = true;
	//		cout << "UUU F1 as impfor-inv" << endl;
			if (foundD2 && foundA) {
				break;
			}
		}
	}
	if (foundD1 && foundD2 && foundA) {
		CartesianGeometry::instance()->seedRadians(_pAtom1->getCoor(), _pAtom2->getCoor(), _pAtom3->getCoor(), d1, d2, a);
		_pAtom1->setHasCoordinates();
		_pAtom2->setHasCoordinates();
		_pAtom3->setHasCoordinates();
	//	cout << "UUU seeded " << _pAtom1->getCoor() << _pAtom2->getCoor() << _pAtom3->getCoor() << " " << d1 << " " << d2 << " " << a << endl;
		return true;
	} else {
		return false;
	}
}

bool System::findIcAtoms(Atom *& _pAtom1, Atom *& _pAtom2, Atom *& _pAtom3, Atom *& _pAtom4, string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name) {

	cerr << "DEPRECATED FUNCTION bool System::findIcAtoms(Atom *& _pAtom1, Atom *& _pAtom2, Atom *& _pAtom3, Atom *& _pAtom4, string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name)" << endl;

	/*****************************************************************
	 *  Given chain, residue and atom info, find the pointer to the
	 *  four atoms in an IC table
	 *
	 *  Used to add IC entries
	 *
	 *  This function receives Atom pointer variables by reference (*&),
	 *  or else the address would be reassigned to the found atom but the
	 *  container (the Atom pointer) would be a local copy and would not
	 *  return the required Atom address
	 *****************************************************************/

	bool A1_exists = false;
	bool A4_exists = false;

	if (exists(_1_chain, _1_resNumIcode, _1_name)) {
		_pAtom1 = &getLastFoundAtom();
		A1_exists = true;
	} else {
		_pAtom1 = NULL;
	}

	if (exists(_4_chain, _4_resNumIcode, _4_name)) {
		_pAtom4 = &getLastFoundAtom();
		A4_exists = true;
	} else {
		_pAtom4 = NULL;
	}

	if (!A1_exists && !A4_exists) {
		// at least 1 or 4 must exist
		return false;
	}

	if (exists(_2_chain, _2_resNumIcode, _2_name)) {
		_pAtom2 = &getLastFoundAtom();
	} else {
		// 2 must exist
		_pAtom2 = NULL;
		return false;
	}

	if (exists(_3_chain, _3_resNumIcode, _3_name)) {
		_pAtom3 = &getLastFoundAtom();
	} else {
		// 3 must exist
		_pAtom3 = NULL;
		return false;
	}

	return true;
}

unsigned int System::getPositionIndex(const Position * _pPos) const {
	for (vector<Position*>::const_iterator k=positions.begin(); k!=positions.end(); k++) {
		if (_pPos == *k) {
			return k-positions.begin();
		}
	}
	cerr << "ERROR 44193: Position not found in unsigned int System::getPositionIndex(Position * _pPos)" << endl;
	exit(44193);
}

unsigned int System::assignCoordinates(const AtomPointerVector & _atoms, bool checkIdentity) {
	// only set coordinates for existing matching atoms in the system, ignore the rest of the atoms
	// returns the number of atoms assigned	
	
	unsigned int counter = 0;
	
	string chainId = "";
	//string resNumAndIcode = "";
	string name = "";
	string identity = "";
	int resnum = 0;
	string icode;
	//char resNumAndIcode [1000];
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		chainId = (*k)->getChainId();
		//sprintf(resNumAndIcode, "%d%s", (*k)->getResidueNumber(), (*k)->getResidueIcode().c_str());
		resnum = (*k)->getResidueNumber();
		icode = (*k)->getResidueIcode();
		name = (*k)->getName();
		identity = (*k)->getResidueName();

		Atom *pAtom = NULL;
		if (checkIdentity) {
			// apply coordinates only to atoms that have the same residue name
			//if (exists(chainId, resNumAndIcode, name, identity)) {
			if (atomExists(chainId, resnum, icode, identity, name)) {
				pAtom = &(getLastFoundAtom());
				pAtom->setCoor((*k)->getCoor());
				counter++;
			} 
		} else {
			// Apply to all idenitites that have an atom with this name
			//if (exists(chainId, resNumAndIcode, name)) {
			if (positionExists(chainId, resnum, icode)) {

				//Position &p = getPosition(chainId,resNumAndIcode);
				Position &p = getLastFoundPosition();
				for (uint i = 0; i < p.getNumberOfIdentities();i++){
					
					if (p.getIdentity(i).atomExists(name)){
						//pAtom = &(p.getIdentity(i)(name));
						pAtom = &(p.getIdentity(i).getLastFoundAtom());
						pAtom->setCoor((*k)->getCoor());
						counter++;
					}
				}
			}
		}

	}
	return counter;
}
		
	


void System::setLinkedPositions(vector<vector<string> > &_linkedPositions){

	for (uint v = 0; v < _linkedPositions.size();v++){

		vector<Position *> positionsLinkedHere;
		for (uint t = 0; t < _linkedPositions[v].size();t++){
			string chain = "";
			int resnum = 0;
			string icode = "";
			bool OK = MslTools::parsePositionId(_linkedPositions[v][t], chain, resnum, icode, 0);
			if (!OK) {
				cerr << "DEPRECATED USE OF UNDERSCORE SEPARATED IDENTIFIERS (I.E. \"A_37\"), USE COMMA SEPARATION (\"A,37\") in void System::setLinkedPositions(vector<vector<string> > &_linkedPositions)" << endl;
				vector<string> tmp = MslTools::tokenizeAndTrim(_linkedPositions[v][t],"_");
				string newId;
				if (tmp.size() > 1) {
					newId = tmp[0] + "," + tmp[1];
				}
				OK = MslTools::parsePositionId(newId, chain, resnum, icode, 0);

			}
			if (OK && positionExists(chain, resnum, icode)){
				positionsLinkedHere.push_back(&getPosition(chain, resnum, icode));
			}else {
				cerr << "ERROR 2222 you wanted position "<< _linkedPositions[v][t] <<" to have some linking properties but it doesn't exist"<<endl;
				exit(2222);
			}
		}
		if (positionsLinkedHere.size() == 0){
			cerr << "ERROR 2222 linkedPosition "<<v<<" had a problem finding a position in the system"<<endl;
			exit(2222);
		}

		// Set first position in list as MASTER
		positionsLinkedHere[0]->setLinkedPositionType(Position::MASTER);
		for (uint p = 1; p < positionsLinkedHere.size();p++){

			// Set all other linked positions at the index v, to SLAVEn
			cout << "Setting position: "<<positionsLinkedHere[p]->getChainId()<<" "<<positionsLinkedHere[p]->getResidueNumber()<<" to SLAVE!\n";

			positionsLinkedHere[p]->setLinkedPositionType(Position::SLAVE);
			positionsLinkedHere[p]->addLinkedPosition(*positionsLinkedHere[0]);
			positionsLinkedHere[0]->addLinkedPosition(*positionsLinkedHere[p]);

			cout << "TYPE: "<<positionsLinkedHere[p]->getLinkedPositionType()<<endl;


		}
		
	}
}

string System::toString() const {
	
	PolymerSequence polSeq;
	stringstream ss;
	polSeq.setSequence(activeAndInactiveAtoms);
	ss << polSeq;
	return ss.str();
}
