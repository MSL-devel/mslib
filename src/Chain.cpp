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

#include "Chain.h"
#include "System.h"

Chain::Chain() {
	setup("A");
}

Chain::Chain(string _chainId) {
	setup(_chainId);
}

Chain::Chain(const vector<Residue> & _residues, string _chainId) {
	setup(_chainId);
	addResidues(_residues);

}

Chain::Chain(const AtomVector & _atoms, string _chainId) {
	setup(_chainId);
	addAtoms(_atoms);

}

Chain::Chain(const Chain & _chain) {
	pParentSystem = NULL;
	noUpdateIndex_flag = false;
	copy(_chain);
}

Chain::~Chain() {
	deletePointers();
}

void Chain::operator=(const Chain & _chain) {
	copy(_chain);
}

void Chain::setup(string _chainId) {
	pParentSystem = NULL;
	chainId = _chainId;
	noUpdateIndex_flag = false;

}

void Chain::copy(const Chain & _chain) {
	deletePointers();
	chainId = _chain.chainId;
	for (vector<Position*>::const_iterator k=_chain.positions.begin(); k!=_chain.positions.end(); k++) {
		positions.push_back(new Position(**k));
		positions.back()->setParentChain(this);
		char c [1000];
		sprintf(c, "%06d%1s", positions.back()->getResidueNumber(), positions.back()->getResidueIcode().c_str());
		positionMap[positions.back()->getResidueNumber()][positions.back()->getResidueIcode()] = positions.back();
	}
	updateIndexing();
	updateAllAtomIndexing();
	updateSystemMap();
}

void Chain::deletePointers() {
	positionMap.clear();
	activeAtoms.clear();
	activeAndInactiveAtoms.clear();
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		delete *k;
	}
	positions.clear();
	updateSystemActiveAtomList();
	updateSystemAllAtomList();
}

void Chain::addResidue(AtomVector _atoms, string _name) {
	int resNum = 1;
	if (positions.size() > 0) {
		resNum = positions.back()->getResidueNumber() + 1;
	}
	addResidue(_atoms, _name, resNum, "");
}

void Chain::addResidue(AtomVector _atoms, string _name, unsigned int _resNum, string _iCode) {
	addResidue(Residue(_atoms, _name, _resNum, _iCode), _resNum, _iCode);
}

void Chain::addResidue(const Residue & _residue) {
	addResidue(_residue, _residue.getResidueNumber(), _residue.getResidueIcode());
}

void Chain::addResidue(const Residue & _residue, unsigned int _resNum, string _iCode) {
	// this function is called add residue but it really adds a position with
	// one identity (the Residue _residue)
	positions.push_back(new Position(_residue, _resNum, _iCode));
	positions.back()->setChainId(chainId);
	positions.back()->setParentChain(this);
	char c [1000];
	sprintf(c, "%06d%1s", _resNum, _iCode.c_str());
	positionMap[_resNum][_iCode] = positions.back();

	// update the list of active atoms
	activeAtoms.insert(activeAtoms.end(), positions.back()->getAtoms().begin(), positions.back()->getAtoms().end());
	activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions.back()->getAllAtoms().begin(), positions.back()->getAllAtoms().end());

	updateSystemActiveAtomList();
	updateSystemAllAtomList();
}

void Chain::addResidues(const vector<Residue> & _residues) {
	for (vector<Residue>::const_iterator k=_residues.begin(); k!=_residues.end(); k++) {
		addResidue(*k);
	}
}

bool Chain::removeResidue(int _resNum, string _iCode) {
	map<int, map<string, Position*> >::iterator foundPosition=positionMap.find(_resNum);
	if (foundPosition!=positionMap.end()) {
		map<string, Position*>::iterator foundPosition2=foundPosition->second.find(_iCode);
		if (foundPosition2!=foundPosition->second.end()) {

			// erase from the atoms list
			for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
				if ((*foundPosition2).second == *k) {
				
					// deallocate from memory and remove from list
					delete *k;
					positions.erase(k);
					foundPosition->second.erase(foundPosition2);
					if (foundPosition->second.size() == 0) {
						positionMap.erase(foundPosition);
					}
					updateIndexing(); // update the active atoms	
					updateAllAtomIndexing();
					return true;
				}
			}
		}
	}
	return false;
}

bool Chain::addIdentityToPosition(AtomVector _atoms, string _name, unsigned int _resNum, string _iCode) {
	return addIdentityToPosition(Residue(_atoms, _name, _resNum, _iCode), _resNum, _iCode);
}

bool Chain::addIdentityToPosition(const Residue & _residue, unsigned int _resNum, string _iCode) {
	map<int, map<string, Position*> >::iterator foundPosition=positionMap.find(_resNum);
	if (foundPosition==positionMap.end()) {
		return false;
	} else {
		map<string, Position*>::iterator foundPosition2=foundPosition->second.find(_iCode);
		if (foundPosition2==foundPosition->second.end()) {
			return false;
		} else {
			positionMap[_resNum][_iCode]->addIdentity(_residue);
			return true;
		}
	}
}

void Chain::removeAllResidues() {
	deletePointers();
}

void Chain::addAtoms(const AtomVector & _atoms) {

	/***********************************************
	 *  This function splits the AtomVector in a number
	 *  of AtomVector objects by residue number and insertion
	 *  code and then calls the addAtoms(const AtomVector & _atoms)
	 *  function of the positions to take care of the rest
	 *
	 *  Note: the chain Id of the atoms are ignored
	 ***********************************************/

	noUpdateIndex_flag = true;

	// store the order of the positions so that it will be preserved
	vector<int> resNumOrder;
	vector<string> iCodeOrder;
	map<int, map<string, AtomVector> > dividedInPositions2;
	for (AtomVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		if (dividedInPositions2.find((*k)->getResidueNumber()) == dividedInPositions2.end() || dividedInPositions2[(*k)->getResidueNumber()].find((*k)->getResidueIcode()) == dividedInPositions2[(*k)->getResidueNumber()].end()) {
			resNumOrder.push_back((*k)->getResidueNumber());
			iCodeOrder.push_back((*k)->getResidueIcode());
		}
		dividedInPositions2[(*k)->getResidueNumber()][(*k)->getResidueIcode()].push_back(*k);
	}


	/*
	for (map<int, map<string, AtomVector> >::iterator k=dividedInPositions2.begin(); k!=dividedInPositions2.end(); k++) {
		for (map<string, AtomVector>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			map<int, map<string, Position*> >::iterator foundPosition=positionMap.find(k->first);
			bool found = false;
			if (foundPosition!=positionMap.end()) {
				map<string, Position*>::iterator foundPosition2=foundPosition->second.find(l->first);
				if (foundPosition2!=foundPosition->second.end()) {
					/ ***********************************************
					 *  A position with the resnum/icode already EXISTS, add
					 *  the atoms to it
					 *********************************************** /
					(*foundPosition2).second->addAtoms(l->second);
					found = true;
				}
			}
			if (!found) {
				/ ***********************************************
				 *  A position with the resnum/icode DOES NOT EXIST, 
				 *  create a new position first and add
				 *  the atoms to it
				 *
				 *  l = iterator, pointer to element of map
				 *  *l = element of map
				 *  l->second  = an AtomVector
				 *  *(l->second.begin()) = Atom *
				 *
				 *********************************************** /
				Atom * tmpAtom = *(l->second.begin());
				positions.push_back(new Position( tmpAtom->getResidueNumber(), tmpAtom->getResidueIcode()));
				positions.back()->setParentChain(this);
				positionMap[tmpAtom->getResidueNumber()][tmpAtom->getResidueIcode()] = positions.back();
				positions.back()->addAtoms(l->second);

				AtomVector active =  positions.back()->getAtoms();
				activeAtoms.insert(activeAtoms.end(), positions.back()->getAtoms().begin(), positions.back()->getAtoms().end());
				activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions.back()->getAllAtoms().begin(), positions.back()->getAllAtoms().end());

			}
		}
	}
	*/
	for (unsigned int i=0; i<resNumOrder.size(); i++) {
		map<int, map<string, Position*> >::iterator foundPosition=positionMap.find(resNumOrder[i]);
		bool found = false;
		if (foundPosition!=positionMap.end()) {
			map<string, Position*>::iterator foundPosition2=foundPosition->second.find(iCodeOrder[i]);
			if (foundPosition2!=foundPosition->second.end()) {
				/***********************************************
				 *  A position with the resnum/icode already EXISTS, add
				 *  the atoms to it
				 ***********************************************/
				(*foundPosition2).second->addAtoms(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]]);
				found = true;
			}
		}
		if (!found) {
			/***********************************************
			 *  A position with the resnum/icode DOES NOT EXIST, 
			 *  create a new position first and add
			 *  the atoms to it
			 *
			 *  l = iterator, pointer to element of map
			 *  *l = element of map
			 *  l->second  = an AtomVector
			 *  *(l->second.begin()) = Atom *
			 *
			 ***********************************************/
			Atom * tmpAtom = *(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]].begin());
			/**********************************************************
			 * changed to insert residues according to the order of their
			 * residue number (and icode)
			 **********************************************************/
			int tmpResnum = tmpAtom->getResidueNumber();
			string tmpIcode = tmpAtom->getResidueIcode();
			vector<Position*>::iterator insertPosition = positions.end();
			unsigned int index = positions.size();
			if (index != 0 && MslTools::sortByResnumIcodeAscending(tmpResnum, tmpIcode, positions.back()->getResidueNumber(), positions.back()->getResidueIcode())) {
				// positions is not empty and the resnum and icode do not place the residue after the last
				// position: search for appropriate placement
				for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
					if (MslTools::sortByResnumIcodeAscending(tmpResnum, tmpIcode, (*k)->getResidueNumber(), (*k)->getResidueIcode())) {
						insertPosition = k;
						index = insertPosition - positions.begin();
						//cout << "UUU " << tmpResnum << " " << tmpIcode << " is before " << (*k)->getResidueNumber() << " " << (*k)->getResidueIcode() << endl;
						break;
					}
				}
			//} else {
			//	cout << "UUU placed " << tmpResnum << " " << tmpIcode << " position last" << endl;
			}
			positions.insert(insertPosition, new Position( tmpResnum, tmpIcode));
		//	cout << "  UUU Inserted " << tmpResnum << " " << tmpIcode << " at " << index << " of " << positions.size() -1 << endl;
			positions[index]->setParentChain(this);
			positionMap[tmpAtom->getResidueNumber()][tmpAtom->getResidueIcode()] = positions[index];
			positions[index]->addAtoms(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]]);
			AtomVector active = positions[index]->getAtoms();
			activeAtoms.insert(activeAtoms.end(), positions[index]->getAtoms().begin(), positions[index]->getAtoms().end());
			activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions[index]->getAllAtoms().begin(), positions[index]->getAllAtoms().end());
			/**********************************************************
			 * OLD CODE, before the residues were added at the end
			 * no matter what their resnumber was
			 **********************************************************
			positions.push_back(new Position( tmpAtom->getResidueNumber(), tmpAtom->getResidueIcode()));
			positions.back()->setParentChain(this);
			positionMap[tmpAtom->getResidueNumber()][tmpAtom->getResidueIcode()] = positions.back();
			positions.back()->addAtoms(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]]);

			AtomVector active =  positions.back()->getAtoms();
			activeAtoms.insert(activeAtoms.end(), positions.back()->getAtoms().begin(), positions.back()->getAtoms().end());
			activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions.back()->getAllAtoms().begin(), positions.back()->getAllAtoms().end());
			 **********************************************************/

		}
	}

	noUpdateIndex_flag = false;
	updateIndexing();
	updateAllAtomIndexing();

}


// the name space defines if the atoms and residues are named with
// PDB or charmm names (or whatever else)
void Chain::setNameSpace(string _nameSpace) {
	if (pParentSystem != NULL) {
		pParentSystem->setNameSpace(_nameSpace);
	} else {
		nameSpace = _nameSpace;
	}
}

string Chain::getNameSpace() const {
	if (pParentSystem != NULL) {
		return pParentSystem->getNameSpace();
	} else {
		return nameSpace;
	}
}

void Chain::updatePositionMap(Position * _position) {
	// if a residue changes its name it calls this function
	// to update the map
	for (map<int, map<string, Position*> >::iterator k=positionMap.begin(); k!=positionMap.end(); k++) {
		for (map<string, Position*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			if (l->second == _position) {
				if (k->first != _position->getResidueNumber() || l->first != _position->getResidueIcode()) {
					k->second.erase(l);
					positionMap.erase(k);
					positionMap[_position->getResidueNumber()][_position->getResidueIcode()] = _position;
				}
			}
		}
	}
}

void Chain::updateSystemMap(){
	if (pParentSystem != NULL) {
		pParentSystem->updateChainMap(this);
	}
}

void Chain::updateSystemActiveAtomList() {
	/******************************************************
	 * this should be called if
	 *   - there are atom changes in the active residue (addition or removal)
	 *   - the active residue changes (i.e. LEU to VAL)
	 *
	 * this is currently a dumb implementation:
	 * the chain rebuilds the index from scratch
	 * - TODO: find a faster way to just swap the atoms as
	 *   needed
	 ******************************************************/
	if (pParentSystem != NULL) {
		pParentSystem->updateIndexing();
	}
}

void Chain::updateSystemAllAtomList() {
	/******************************************************
	 * this should be called if
	 *   - there are atom changes in the number of total atoms, active or
	 *     inactive
	 *   - NOT FOR active residue changes (i.e. LEU to VAL)
	 *
	 * this is currently a dumb implementation:
	 * the chain rebuilds the index from scratch
	 * - TODO: find a faster way to just swap the atoms as
	 *   needed
	 ******************************************************/
	if (pParentSystem != NULL) {
		pParentSystem->updateAllAtomIndexing();
	}
}

/*
void Chain::swapInActiveList(Position * _position, AtomVector & _atoms) {
	
	// NEW METHOD
	map<Position *, ResidueAtoms>::iterator found = residueLookupMap.find(_position);
	if (found != residueLookupMap.end()) {
		if (found->second.size > 0) {
			activeAtoms.erase(activeAtoms.begin() + found->second.start, activeAtoms.begin() + found->second.start + found->second.size - 1);
		}
		if (_atoms.size() > 0) {
			activeAtoms.insert(activeAtoms.begin() + found->second.start, _atoms.begin(), _atoms.end());
		}
	}
	if (pParentSystem != NULL) {
		swapInActiveList(_position, _atoms);
	}
	
}
*/

void Chain::updateIndexing() {
	if (noUpdateIndex_flag) {
		return;
	}
	activeAtoms.clear();
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		activeAtoms.insert(activeAtoms.end(), (*k)->getAtoms().begin(), (*k)->getAtoms().end());

	}

	updateSystemActiveAtomList();
}

void Chain::updateAllAtomIndexing() {
	if (noUpdateIndex_flag) {
		return;
	}
	activeAndInactiveAtoms.clear();
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), (*k)->getAllAtoms().begin(), (*k)->getAllAtoms().end());
	}

	updateSystemAllAtomList();
}

bool Chain::exists(int _resNum) {
	map<int, map<string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition = found->second.find("");
		return foundPosition != found->second.end();
	} else {
		return false;
	}
}

bool Chain::exists(string _resNumAndIcode) {
	int resNum;
	string iCode;
	MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode);
	map<int, map<string, Position*> >::iterator found=positionMap.find(resNum);
	if (found != positionMap.end()) {
		foundPosition = found->second.find(iCode);
		return foundPosition != found->second.end();
	} else {
		return false;
	}
}
bool Chain::exists(int _resNum, string _name) {
	map<int, map<string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find("");
		if (foundPosition != found->second.end()) {
			return foundPosition->second->exists(_name);
		}
	}
	return false;
}
bool Chain::exists(string _resNumAndIcode, string _name) {
	int resNum=0;
	string iCode;
	MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode);
	map<int, map<string, Position*> >::iterator found=positionMap.find(resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find(iCode);
		if (foundPosition != found->second.end()) {
			return foundPosition->second->exists(_name);
		}
	}
	return false;
}
 
bool Chain::exists(int _resNum, string _name, string _identity) {
	map<int, map<string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find("");
		if (foundPosition != found->second.end()) {
			return foundPosition->second->exists(_name, _identity);
		}
	}
	return false;
}
bool Chain::exists(string _resNumAndIcode, string _name, string _identity) {
	int resNum=0;
	string iCode;
	MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode);
	map<int, map<string, Position*> >::iterator found=positionMap.find(resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find(iCode);
		if (foundPosition != found->second.end()) {
			return foundPosition->second->exists(_name, _identity);
		}
	}
	return false;
}

unsigned int Chain::getPositionIndex(const Position * _pPos) const {
	if (pParentSystem != NULL) {
		return pParentSystem->getPositionIndex(_pPos);
	} else {
		for (vector<Position*>::const_iterator k=positions.begin(); k!=positions.end(); k++) {
			if (_pPos == *k) {
				return k-positions.begin();
			}
		}
	}
	cerr << "ERROR 34193: Position not found in unsigned int Chain::getPositionIndex(Position * _pPos)" << endl;
	exit(34193);
}

