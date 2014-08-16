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


#include "Chain.h"
#include "System.h"
#include "PolymerSequence.h"

using namespace MSL;
using namespace std;


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

Chain::Chain(const AtomPointerVector & _atoms, string _chainId) {
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
	//	char c [1000];
	//	sprintf(c, "%06d%1s", positions.back()->getResidueNumber(), positions.back()->getResidueIcode().c_str());
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
	// avoid calling updates for deleted positions
	noUpdateIndex_flag = true;
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		delete *k;
		*k = NULL;
	}
	positions.clear();
	// reset the flag and run the updates
	noUpdateIndex_flag = false;
	updateIndexing();
	updateAllAtomIndexing();
	updateSystemActiveAtomList();
	updateSystemAllAtomList();
}

void Chain::addResidue(AtomPointerVector _atoms, string _name) {
	int resNum = 1;
	if (positions.size() > 0) {
		resNum = positions.back()->getResidueNumber() + 1;
	}
	addResidue(_atoms, _name, resNum, "");
}

void Chain::addResidue(AtomPointerVector _atoms, string _name, unsigned int _resNum, string _iCode) {
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
	activeAtoms.insert(activeAtoms.end(), positions.back()->getAtomPointers().begin(), positions.back()->getAtomPointers().end());
	activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions.back()->getAllAtomPointers().begin(), positions.back()->getAllAtomPointers().end());

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

bool Chain::addIdentityToPosition(AtomPointerVector _atoms, string _name, unsigned int _resNum, string _iCode) {
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

void Chain::addAtoms(const AtomPointerVector & _atoms, bool _keepOrder) {

	/***********************************************
	 *  This function splits the AtomPointerVector in a number
	 *  of AtomPointerVector objects by residue number and insertion
	 *  code and then calls the addAtoms(const AtomPointerVector & _atoms)
	 *  function of the positions to take care of the rest
	 *
	 *  Note: the chain Id of the atoms are ignored
	 ***********************************************/

	noUpdateIndex_flag = true;

	// store the order of the positions so that it will be preserved
	vector<int> resNumOrder;
	vector<string> iCodeOrder;
	map<int, map<string, AtomPointerVector> > dividedInPositions2;
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		if (dividedInPositions2.find((*k)->getResidueNumber()) == dividedInPositions2.end() || dividedInPositions2[(*k)->getResidueNumber()].find((*k)->getResidueIcode()) == dividedInPositions2[(*k)->getResidueNumber()].end()) {
			resNumOrder.push_back((*k)->getResidueNumber());
			iCodeOrder.push_back((*k)->getResidueIcode());
		}
		dividedInPositions2[(*k)->getResidueNumber()][(*k)->getResidueIcode()].push_back(*k);
	}

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
			 *  l->second  = an AtomPointerVector
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
			if (index != 0 && !_keepOrder && MslTools::sortByResnumIcodeAscending(tmpResnum, tmpIcode, positions.back()->getResidueNumber(), positions.back()->getResidueIcode())) {
				// new position does not go at the end: search for appropriate placement
				for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
					if (MslTools::sortByResnumIcodeAscending(tmpResnum, tmpIcode, (*k)->getResidueNumber(), (*k)->getResidueIcode())) {
						insertPosition = k;
						index = insertPosition - positions.begin();
						break;
					}
				}
			}
			positions.insert(insertPosition, new Position( tmpResnum, tmpIcode));
			positions[index]->setParentChain(this);
			positionMap[tmpAtom->getResidueNumber()][tmpAtom->getResidueIcode()] = positions[index];
			positions[index]->addAtoms(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]]);
			AtomPointerVector active = positions[index]->getAtomPointers();
			activeAtoms.insert(activeAtoms.end(), positions[index]->getAtomPointers().begin(), positions[index]->getAtomPointers().end());
			activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions[index]->getAllAtomPointers().begin(), positions[index]->getAllAtomPointers().end());
			/**********************************************************
			 * OLD CODE, before the residues were added at the end
			 * no matter what their resnumber was
			 **********************************************************
			positions.push_back(new Position( tmpAtom->getResidueNumber(), tmpAtom->getResidueIcode()));
			positions.back()->setParentChain(this);
			positionMap[tmpAtom->getResidueNumber()][tmpAtom->getResidueIcode()] = positions.back();
			positions.back()->addAtoms(dividedInPositions2[resNumOrder[i]][iCodeOrder[i]]);

			AtomPointerVector active =  positions.back()->getAtomPointers();
			activeAtoms.insert(activeAtoms.end(), positions.back()->getAtomPointers().begin(), positions.back()->getAtomPointers().end());
			activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), positions.back()->getAllAtomPointers().begin(), positions.back()->getAllAtomPointers().end());
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
//					positionMap.erase(k);
					if (k->second.size() == 0) {
						positionMap.erase(k);
					}
					positionMap[_position->getResidueNumber()][_position->getResidueIcode()] = _position;
					return;
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

void Chain::renumberChain(int _start) {
	for (unsigned int i=0; i<positions.size(); i++) {
		positions[i]->renumberNoUpdate(_start + i);
	}
	positionMap.clear();
	for (std::vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		//positionMap[positions.back()->getResidueNumber()][positions.back()->getResidueIcode()] = *k;
		positionMap[(*k)->getResidueNumber()][(*k)->getResidueIcode()] = *k;
	}
	updateSystemMap();
}

/*
void Chain::swapInActiveList(Position * _position, AtomPointerVector & _atoms) {
	
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
		activeAtoms.insert(activeAtoms.end(), (*k)->getAtomPointers().begin(), (*k)->getAtomPointers().end());

	}

	updateSystemActiveAtomList();
}

void Chain::updateAllAtomIndexing() {
	if (noUpdateIndex_flag) {
		return;
	}
	activeAndInactiveAtoms.clear();
	for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), (*k)->getAllAtomPointers().begin(), (*k)->getAllAtomPointers().end());
	}

	updateSystemAllAtomList();
}

/*
bool Chain::exists(int _resNum) {
	cerr << "DEPRECATED: bool Chain::exists(int _resNum), use Chain::positionExist(string _positionId)" << endl;
	map<int, map<string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition = found->second.find("");
		return foundPosition != found->second.end();
	} else {
		return false;
	}
}

bool Chain::exists(string _resNumAndIcode) {
	cerr << "DEPRECATED: bool Chain::exists(string _resNumAndIcode), use Chain::positionExist(string _positionId)" << endl;
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
	cerr << "DEPRECATED: bool Chain::exists(int _resNum, string _name), use Chain::atomExist(string _atomId)" << endl;
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
	cerr << "DEPRECATED: bool Chain::exists(string _resNumAndIcode, string _name), use Chain::atomExist(string _atomId)" << endl;
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
	cerr << "DEPRECATED: bool Chain::exists(int _resNum, string _name, string _identity), use Chain::identityExist(string _identityId)" << endl;
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
	cerr << "DEPRECATED: bool Chain::exists(string _resNumAndIcode, string _name, string _identity), use Chain::identityExist(string _identityId)" << endl;
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
*/

int Chain::getPositionIndex(const Position * _pPos) const {

	for (vector<Position*>::const_iterator k=positions.begin(); k!=positions.end(); k++) {
	  if (_pPos == *k) {
	    return k-positions.begin();
	  }
	}

	cerr << "ERROR 34193: Position not found in unsigned int Chain::getPositionIndex(Position * _pPos)" << endl;
	exit(34193);
}

int Chain::getReversePositionIndex(const Position * _pPos) const {
	if (pParentSystem != NULL){
		return (getPositionIndex(_pPos) - positionSize());
	}

	cerr << "ERROR 34194: Position not found in unsigned int Chain::getReversePositionIndexInSystem(Position * _pPos)" << endl;
	exit(34194);
}

int Chain::getPositionIndexInSystem(const Position * _pPos) const {

	if (pParentSystem != NULL) {
		return pParentSystem->getPositionIndex(_pPos);
	} 

	cerr << "ERROR 34195: Position not found in unsigned int Chain::getPositionIndexInSystem(Position * _pPos)" << endl;
	exit(34195);
}


int Chain::getReversePositionIndexInSystem(const Position * _pPos) const {
  if (pParentSystem != NULL){
    return (getPositionIndexInSystem(_pPos) - pParentSystem->positionSize());
  }

  cerr << "ERROR 34196: Position not found in unsigned int Chain::getReversePositionIndexInSystem(Position * _pPos)" << endl;
  exit(34196);
}


string Chain::toString() const {
	
	PolymerSequence polSeq;
	stringstream ss;
	polSeq.setSequence(activeAndInactiveAtoms);
	ss << polSeq;
	return ss.str();
}
