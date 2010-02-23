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

#include "Atom.h"
#include "AtomGroup.h"
#include "AtomContainer.h"
#include "IcEntry.h"

using namespace MSL;
using namespace std;



Atom::Atom() : Selectable<Atom>(this) {
	setup(CartesianPoint(0.0, 0.0, 0.0), "", "");
}

Atom::Atom(const string _name, string _element) : Selectable<Atom>(this){
	setup(CartesianPoint(0.0, 0.0, 0.0), _name, _element);
}

Atom::Atom(string _name, Real _x, Real _y, Real _z, string _element) : Selectable<Atom>(this){
	setup(CartesianPoint(_x, _y, _z), _name, _element);
}

Atom::Atom(string _name, const CartesianPoint & _p, string _element) : Selectable<Atom>(this){
	setup(_p, _name, _element);
}

Atom::Atom(const Atom & _atom) : Selectable<Atom>(this){
	pParentGroup = NULL; // note: the copy constructor of a Residue needs to set itself as a parent in the new atoms
	pParentContainer = NULL;
	groupNumber = 0;
	copy(_atom);

	addSelectableFunctions();
}

Atom::~Atom() {
	//removeBonds();
	deletePointers();
}

void Atom::deletePointers() {
	for (vector<CartesianPoint*>::iterator k=pCoorVec.begin(); k != pCoorVec.end(); k++) {
		delete *k;
	}
	pCoorVec.clear();
	currentCoorIterator = pCoorVec.end();
	clearSavedCoor();
}

/*
void Atom::removeBonds() {
	for (map<Atom*, bool>::iterator bondItr=bonds.begin(); bondItr!=bonds.end(); bondItr++) {
		bondItr->first->setBondedTo(this, false);
	}
}
*/

void Atom::clearSavedCoor() {
	for (map<string, CartesianPoint*>::iterator k=savedCoor.begin(); k!=savedCoor.end(); k++) {
		delete k->second;
	}
	savedCoor.clear();
}

void Atom::reset() {
	pParentGroup = NULL;
	pParentContainer = NULL;
	hasCoordinates = true;

	deletePointers();
	name = "";
	residueNumber = 1;
	residueName = "";
	residueIcode = "";
	chainId = "A";
	element = "";
	charge = 0.0;
	radius = -1.0;
	type = "";
	tempFactor = 0.0;
	sasa = 0.0;
	segId = "";
	groupNumber = 0;
	nameSpace = "";
//	bonds.clear();
	oneThreeAtoms.clear();
	oneFourAtoms.clear();
}

void Atom::setup(CartesianPoint _point, string _name, string _element) {
	reset();
	name = _name;
	element = _element;
	pCoorVec.push_back(new CartesianPoint(_point));
	currentCoorIterator = pCoorVec.begin();
	groupNumber = 0;

	addSelectableFunctions();

	// Every atom should be marked as "all"
	setSelectionFlag("all",true);
}

void Atom::operator=(const Atom & _atom) {
	reset();
	copy(_atom);
}

void Atom::copy(const Atom & _atom) {
	// call the functions so that we won't copy the local residue/chain values 
	// but the value of the parent residue and chain if those exists
	//name = _atom.name; 
	setName(_atom.name); // this updates the map in the Residue
	if (pParentGroup == NULL) {
		residueNumber = _atom.getResidueNumber();
		residueName = _atom.getResidueName();
		residueIcode = _atom.getResidueIcode();
		chainId = _atom.getChainId();
		groupNumber = _atom.getGroupNumber();
	}
	element = _atom.element;
	charge = _atom.charge;
	radius = _atom.radius;
	type = _atom.type;
	segId = _atom.segId;
	tempFactor = _atom.tempFactor;
	sasa = _atom.sasa;


	// remove all coordinates and add the new ones
	deletePointers();
	for (vector<CartesianPoint*>::const_iterator k=_atom.pCoorVec.begin(); k!=_atom.pCoorVec.end(); k++) {
		pCoorVec.push_back(new CartesianPoint(**k));
	}
	currentCoorIterator = pCoorVec.begin() + (_atom.currentCoorIterator - _atom.pCoorVec.begin());
	hasCoordinates = _atom.hasCoordinates;

	// Every atom should be marked as "all"
	setSelectionFlag("all",true);

}


/*********************************************************
 *   SETTING THE RESIDUE / CHAIN PROPERTIES FOR THE ATOM 
 *   SETS THEM FOR THE WHOLE RESIDUE / CHAIN (via the Group)
 *  
 *   Is this the obvious behavior?
 *  
 *   Alternatively if it sets it only for the atoms, the
 *   setting is useless if a group exists
 *********************************************************/
void Atom::setResidueName(string _resname) {
	if (pParentGroup != NULL) {
		pParentGroup->setResidueName(_resname);
	} else {
		residueName = _resname;
	}
}

string Atom::getResidueName() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getResidueName();
	} else {
		return residueName;
	}
}

void Atom::setResidueNumber(int _resnum) {
	if (pParentGroup != NULL) {
		pParentGroup->setResidueNumber(_resnum);
	} else {
		residueNumber = _resnum;
		if (pParentContainer != NULL) {
			updateContainerMap();
		}
	}
}

int Atom::getResidueNumber() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getResidueNumber();
	} else {
		return residueNumber;
	}
}

void Atom::setResidueIcode(string _icode) {
	if (pParentGroup != NULL) {
		pParentGroup->setResidueIcode(_icode);
	} else {
		residueIcode = _icode;
		if (pParentContainer != NULL) {
			updateContainerMap();
		}
	}
}

string Atom::getResidueIcode() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getResidueIcode();
	} else {
		return residueIcode;
	}
}

void Atom::setChainId(string _chainId) {
	if (pParentGroup != NULL) {
		pParentGroup->setChainId(_chainId);
	} else {
		chainId = _chainId;
		if (pParentContainer != NULL) {
			updateContainerMap();
		}
	}
}


string Atom::getChainId() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getChainId();
	} else {
		return chainId;
	}
}

// the name space defines if the atoms and residues are named with
// PDB or charmm names (or whatever else)
void Atom::setNameSpace(string _nameSpace) {
	if (pParentGroup != NULL) {
		pParentGroup->setNameSpace(_nameSpace);
	} else {
		nameSpace = _nameSpace;
	}
}

string Atom::getNameSpace() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getNameSpace();
	} else {
		return nameSpace;
	}
}

/*********************************************************/

void Atom::addSelectableFunctions(){
	addStringFunction("NAME", &Atom::getName);
	addStringFunction("RESN", &Atom::getResidueName);
	addIntFunction("RESI", &Atom::getResidueNumber);
	addStringFunction("ICODE", &Atom::getResidueIcode);
	addStringFunction("CHAIN", &Atom::getChainId);


	addRealFunction("X", &Atom::getX);
	addRealFunction("Y", &Atom::getY);
	addRealFunction("Z", &Atom::getZ);

	addBoolFunction("HASCRD", &Atom::hasCoor);

}

// BUILDING FROM IC ========================================================

void Atom::addIcEntry(IcEntry * _ic) {
	// ADD A CHECK THAT _ic a1 has the same address of 'this'
	icEntries.push_back(_ic);
}

/*
bool Atom::buildFromIc(const vector<IcEntry*> & _exclude) {
	if (hasCoordinates) {
		// is already built
		cout << "UUU " << *this << " was ALREADY built" << endl;
		return true;
	}
	// try to build until you find an ic entry that is successfull
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		cout << "UUU trying ic " << *k << endl;
		bool found = false;
		for (vector<IcEntry*>::const_iterator l=_exclude.begin(); l<_exclude.end(); l++) {
			cout << "UUU check exclude " << *k << "!=" << *l << endl;
			if (*k == *l) {
				cout << "UUU found!" << endl;
				found = true;
				break;
			}
		}
		if (found) {
			continue;
		}
		cout << "UUU try IC " << **k << endl;
		if ((*k)->build(this, _exclude)) {
			cout << "UUU " << *this << " we BUILT it" << endl;
			return true;
		}
	}
	cout << "UUU " << *this << " we could NOT build" << endl;
	return false;
}
*/

bool Atom::buildFromIc(bool _onlyFromActive) {
	if (hasCoordinates) {
		// is already built
		return true;
	}
	// try to build until you find an ic entry that is successfull
	map<Atom*, bool> exclude;
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		if ((*k)->build(this, exclude, _onlyFromActive)) {
			return true;
		}
	}
	return false;
}

bool Atom::buildFromIc(map<Atom*, bool> & _exclude, bool _onlyFromActive) {
	if (hasCoordinates) {
		// is already built
		return true;
	}
	// try to build until you find an ic entry that is successfull
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		if (_exclude.find(this) != _exclude.end()) {
			continue;
		}
		if ((*k)->build(this, _exclude, _onlyFromActive)) {
			return true;
		}
	}
	return false;
}

void Atom::removeIcEntry(IcEntry * _ic) {
	/**************************************************************
	 * Remove element if the address of the pointer is the same
	 *
	 * This function is called by the distructor of the IcEntry object
	 * so that the atom removes the pointer to it before it gets
	 * removed from memory
	 **************************************************************/
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		if (*k == _ic) {
			icEntries.erase(k);
		}
	}
}



void Atom::removeAltConformation(unsigned int _i) {
	if (_i>=pCoorVec.size()) {
		return;
	}
	if (pCoorVec.size() > 1) {
		unsigned int current = currentCoorIterator - pCoorVec.begin();
		vector<CartesianPoint*>::iterator k=pCoorVec.begin()+_i;

		// make sure we point to the same coordinate (unless we deleted it)
		if (_i == current) {
			// we deleted the current conformation, point
			// to the first conformation
			currentCoorIterator = pCoorVec.begin();
		} else if (_i > current) {
			currentCoorIterator = pCoorVec.begin() + current;
		} else {
			currentCoorIterator = pCoorVec.begin() + current - 1;
		}

		delete *k;
		pCoorVec.erase(k);
	}

}

bool Atom::getActive() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getActive();
	} else {
		return true;
	}
}

void Atom::updateResidueMap() {
	if (pParentGroup != NULL) {
		pParentGroup->updateResidueMap(this);
	}
}

void Atom::updateContainerMap() {
	if (pParentContainer != NULL) {
		pParentContainer->updateAtomMap(this);
	}
}

unsigned int Atom::getGroupNumber() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getGroupNumber();
	} else {
		return groupNumber;
	}
}

void Atom::setBoundTo(Atom * _pAtom) {
	boundAtoms[_pAtom]; // create 1-2
	_pAtom->boundAtoms[this]; // create reciprocal 1-2

	// to all 1-2 atoms of this atom, add a 1-3 through _pAtom
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k=_pAtom->boundAtoms.begin(); k!=_pAtom->boundAtoms.end(); k++) {
		if (k->first == this) {
			continue;
		}
		if (this->getParentPosition() == k->first->getParentPosition() && this->getParentResidue() != k->first->getParentResidue()) {
			// exclude if the atoms that are on different identities of the same position
			continue;
		}
		boundAtoms[_pAtom][k->first]; // create 1-3
		k->first->boundAtoms[_pAtom][this]; // create reciprocal 1-3
		oneThreeAtoms[k->first][_pAtom] = true;
		k->first->oneThreeAtoms[this][_pAtom] = true;


		// to all 1-2 atoms of the 1-2 atom, add a 1-4 through _pAtom and its 1-2 atom
		for (map<Atom*, map<Atom*, bool> >::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			if (l->first == this || l->first == _pAtom) {
				continue;
			}
			if (this->getParentPosition() == l->first->getParentPosition() && this->getParentResidue() != l->first->getParentResidue()) {
				// exclude if the atoms that are on different identities of the same position
				continue;
			}
			boundAtoms[_pAtom][k->first][l->first] = true; // create 1-4
			l->first->boundAtoms[k->first][_pAtom][this] = true; // create reciprocal 1-4
			oneFourAtoms[l->first][_pAtom][k->first] = true;
			l->first->oneFourAtoms[this][k->first][_pAtom] = true;
		}
	}
	// to all 1-2 atoms of _pAtom, add a 1-3 through this atom
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		if (k->first == _pAtom) {
			continue;
		}
		if (_pAtom->getParentPosition() == k->first->getParentPosition() && _pAtom->getParentResidue() != k->first->getParentResidue()) {
			// exclude if the atoms that are on different identities of the same position
			continue;
		}
		_pAtom->boundAtoms[this][k->first]; // create 1-3
		k->first->boundAtoms[this][_pAtom]; // create reciprocal 1-3
		_pAtom->oneThreeAtoms[k->first][this] = true;
		k->first->oneThreeAtoms[_pAtom][this] = true;


		// to all 1-2 atoms of the 1-2 atom, add a 1-4 through _pAtom and its 1-2 atom
		for (map<Atom*, map<Atom*, bool> >::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			if (l->first == this || l->first == _pAtom) {
				continue;
			}
			if (_pAtom->getParentPosition() == l->first->getParentPosition() && _pAtom->getParentResidue() != l->first->getParentResidue()) {
				// exclude if the atoms that are on different identities of the same position
				continue;
			}
			_pAtom->boundAtoms[this][k->first][l->first] = true; // create 1-4
			l->first->boundAtoms[k->first][this][_pAtom] = true; // create reciprocal 1-4
			_pAtom->oneFourAtoms[l->first][this][k->first] = true;
			l->first->oneFourAtoms[_pAtom][k->first][this] = true;
		}
	}
	
	// finally, combinatorially add 1-4 relationships between all 1-2 atoms of this atom and _pAtom
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		if (k->first == _pAtom) {
			continue;
		}
		if (_pAtom->getParentPosition() == k->first->getParentPosition() && _pAtom->getParentResidue() != k->first->getParentResidue()) {
			// exclude if the atoms that are on different identities of the same position
			continue;
		}
		for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator l=_pAtom->boundAtoms.begin(); l!=_pAtom->boundAtoms.end(); l++) {
			if (l->first == this || l->first == k->first) {
				continue;
			}
			if (this->getParentPosition() == l->first->getParentPosition() && this->getParentResidue() != l->first->getParentResidue()) {
				// exclude if the atoms that are on different identities of the same position
				continue;
			}
			if (k->first->getParentPosition() == l->first->getParentPosition() && k->first->getParentResidue() != l->first->getParentResidue()) {
				// exclude if the atoms that are on different identities of the same position
				continue;
			}
			k->first->boundAtoms[this][_pAtom][l->first] = true; // create 1-4 beetween the 1-2 or this and the 1-2 of _pAtom
			l->first->boundAtoms[_pAtom][this][k->first] = true; // create reciprocal 1-4
			k->first->oneFourAtoms[l->first][this][_pAtom] = true;
			l->first->oneFourAtoms[k->first][_pAtom][this] = true;
		}
	}
}

Residue * Atom::getParentResidue() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getParentResidue();
	} else {
		return NULL;
	}
}

Position * Atom::getParentPosition() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getParentPosition();
	} else {
		return NULL;
	}
}

Chain * Atom::getParentChain() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getParentChain();
	} else {
		return NULL;
	}
}

System * Atom::getParentSystem() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getParentSystem();
	} else {
		return NULL;
	}
}

unsigned int Atom::getIdentityIndex() {
	if (pParentGroup != NULL) {
		return pParentGroup->getIdentityIndex();
	} else {
		return 0;
	}
}
CartesianPoint Atom::getGroupGeometricCenter(unsigned int _stamp) {
	if (pParentGroup != NULL) {
		return pParentGroup->getGeometricCenter(_stamp);
	} else {
		return *(*currentCoorIterator);
	}
}

set<Atom*> Atom::findLinkedAtoms(const set<Atom*> & _excluded) {
	// Find all atoms that are connected to this atom, except those going through the exclusion list.
	// This may for example find the end of a side chain or half of the protein

//	// get all the atoms that are directly bonded to _pAtom
//	map<Atom*, bool> bonded = _pAtom->getBonds();

	set<Atom*> linked;
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		if (_excluded.find(k->first) != _excluded.end()) {
			continue;
		}
		linked.insert(k->first);
		set<Atom*> localExcluded = _excluded;
		localExcluded.insert(this);
		set<Atom*> linkedToBonded = k->first->findLinkedAtoms(localExcluded);
		linked.insert(linkedToBonded.begin(), linkedToBonded.end());
	}
	return linked;
/*
	for (map<Atom*, bool>::iterator k=bonded.begin(); k!=bonded.end(); k++) {
		if (_excluded.find(k->first) != _excluded.end() || _list.find(k->first) != _list.end()) {
			// already in the list or excluded atoms
			continue;
		} else {
			_list[k->first] = true;
			// add the atoms bonded to this atom to the list recursively
			findLinkedAtoms(k->first, _excluded, _list);
		}
	}
*/
}

