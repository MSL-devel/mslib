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



#include "Atom.h"
#include "AtomGroup.h"
#include "AtomContainer.h"
#include "IcEntry.h"
#include "IcTable.h"

using namespace MSL;
using namespace std;


#include "MslOut.h"
static MslOut MSLOUT("Atom");

Atom::Atom() : Selectable<Atom>(this) {
	setup(CartesianPoint(0.0, 0.0, 0.0), "", "");
}

Atom::Atom(const string _atomId, string _element) : Selectable<Atom>(this){
	setup(CartesianPoint(0.0, 0.0, 0.0), _atomId, _element);
}

Atom::Atom(string _atomId, Real _x, Real _y, Real _z, string _element) : Selectable<Atom>(this){
	setup(CartesianPoint(_x, _y, _z), _atomId, _element);
}

Atom::Atom(string _atomId, const CartesianPoint & _p, string _element) : Selectable<Atom>(this){
	setup(_p, _atomId, _element);
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
	setUnboundFromAll(); // this will also remove the bonds to this atom from the atoms bound to it
	removeFromIc(); // call the IcTable (if present) and make sure this atom is removed
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

void Atom::clearSavedCoor(string _coordName) {
	if (_coordName == "") {
		// name left blank, erase all
		for (map<string, CartesianPoint*>::iterator k=savedCoor.begin(); k!=savedCoor.end(); k++) {
			delete k->second;
		}
		savedCoor.clear();
		for (map<string, vector<CartesianPoint*> >::iterator k=savedAltCoor.begin(); k!=savedAltCoor.end(); k++) {
			for (vector<CartesianPoint*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
				delete *l;
			}
		}
		savedAltCoor.clear();
		for (map<string, vector<CartesianPoint*> >::iterator k=savedHiddenCoor.begin(); k!=savedHiddenCoor.end(); k++) {
			for (vector<CartesianPoint*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
				delete *l;
			}
		}
		savedAltCoorCurrent.clear();
		savedHiddenCoor.clear();
		savedHiddenCoorIndeces.clear();
	} else {
		// name given, erase only specific entry
		map<string, CartesianPoint*>::iterator f1 = savedCoor.find(_coordName);
		if (f1 != savedCoor.end()) {
			delete f1->second;
			savedCoor.erase(f1);
		}
		map<string, vector<CartesianPoint*> >::iterator f2 = savedAltCoor.find(_coordName);
		if (f2 != savedAltCoor.end()) {
			for (vector<CartesianPoint*>::iterator l=f2->second.begin(); l!=f2->second.end(); l++) {
				delete *l;
			}
			savedAltCoor.erase(f2);
			map<string, unsigned int>::iterator f3 = savedAltCoorCurrent.find(_coordName);
			savedAltCoorCurrent.erase(f3);
		}
		f2 = savedHiddenCoor.find(_coordName);
		if (f2 != savedHiddenCoor.end()) {
			for (vector<CartesianPoint*>::iterator l=f2->second.begin(); l!=f2->second.end(); l++) {
				delete *l;
			}
			savedHiddenCoor.erase(f2);
			std::map<std::string, std::vector<unsigned int> >::iterator f4 = savedHiddenCoorIndeces.find(_coordName);
			savedHiddenCoorIndeces.erase(f4);
		}
	}
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

void Atom::setup(CartesianPoint _point, string _atomId, string _element) {
	reset();
	
	string chain = "";
	int resnum = 1;
	string icode = "";
	string identity = "";
	string atomname = "";
	if (MslTools::parseAtomId(_atomId, chain, resnum, icode, atomname, 2)) {
		chainId = chain;
		residueNumber = resnum;
		residueIcode = icode;
		name = atomname;
	} else if (MslTools::parseAtomOfIdentityId(_atomId, chain, resnum, icode, identity, atomname, 3)) {
		chainId = chain;
		residueNumber = resnum;
		residueIcode = icode;
		residueName = identity;
		name = atomname;
	} else {
		name = _atomId;
	}
	element = _element;
	pCoorVec.push_back(new CartesianPoint(_point));
	currentCoorIterator = pCoorVec.begin();
	groupNumber = 0;

	addSelectableFunctions();

	// Every atom should be marked as "all"
	setSelectionFlag("all",true);

	minIndex = -1;
	toStringFormat = 0;
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

	minIndex = _atom.minIndex;
	toStringFormat = _atom.toStringFormat;
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
	//map<Atom*, bool> exclude;
	map<IcEntry*, bool> exclude;
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		//cout << "   UUUA " << getAtomId() << " try IC entry " << **k << endl;
		if ((*k)->build(this, exclude, _onlyFromActive)) {
			//cout << "   UUUA succesful" << endl;
			return true;
		}
		//cout << "   UUUA NOT succesful" << endl;
	}
	return false;
}

/*
bool Atom::buildFromIc(map<Atom*, bool> & _exclude, bool _onlyFromActive) {
	if (hasCoordinates) {
		// is already built
		return true;
	}
	// try to build until you find an ic entry that is successfull
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		cout << "   UUUAe " << getAtomId() << " try IC entry " << **k << endl;
		if (_exclude.find(this) != _exclude.end()) {
			cout << "   UUUAe excluded" << endl;
			continue;
		}
		if ((*k)->build(this, _exclude, _onlyFromActive)) {
			cout << "   UUUAe succesful" << endl;
			return true;
		}
		cout << "   UUUAe NOT succesful" << endl;
	}
	return false;
}
*/

bool Atom::buildFromIc(map<IcEntry*, bool> & _exclude, bool _onlyFromActive) {
	if (hasCoordinates) {
		// is already built
		return true;
	}
	// try to build until you find an ic entry that is successfull
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		//cout << "   UUUAe " << getAtomId() << " try IC entry " << **k << endl;
		if (_exclude.find(*k) != _exclude.end()) {
			//cout << "   UUUAe excluded" << endl;
			continue;
		}
		if ((*k)->build(this, _exclude, _onlyFromActive)) {
			//cout << "   UUUAe succesful" << endl;
			return true;
		}
		//cout << "   UUUAe NOT succesful" << endl;
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
			return;
		}
	}
}

void Atom::printIcEntries() const {
	for (vector<IcEntry*>::const_iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		cout << **k << endl;
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
bool Atom::getHidden() const {
	if (pParentGroup != NULL) {
		return pParentGroup->getHidden();
	} else {
		return false;
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

void Atom::setUnboundFromAll(bool _propagate) {
	// return if there are no bonds
	if (boundAtoms.size() == 0) {
		return;
	}
	if (_propagate) {
		// remove the bond to this atoms from all bound atoms
		for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k = boundAtoms.begin(); k!=boundAtoms.end(); k++) {
			k->first->setUnboundFrom(this, false);
		}
	}
	boundAtoms.clear();
	oneThreeAtoms.clear();
	oneFourAtoms.clear();

}

void Atom::setUnboundFrom(Atom * _pAtom, bool _propagate) {

	map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator found = boundAtoms.find(_pAtom);
	if (found != boundAtoms.end()) {
		// atom is bound: erase it

		// find all atoms that are connected through _pAtom
		vector<Atom *> thirdAtoms;
		for (map<Atom*, map<Atom*, bool> >::iterator k = found->second.begin(); k!=found->second.end(); k++) {
			thirdAtoms.push_back(k->first);
		}


		// remove all the 1-3 and 1-4 references through _pAtom
		for (vector<Atom*>::iterator third=thirdAtoms.begin(); third!=thirdAtoms.end(); third++) {
			purge14mid(_pAtom, *third);
			purge13(_pAtom, *third);
		}

		// erase the bond (it was most likely erased by the purge functions but just in case)
		boundAtoms.erase(found);


		// remove all the 1-3 and 1-4 through _pAtom from the atoms bound to this
		for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k = boundAtoms.begin(); k!=boundAtoms.end(); k++) {

			k->first->purge14mid(this, _pAtom);
			k->first->purge13(this, _pAtom);
			// and remove all the 1-4 through this and _pAtom from the atoms bound to those atoms
			for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator kk = k->first->boundAtoms.begin(); kk!=k->first->boundAtoms.end(); kk++) {
				kk->first->purge14end(this, _pAtom);
			}
		}
		// call the same function on the other atoms (the _propagate=false means to call this back)
		if (_propagate) {
			_pAtom->setUnboundFrom(this, false);
		}

	}

}

void Atom::purge13(Atom * _pAtom2, Atom * _pAtom3) {
	/*********************************************************
	 *  Remove all 1-3 references that are 
	 *
	 *    this -- _pAtom2 -- _pAtom3
	 *
	 *  corresponding to
	 *    boundAtoms[_pAtom2][_pAtom3]
	 *    oneThreeAtoms[_pAtom3][_pAtom2]
	 *********************************************************/
	map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator found1=boundAtoms.find(_pAtom2);
	if (found1 != boundAtoms.end()) {
		map<Atom*, map<Atom*, bool> >::iterator found2 = found1->second.find(_pAtom3);
		if (found2 != found1->second.end()) {
			found1->second.erase(found2);
		}
	}
	map<Atom*, map<Atom*, bool> >::iterator found3 = oneThreeAtoms.find(_pAtom3);
	if (found3 != oneThreeAtoms.end()) {
		map<Atom*, bool>::iterator found4 = found3->second.find(_pAtom2);
		if (found4 != found3->second.end()) {
			found3->second.erase(found4);
		}
		if (found3->second.size() == 0) {
			oneThreeAtoms.erase(found3);
		}
	}
}

void Atom::purge14mid(Atom * _pAtom2, Atom * _pAtom3) {
	/*********************************************************
	 *  Remove all 1-4 references that are 
	 *
	 *    this -- _pAtom2 -- _pAtom3 -- any atom
	 *
	 *  corresponding to
	 *    boundAtoms[_pAtom2][_pAtom3][any]
	 *    oneFourAtoms[any][_pAtom2][_pAtom3]
	 *********************************************************/
	vector<Atom*> fourthAtoms;
	map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator found1=boundAtoms.find(_pAtom2);
	if (found1 != boundAtoms.end()) {
		map<Atom*, map<Atom*, bool> >::iterator found2 = found1->second.find(_pAtom3);
		if (found2 != found1->second.end()) {
			for (map<Atom*, bool>::iterator k=found2->second.begin(); k!=found2->second.end(); k++) {
				fourthAtoms.push_back(k->first);
			}
			found1->second.erase(found2);
		}
	}

	for (vector<Atom*>::iterator k=fourthAtoms.begin(); k!=fourthAtoms.end(); k++) {
		map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator found3 = oneFourAtoms.find(*k);
		if (found3 != oneFourAtoms.end()) {
			map<Atom*, map<Atom*, bool> >::iterator found4 = found3->second.find(_pAtom2);
			if (found4 != found3->second.end()) {
				map<Atom*, bool>::iterator found5 = found4->second.find(_pAtom3);
				if (found5 != found4->second.end()) {
					found4->second.erase(found5);
				}
				if (found4->second.size() == 0) {
					found3->second.erase(found4);
				}
			}
			if (found3->second.size() == 0) {
				oneFourAtoms.erase(found3);
			}
		}
	}
}

void Atom::purge14end(Atom * _pAtom3, Atom * _pAtom4) {
	/*********************************************************
	 *  Remove all 1-4 references that are 
	 *
	 *    this -- any atom -- _pAtom3 -- _pAtom4
	 *
	 *  corresponding to
	 *    boundAtoms[any][_pAtom3][_pAtom4]
	 *    oneFourAtoms[_pAtom4][any][_pAtom3]
	 *********************************************************/
	for (map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		map<Atom*, map<Atom*, bool> >::iterator found1 = k->second.find(_pAtom3);
		if (found1 != k->second.end()) {
			map<Atom*, bool>::iterator found2 = found1->second.find(_pAtom4);
			if (found2 != found1->second.end()) {
				found1->second.erase(found2);
			}
		}
	}
	map<Atom*, map<Atom*, map<Atom*, bool> > >::iterator found3 = oneFourAtoms.find(_pAtom4);
	if (found3 != oneFourAtoms.end()) {
		for (map<Atom*, map<Atom*, bool> >::iterator k=found3->second.begin(); k!=found3->second.end();) {
			map<Atom*, bool>::iterator found4 = k->second.find(_pAtom3);
			if (found4 != k->second.end()) {
				k->second.erase(found4);
			}
			map<Atom*, map<Atom*, bool> >::iterator tmp = k;
			k++;
			if (tmp->second.size() == 0) {
				found3->second.erase(tmp);
			}
		}
		if (found3->second.size() == 0) {
			oneFourAtoms.erase(found3);
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
CartesianPoint& Atom::getGroupGeometricCenter(unsigned int _stamp) {
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

void Atom::copyAllCoor(const Atom _a) {
	//clearSavedCoor(); // let's get rid of any save coordinates

	//cout << "UUU atom has " << pCoorVec.size() << " coors, other atom has " << _a.pCoorVec.size() << " coors" << endl;

	// add elements if needed
	while (pCoorVec.size() < _a.pCoorVec.size()) {
		pCoorVec.push_back(new CartesianPoint(0.0, 0.0, 0.0));
		//cout << "UUU pushed new point" << endl;
	}
	//cout << "UUU added coordinates, now pCoorVec has " << pCoorVec.size() << " coordinates" << endl;
	// delete elements if needed
	if (pCoorVec.size() > _a.pCoorVec.size()) {
		for (unsigned int i=_a.pCoorVec.size(); i<pCoorVec.size(); i++) {
			delete pCoorVec[i];
			//cout << "UUU deleted " << i << "-th pointer" << endl;
		}
		pCoorVec.erase(pCoorVec.begin()+_a.pCoorVec.size(), pCoorVec.end());
	}
	//cout << "UUU erased coordinates, now pCoorVec has " << pCoorVec.size() << " coordinates" << endl;
	// assign coordinates
	for (unsigned int i=0; i<_a.pCoorVec.size(); i++) {
		*(pCoorVec[i]) = *(_a.pCoorVec[i]);
	}
	// set the current coordinate (currentIterator)
	currentCoorIterator = pCoorVec.begin() + _a.getActiveConformation();
	//cout << "UUU the current coordinates are " << getActiveConformation() << endl;
}

bool Atom::hideAltCoorAbsIndex(unsigned int _absoluteIndex) {
	// hide a coor based on relative index
	if (pCoorVec.size() == 1) {
		// cannot hide all coors
		return false;
	} else if (_absoluteIndex > pCoorVec.size() + pHiddenCoorVec.size()) {
		// the index is too big
		return false;
	}
	// calculate the relative index of the hidden position
	unsigned int relIndex;
	unsigned int indexInHidden;
	if (hiddenCoorIndeces.size() > 0 && _absoluteIndex < hiddenCoorIndeces[0]) {
		// it is at the beginning of the hidden positions
		relIndex = _absoluteIndex;
		indexInHidden = 0;
	} else {
		// default it to after all hidden
		relIndex = _absoluteIndex - hiddenCoorIndeces.size();
		indexInHidden = hiddenCoorIndeces.size();
		// check if it is in the middle
		for (unsigned int i=1; i<hiddenCoorIndeces.size(); i++) {
			if (_absoluteIndex < hiddenCoorIndeces[i]) {
				relIndex = _absoluteIndex - i;
				indexInHidden = i;
			}
		}
	}

	return hideAltCoors(_absoluteIndex, relIndex, indexInHidden);
}

bool Atom::hideAltCoorRelIndex(unsigned int _relativeIndex) {
	// hide a coor based on relative index
	if (pCoorVec.size() == 1) {
		// cannot hide all coors
		return false;
	} else if (_relativeIndex > pCoorVec.size()) {
		// the index is too big
		return false;
	}
	// calculate the absolute index of the hidden position
	unsigned int absIndex;
	unsigned int indexInHidden;
	if (hiddenCoorIndeces.size() > 0 && _relativeIndex < hiddenCoorIndeces[0]) {
		// it is at the beginning of the hidden positions
		absIndex = _relativeIndex;
		indexInHidden = 0;
	} else {
		// default it to after all hidden
		absIndex = _relativeIndex;
		indexInHidden = hiddenCoorIndeces.size();
		// check if it is in the middle
		for (unsigned int i=0; i<hiddenCoorIndeces.size(); i++) {
			if (hiddenCoorIndeces[i] <= absIndex) {
				absIndex++;
				indexInHidden = i+1;
			}
		}
	}

	return hideAltCoors(absIndex, _relativeIndex, indexInHidden);
}

void Atom::convertRelToAbs(unsigned int _relativeIndex, unsigned int & _absoluteIndex, unsigned int & _indexInHidden) const {
	// calculate the absolute index of the hidden position
	if (hiddenCoorIndeces.size() == 0) {
		// this conversion does not make sense if there are no hidden
		_absoluteIndex = 0;
		_indexInHidden = 0;
		return;
	}
	if (_relativeIndex < hiddenCoorIndeces[0]) {
		// it is at the beginning of the hidden positions
		_absoluteIndex = _relativeIndex;
		_indexInHidden = 0;
	} else {
		// default it to after all hidden
		_absoluteIndex = _relativeIndex + hiddenCoorIndeces.size();
		_indexInHidden = hiddenCoorIndeces.size();
		// check if it is in the middle
		for (unsigned int i=1; i<hiddenCoorIndeces.size(); i++) {
			if (_relativeIndex < hiddenCoorIndeces[i]) {
				_absoluteIndex = _relativeIndex + i;
				_indexInHidden = i;
			}
		}
	}
}

void Atom::convertAbsToRel(unsigned int _absoluteIndex, unsigned int & _relativeIndex, unsigned int & _indexInHidden) const {
	// calculate the relative index of the hidden position
	if (hiddenCoorIndeces.size() == 0) {
		// this conversion does not make sense if there are no hidden
		_relativeIndex = 0;
		_indexInHidden = 0;
		return;
	}
	if (hiddenCoorIndeces.size() > 0 && _absoluteIndex < hiddenCoorIndeces[0]) {
		// it is at the beginning of the hidden positions
		_relativeIndex = _absoluteIndex;
		_indexInHidden = 0;
	} else {
		// default it to after all hidden
		_relativeIndex = _absoluteIndex - hiddenCoorIndeces.size();
		_indexInHidden = hiddenCoorIndeces.size();
		// check if it is in the middle
		for (unsigned int i=1; i<hiddenCoorIndeces.size(); i++) {
			if (_absoluteIndex < hiddenCoorIndeces[i]) {
				_relativeIndex = _absoluteIndex - i;
				_indexInHidden = i;
			}
		}
	}
}


bool Atom::hideAllAltCoorsButOneRelIndex(unsigned int _keepThisIndex) {
	// turns all alt coor of except one, expressed as relative index
	unsigned int absIndex;
	unsigned int indexInHidden;
	convertRelToAbs(_keepThisIndex, absIndex, indexInHidden);
	return hideAllAltCoorsButOneAbsIndex(absIndex);
}
bool Atom::hideAllAltCoorsButOneAbsIndex(unsigned int _keepThisIndex) {
	if (_keepThisIndex >= pCoorVec.size() + pHiddenCoorVec.size()) {
		return false;
	}
	// turns all alt coor of except one, expressed as absolute index
	vector<CartesianPoint*> tmpHidden;
	vector<CartesianPoint*> tmpCoor;
	vector<unsigned int> tmpIndeces;
	unsigned int hI = 0;
	bool found = false;
	for (unsigned int i=0; i<pCoorVec.size(); i++) {
		// add from the hidden
		while (hI < hiddenCoorIndeces.size() && hiddenCoorIndeces[hI] == tmpHidden.size()+tmpCoor.size()) {
			if (!found && tmpHidden.size() == _keepThisIndex) {
				// coordinate to be keep unhidden
				tmpCoor.push_back(pHiddenCoorVec[hI]);
				found = true;
			} else {
				tmpHidden.push_back(pHiddenCoorVec[hI]);
				tmpIndeces.push_back(tmpCoor.size()+tmpHidden.size()-1);
			}
			hI++;
		}
		if (!found && tmpHidden.size() == _keepThisIndex) {
			// coordinate to be keep unhidden
			tmpCoor.push_back(pCoorVec[i]);
			found = true;
		} else {
			tmpHidden.push_back(pCoorVec[i]);
			tmpIndeces.push_back(tmpCoor.size()+tmpHidden.size()-1);
		}
	}
	// add any remining hidden
	while (hI < hiddenCoorIndeces.size()) {
		if (!found && tmpHidden.size() == _keepThisIndex) {
			// coordinate to be keep unhidden
			tmpCoor.push_back(pHiddenCoorVec[hI]);
			found = true;
		} else {
			tmpHidden.push_back(pHiddenCoorVec[hI]);
			tmpIndeces.push_back(tmpCoor.size()+tmpHidden.size()-1);
		}
		hI++;
	}
	pCoorVec = tmpCoor;
	pHiddenCoorVec = tmpHidden;
	hiddenCoorIndeces = tmpIndeces;
	setActiveConformation(0);
	return true;
	
}

bool Atom::hideAltCoors(unsigned int _absoluteIndex, unsigned int _relativeIndex, unsigned int _indexInHiddenn) {
	// private function that takes care of both absolute and relative hides

	// recalculate the active conformation
	unsigned int active = getActiveConformation();
	if (active > _relativeIndex || active == pCoorVec.size()) {
		active--;
	}

	// move the pCoor to the hidden vector
	pHiddenCoorVec.insert(pHiddenCoorVec.begin()+_indexInHiddenn, pCoorVec[_relativeIndex]);
	// record the absolute index of the hidden coordinate
	hiddenCoorIndeces.insert(hiddenCoorIndeces.begin()+_indexInHiddenn, _absoluteIndex);
	// alter the pCoorVec
	pCoorVec.erase(pCoorVec.begin()+_relativeIndex);
	setActiveConformation(active);
	return true;

}

bool Atom::unhideAltCoorAbsIndex(unsigned int _absoluteIndex) {
	if (_absoluteIndex >= pCoorVec.size() + pHiddenCoorVec.size()) {
		return false;
	}
	if (pHiddenCoorVec.size() == 0) {
		// nothing to do
		return true;
	}
	// unhide a specific coor based on absolute index
	unsigned int currHidden = 0;
	unsigned int curr = getActiveConformation();
	for (unsigned int i=0; i<pCoorVec.size(); i++) {
		while (currHidden < hiddenCoorIndeces.size() && hiddenCoorIndeces[currHidden] == i+currHidden) {
			if (i+currHidden == _absoluteIndex) {
				pCoorVec.insert(pCoorVec.begin()+i, pHiddenCoorVec[currHidden]);
				pHiddenCoorVec.erase(pHiddenCoorVec.begin()+currHidden);
				hiddenCoorIndeces.erase(hiddenCoorIndeces.begin()+currHidden);
				if (i <= curr) {
					curr++;
				}
				setActiveConformation(curr);
				return true;
			}
			currHidden++;
		}
		if (i+currHidden == _absoluteIndex) {
			// found, was already not hidden
			return true;
		}
	}
	for (unsigned int i=currHidden; i<hiddenCoorIndeces.size(); i++) {
		if (i+pCoorVec.size() == _absoluteIndex) {
			pCoorVec.insert(pCoorVec.end(), pHiddenCoorVec[i]);
			pHiddenCoorVec.erase(pHiddenCoorVec.begin()+i);
			hiddenCoorIndeces.erase(hiddenCoorIndeces.begin()+i);
			setActiveConformation(curr);
			cout << "UUU1 found " << i << " + " << pCoorVec.size() << " = " << _absoluteIndex << endl;
			return true;
		}
	}
	return false;
}

bool Atom::hideAllAltCoorsButFirstN(unsigned int _numberToKeepAbsIndex) {
	if (_numberToKeepAbsIndex == 0) {
		return false;
	}
	// turns all alt coor of except the first N, expressed as absolute index
	vector<CartesianPoint*> tmpHidden;
	vector<CartesianPoint*> tmpCoor;
	vector<unsigned int> tmpIndeces;

	unsigned int currHidden = 0;
	unsigned int curr = getActiveConformation();
	for (unsigned int i=0; i<pCoorVec.size(); i++) {
		while (currHidden < hiddenCoorIndeces.size() && hiddenCoorIndeces[currHidden] == i+currHidden) {
			// insert any hidden alt confs
			if (i+currHidden < _numberToKeepAbsIndex) {
				// visible
				tmpCoor.push_back(pHiddenCoorVec[currHidden]);
			} else {
				//hidden
				tmpHidden.push_back(pHiddenCoorVec[currHidden]);
				tmpIndeces.push_back(i+currHidden);
			}
			currHidden++;
		}
		// insert the non-hidden alt conf
		if (i+currHidden < _numberToKeepAbsIndex) {
			tmpCoor.push_back(pCoorVec[i]);
		} else {
			tmpHidden.push_back(pCoorVec[i]);
			tmpIndeces.push_back(i+currHidden);
		}
		// recalculate the current alt conf
		if (i == curr) {
			if (i+currHidden >= _numberToKeepAbsIndex) {
				// the original is now hidden, set it to 0
				curr = 0;
			} else {
				// stay with the original
				curr = i+currHidden;
			}
		}
	}
	// insert any terminal hidden alt conf
	for (unsigned int i=currHidden; i<hiddenCoorIndeces.size(); i++) {
		if (i+pCoorVec.size() < _numberToKeepAbsIndex) {
			tmpCoor.push_back(pHiddenCoorVec[i]);
		} else {
			tmpHidden.push_back(pHiddenCoorVec[i]);
			tmpIndeces.push_back(i+pCoorVec.size());
		}
	}
	pCoorVec = tmpCoor;
	pHiddenCoorVec = tmpHidden;
	hiddenCoorIndeces = tmpIndeces;
	setActiveConformation(curr);
	if (_numberToKeepAbsIndex > pCoorVec.size()) {
		return false;
	} else {
		return true;	
	}
}

bool Atom::unhideAllAltCoors() {
	// turns all alt coor of except the first N, expressed as absolute index
	vector<CartesianPoint*> tmpCoor;

	unsigned int currHidden = 0;
	unsigned int curr = getActiveConformation();
	for (unsigned int i=0; i<pCoorVec.size(); i++) {
		while (currHidden < hiddenCoorIndeces.size() && hiddenCoorIndeces[currHidden] == i+currHidden) {
			// insert any hidden alt confs
			tmpCoor.push_back(pHiddenCoorVec[currHidden]);
			currHidden++;
		}
		// insert the non-hidden alt conf
		tmpCoor.push_back(pCoorVec[i]);
		// recalculate the current alt conf
		if (i == curr) {
			// stay with the original
			curr = i+currHidden;
		}
	}
	// insert any terminal hidden alt conf
	for (unsigned int i=currHidden; i<hiddenCoorIndeces.size(); i++) {
		tmpCoor.push_back(pHiddenCoorVec[i]);
	}
	pCoorVec = tmpCoor;
	pHiddenCoorVec.clear();
	hiddenCoorIndeces.clear();
	setActiveConformation(curr);
	return true;	
}

std::string Atom::toString() const {
	std::string qm = " ";
	if (!hasCoordinates) {qm = "?";};
	std::string act = "+";
	if (!getActive()) {
		act="-";
	}
	char c [100];
	std::string out;
	if (toStringFormat == 0) {
		// print just the current conformation
		// CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		sprintf(c, "%-4s %-3s %4d%1s %1s [%10.3f %10.3f %10.3f]%1s(conf %3u/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), (*currentCoorIterator)->getX(), (*currentCoorIterator)->getY(), (*currentCoorIterator)->getZ(), qm.c_str(), getActiveConformation()+1, getNumberOfAltConformations(), act.c_str());
		out = (std::string)c;
	} else if (toStringFormat == 1) {
		// print all conformations that are not hidden
		// CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		// CA   LEU   37   *[     3.000      0.000      0.000] (conf   2/  3) +   active conf (*)
		// CA   LEU   37    [     4.000      0.000      0.000] (conf   3/  3) +
		unsigned int curr = getActiveConformation();
		for (unsigned int i=0; i<pCoorVec.size(); i++) {
			string active = " ";
			if (i==curr) {
				active = "*";
			}
			sprintf(c, "%-4s %-3s %4d%1s %1s%1s[%10.3f %10.3f %10.3f]%1s(conf %3u/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), active.c_str(), pCoorVec[i]->getX(), pCoorVec[i]->getY(), pCoorVec[i]->getZ(), qm.c_str(), i+1, getNumberOfAltConformations(), act.c_str());
			out += (std::string)c + (std::string)"\n";
		}
	} else if (toStringFormat == 2) {
		// print all conformations including the hidden ones
		// CA   LEU   37    [     0.000      0.000      0.000] (conf   1/  3) +
		// CA   LEU   37    [     2.000      0.000      0.000] (conf   H/  1) +   hidden conf
		// CA   LEU   37   *[     3.000      0.000      0.000] (conf   2/  3) +
		// CA   LEU   37    [     4.000      0.000      0.000] (conf   3/  3) +
		unsigned int currHidden = 0;
		unsigned int curr = getActiveConformation();
		unsigned int tot = getNumberOfAltConformations();
		for (unsigned int i=0; i<pCoorVec.size(); i++) {
			string active = " ";
			if (i==curr) {
				active = "*";
			}
			while (currHidden < hiddenCoorIndeces.size() && hiddenCoorIndeces[currHidden] == i+currHidden) {
				sprintf(c, "%-4s %-3s %4d%1s %1s [%10.3f %10.3f %10.3f]%1s(conf   H/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), pHiddenCoorVec[currHidden]->getX(), pHiddenCoorVec[currHidden]->getY(), pHiddenCoorVec[currHidden]->getZ(), qm.c_str(), (unsigned int)pHiddenCoorVec.size(), act.c_str());
				out += (std::string)c + (std::string)"\n";
				currHidden++;
			}
			sprintf(c, "%-4s %-3s %4d%1s %1s%1s[%10.3f %10.3f %10.3f]%1s(conf %3u/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), active.c_str(), pCoorVec[i]->getX(), pCoorVec[i]->getY(), pCoorVec[i]->getZ(), qm.c_str(), i+1, tot, act.c_str());
			out += (std::string)c + (std::string)"\n";
		}
		for (unsigned int i=currHidden; i<hiddenCoorIndeces.size(); i++) {
			sprintf(c, "%-4s %-3s %4d%1s %1s [%10.3f %10.3f %10.3f]%1s(conf   H/%3u) %1s", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), pHiddenCoorVec[i]->getX(), pHiddenCoorVec[i]->getY(), pHiddenCoorVec[i]->getZ(), qm.c_str(), (unsigned int)pHiddenCoorVec.size(), act.c_str());
			out += (std::string)c + (std::string)"\n";
		}
		/*
		//SIMPLER PRINTING FOR TESTING
		out += (string)"\n";
		for (unsigned int i=0; i<pCoorVec.size(); i++) {
			string active = " ";
			if (i==curr) {
				active = "*";
			}
			sprintf(c, "%-4s %-3s %4d%1s %1s%1s[%10.3f %10.3f %10.3f]%1s(conf %3u/%3u) %1s %x", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), active.c_str(), pCoorVec[i]->getX(), pCoorVec[i]->getY(), pCoorVec[i]->getZ(), qm.c_str(), i+1, tot, act.c_str(), pCoorVec[i]);
			out += (std::string)c + (std::string)"\n";
		}
		for (unsigned int i=0; i<hiddenCoorIndeces.size(); i++) {
			sprintf(c, "%-4s %-3s %4d%1s %1s [%10.3f %10.3f %10.3f]%1s(conf   H/%3u) %1s %x %u", name.c_str(), getResidueName().c_str(), getResidueNumber(), getResidueIcode().c_str(), getChainId().c_str(), pHiddenCoorVec[i]->getX(), pHiddenCoorVec[i]->getY(), pHiddenCoorVec[i]->getZ(), qm.c_str(), (unsigned int)pHiddenCoorVec.size(), act.c_str(), pHiddenCoorVec[i], hiddenCoorIndeces[i]);
			out += (std::string)c + (std::string)"\n";
		}
		*/

	}

	return out;
}

void Atom::removeFromIc() {
	for (vector<IcEntry*>::iterator k=icEntries.begin(); k<icEntries.end(); k++) {
		if (*k != NULL) {
			IcTable * pParentTable = (*k)->getParentTable();
			if (pParentTable != NULL) {
				pParentTable->removeAtom(this);
				return;
			}
		}
	}
}


bool Atom::isPositionNterminal() const {
	if (pParentGroup != NULL) {
		return pParentGroup->isPositionNterminal();
	}
	return false;
}
bool Atom::isPositionCterminal() const {
	if (pParentGroup != NULL) {
		return pParentGroup->isPositionCterminal();
	}
	return false;
}

