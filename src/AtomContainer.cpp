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


#include "AtomContainer.h"

using namespace MSL;
using namespace std;


AtomContainer::AtomContainer() {
	addAtomsAsAltCoors_flag = false;
	setup();
}

AtomContainer::AtomContainer(const AtomPointerVector & _atoms) {
	addAtomsAsAltCoors_flag = false;
	setup();
	addAtoms(_atoms);
}

AtomContainer::AtomContainer(const AtomContainer & _AC) {
	addAtomsAsAltCoors_flag = false;
	pdbReader = NULL;
	pdbWriter = NULL;
	copy(_AC);
}

AtomContainer::~AtomContainer() {
	deletePointers();
}

void AtomContainer::operator=(const AtomContainer & _AC) {
	deletePointers();
	copy(_AC);
}

void AtomContainer::setup() {
	pdbReader = new PDBReader;
	pdbWriter = new PDBWriter;
	found = atomMap.end();
}

void AtomContainer::setup(std::stringstream& _str) {
        pdbReader = new PDBReader(_str);
        pdbWriter = new PDBWriter;
        found = atomMap.end();
}

void AtomContainer::copy(const AtomContainer & _AC) {
        deletePointers();
        setup();
        addAtoms(_AC.atoms);
}

void AtomContainer::reset() {
	deletePointers();
	setup();
}

void AtomContainer::deletePointers() {
	atomMap.clear();
	atomMapWithIdentities.clear();
	found = atomMap.end();
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
	if (pdbReader != NULL) delete pdbReader;
	if (pdbWriter != NULL) delete pdbWriter;
        pdbReader = NULL; pdbWriter = NULL;
}

void AtomContainer::addAtom(const Atom & _atom) {
	/*
	atoms.push_back(new Atom(_atom));
	atomMap[_atom.getAtomId()] = atoms.back();
	atomMapWithIdentities[_atom.getAtomOfIdentityId()] = atoms.back();
	found = atomMap.end();
	*/

	bool addedCoor = false;
	if (addAtomsAsAltCoors_flag) {
		//cout << "UUU set as true" << endl;
		if (_atom.getResidueName() == "") {
			// the new atom does not have a residue name
			string atomId = _atom.getAtomId();
			found = atomMap.find(atomId);
			if (found != atomMap.end()) {
				// add the coordinate to the atom
				CartesianPoint p(_atom.getX(), _atom.getY(), _atom.getZ());
				(found->second)->addAltConformation(p);
				addedCoor = true;
		//		cout << "UUU0 added alt coor" << endl;
			}
		} else {
			// has a residue name
			string atomId = _atom.getAtomOfIdentityId();
			found = atomMapWithIdentities.find(atomId);
			if (found != atomMapWithIdentities.end()) {
				// add the coordinate to the atom
				CartesianPoint p(_atom.getX(), _atom.getY(), _atom.getZ());
				(found->second)->addAltConformation(p);
				addedCoor = true;
		//		cout << "UUU1 added alt coor" << endl;
			}
		}
	}
	if (!addedCoor) {
		// an atom with the same id wasn't found: add a new atom
		atoms.push_back(new Atom(_atom));
		atomMap[_atom.getAtomId()] = atoms.back();
		atomMapWithIdentities[_atom.getAtomOfIdentityId()] = atoms.back();
		found = atomMap.end();
		//cout << "UUU0 added atom" << endl;
	}
}

void AtomContainer::addAtom(string _atomId, const CartesianPoint & _coor, string _element) {
	addAtom(Atom(_atomId, _coor, _element));
}

void AtomContainer::addAtom(string _atomId, double _x, double _y, double _z, string _element) {
	addAtom(Atom(_atomId, CartesianPoint(_x, _y, _z), _element));
}

void AtomContainer::addAtoms(const AtomPointerVector & _atoms) {
	for (AtomPointerVector::const_iterator k = _atoms.begin(); k != _atoms.end(); k++) {
		addAtom(*(*k));
	}
}

void AtomContainer::insertAtom(const Atom & _atom, unsigned int _skipPositions) {
	Atom * newAtom = new Atom(_atom);
	atoms.insert(atoms.begin()+_skipPositions, new Atom(_atom));

	atomMap[newAtom->getAtomId()] = newAtom;

	atomMapWithIdentities[newAtom->getAtomOfIdentityId()] = newAtom;
	found = atomMap.end();
}

void AtomContainer::insertAtom(string _atomId, const CartesianPoint & _coor, unsigned int _skipPositions) {
	insertAtom(Atom(_atomId, _coor), _skipPositions);
}

void AtomContainer::insertAtoms(const AtomPointerVector & _atoms, unsigned int _skipPositions) {
	AtomPointerVector newAtoms;
	for (AtomPointerVector::const_iterator k = _atoms.begin(); k != _atoms.end(); k++) {
		newAtoms.push_back(new Atom(**k));
		atomMap[(*k)->getAtomId()] = newAtoms.back();
		atomMapWithIdentities[(*k)->getAtomOfIdentityId()] = newAtoms.back();
	}
	found = atomMap.end();
	atoms.insert(atoms.begin()+_skipPositions, newAtoms.begin(), newAtoms.end());
}

void AtomContainer::removeAtom(unsigned int _n) {
	for (map<string, Atom*>::iterator k=atomMap.begin(); k!=atomMap.end(); k++) {
		if (k->second == atoms[_n]) {
			atomMap.erase(k);
			break;
		}
	}
	atoms.erase(atoms.begin()+_n);
}

bool AtomContainer::removeAtom(string _atomId) {
	map<string, Atom*>::iterator foundAtom=atomMap.find(_atomId);

	if (foundAtom!=atomMap.end()) {
		// erase from the atoms list
		for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
			if ((*foundAtom).second == *k) {
				// deallocate from memory and remove from list
				delete *k;
				atoms.erase(k);
				// erase from the map
				return true;
			}
		}
		atomMap.erase(foundAtom);
		found=atomMap.end();
		return true;
	}
	return false;
}


bool AtomContainer::atomExists(string _atomId) {
	string chainid;
	int resnum;
	string icode;
	string atomName;
	bool OK = MslTools::parseAtomId(_atomId, chainid, resnum, icode, atomName);
	if (OK) {
		// it was in "A 37 CA" format (without identity)
		string key = MslTools::getAtomId(chainid, resnum, icode, atomName);
		found = atomMap.find(key);
		return found != atomMap.end();
	}
	string identity;
	OK = MslTools::parseAtomOfIdentityId(_atomId, chainid, resnum, icode, identity, atomName);
	if (OK) {
		// it was in "A 37 ILE CA" format (with the identity)
		string key = MslTools::getAtomOfIdentityId(chainid, resnum, icode, identity, atomName);
		found = atomMapWithIdentities.find(key);
		return found != atomMapWithIdentities.end();
	}
	found = atomMap.end();
	return false;

}

void AtomContainer::updateAtomMap(Atom * _atom) {
	// if an atom changes its name it calls this function
	// to update the map (through its group)
	for (map<string, Atom*>::iterator k=atomMap.begin(); k!=atomMap.end(); k++) {
		if (k->second == _atom) {
			atomMap.erase(k);
			atomMap[_atom->getAtomId()] = _atom;
			break;
		}
	}
	for (map<string, Atom*>::iterator k=atomMapWithIdentities.begin(); k!=atomMapWithIdentities.end(); k++) {
		if (k->second == _atom) {
			atomMapWithIdentities.erase(k);
			atomMapWithIdentities[_atom->getAtomOfIdentityId()] = _atom;
			break;
		}
	}
	found = atomMap.end();
}

Atom & AtomContainer::getAtom(string _atomId) {
	if (atomExists(_atomId)) {
		return *(found->second);
	}
	cerr << "ERROR 38182: atom " << _atomId << " not found in Atom & AtomContainer::getAtom(string _chain_resnum_name)" << endl;
	exit(38182);
}

void AtomContainer::removeAllAtoms(){
  reset();
}
