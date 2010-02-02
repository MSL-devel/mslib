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

#include "AtomContainer.h"

AtomContainer::AtomContainer() {
	setup();
}

AtomContainer::AtomContainer(const AtomVector & _atoms) {
	setup();
	addAtoms(_atoms);
}

AtomContainer::AtomContainer(const AtomContainer & _AC) {
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

void AtomContainer::copy(const AtomContainer & _AC) {
	deletePointers();
}

void AtomContainer::reset() {
	deletePointers();
	setup();
}

void AtomContainer::deletePointers() {

	atomMap.clear();
	found = atomMap.end();
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
	delete pdbReader;
	delete pdbWriter;
}

void AtomContainer::addAtom(const Atom & _atom) {
	atoms.push_back(new Atom(_atom));
	string key = getMapKey(_atom.getChainId(), _atom.getResidueNumber(), _atom.getResidueIcode(), _atom.getName());
	atomMap[key] = atoms.back();
	found = atomMap.end();
}

void AtomContainer::addAtom(string _name, const CartesianPoint & _coor) {
	addAtom(Atom(_name, _coor));
}

void AtomContainer::addAtom(string _name, double _x, double _y, double _z) {
	addAtom(Atom(_name, CartesianPoint(_x, _y, _z)));
}

void AtomContainer::addAtoms(const AtomVector & _atoms) {
	for (AtomVector::const_iterator k = _atoms.begin(); k != _atoms.end(); k++) {
		addAtom(*(*k));
	}
}

void AtomContainer::insertAtom(const Atom & _atom, unsigned int _skipPositions) {
	AtomVector newAtoms;
	newAtoms.push_back(new Atom(_atom));
	found = atomMap.end();
	string key = getMapKey(_atom.getChainId(), _atom.getResidueNumber(), _atom.getResidueIcode(), _atom.getName());
	atoms.insert(atoms.begin()+_skipPositions, newAtoms.begin(), newAtoms.end());
}

void AtomContainer::insertAtom(string _name, const CartesianPoint & _coor, unsigned int _skipPositions) {
	insertAtom(Atom(_name, _coor), _skipPositions);
}

void AtomContainer::insertAtoms(const AtomVector & _atoms, unsigned int _skipPositions) {
	AtomVector newAtoms;
	for (AtomVector::const_iterator k = _atoms.begin(); k != _atoms.end(); k++) {
		newAtoms.push_back(new Atom(**k));
		string key = getMapKey((*k)->getChainId(), (*k)->getResidueNumber(), (*k)->getResidueIcode(), (*k)->getName());
		atomMap[key] = newAtoms.back();
	}
	found = atomMap.end();
	atoms.insert(atoms.begin()+_skipPositions, newAtoms.begin(), newAtoms.end());
}

Atom & AtomContainer::operator()(string _chain_resnum_name) {
	exists(_chain_resnum_name);
	return getLastFoundAtom();
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

bool AtomContainer::removeAtom(string _chain_resnum_name) {
	map<string, Atom*>::iterator foundAtom=atomMap.find(_chain_resnum_name);

	if (foundAtom!=atomMap.end()) {
		// erase from the atoms list
		for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
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


bool AtomContainer::exists(string _chain_resnum_name) {
	vector<string> tokens = MslTools::tokenize(_chain_resnum_name, ",");
	if (tokens.size() != 3) {
		tokens = MslTools::tokenize( _chain_resnum_name, " ");
	}
	if (tokens.size() != 3) {
		found = atomMap.end();
		return false;
	}
	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}
	int resNum = 0;
	string iCode;
	MslTools::splitIntAndString(tokens[1], resNum, iCode);
	string key = getMapKey(tokens[0], resNum, iCode, tokens[2]);
	found = atomMap.find(key);
	return found != atomMap.end();

}

void AtomContainer::updateAtomMap(Atom * _atom) {
	// if an atom changes its name it calls this function
	// to update the map (through its group)
	for (map<string, Atom*>::iterator k=atomMap.begin(); k!=atomMap.end(); k++) {
		if (k->second == _atom) {
			atomMap.erase(k);
			string key = getMapKey(_atom->getChainId(), _atom->getResidueNumber(), _atom->getResidueIcode(), _atom->getName());
			atomMap[key] = _atom;
		}
	}
}

string AtomContainer::getMapKey(string _chainId, int _resNum, string _iCode, string _name) {
	char c [1000];
	sprintf(c, "%1s %04u %1s %4s", _chainId.c_str(), _resNum, _iCode.c_str(), _name.c_str());
	return (string)c;

}

