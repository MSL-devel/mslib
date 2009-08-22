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

#include "IcTable.h"

IcTable::IcTable() {
}

IcTable::IcTable(const IcTable & _ic) {
	//copy(_ic);
}

IcTable::~IcTable() {
}

void IcTable::mapValues(IcEntry * _ic) {
	vector<double> & vals = _ic->getValues();
	Atom * pA1 = _ic->getAtom1();
	Atom * pA2 = _ic->getAtom2();
	Atom * pA3 = _ic->getAtom3();
	Atom * pA4 = _ic->getAtom4();
	if (_ic->isImproper()) {
		bondMap[pA1][pA3].push_back(&(vals[0]));
		bondMap[pA3][pA1].push_back(&(vals[0]));
		angleMap[pA1][pA3][pA2].push_back(&(vals[1]));
		angleMap[pA2][pA3][pA1].push_back(&(vals[1]));
	} else {
		bondMap[pA1][pA2].push_back(&(vals[0]));
		bondMap[pA2][pA1].push_back(&(vals[0]));
		angleMap[pA1][pA2][pA3].push_back(&(vals[1]));
		angleMap[pA3][pA2][pA1].push_back(&(vals[1]));
	}
	dihedralMap[pA1][pA2][pA3][pA4].push_back(&(vals[2]));
	dihedralMap[pA4][pA3][pA2][pA1].push_back(&(vals[2]));
	angleMap[pA2][pA3][pA4].push_back(&(vals[3]));
	angleMap[pA4][pA3][pA2].push_back(&(vals[3]));
	bondMap[pA3][pA4].push_back(&(vals[4]));
	bondMap[pA4][pA3].push_back(&(vals[4]));
}

void IcTable::updateMap(Atom * _pOldAtom, Atom * _pNewAtom) {
	replaceInMap(bondMap, _pOldAtom, _pNewAtom);
	replaceInMap(angleMap, _pOldAtom, _pNewAtom);
	replaceInMap(dihedralMap, _pOldAtom, _pNewAtom);
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom) {
	map<Atom*, map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		// erase and create a new entry with the new atom
		tmpMap = k->second;
		_map.erase(k);
		_map[_pNewAtom] = tmpMap;
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	map<Atom*, map<Atom*, map<Atom*, vector<double*> > > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, map<Atom*, vector<double*> > > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		_map[_pNewAtom] = tmpMap;
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, vector<double*> > > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	map<Atom*, map<Atom*, vector<double*> > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, vector<double*> > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		_map[_pNewAtom] = tmpMap;
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, vector<double*> > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	map<Atom*, vector<double*> >::iterator k = _map.find(_pOldAtom);
	vector<double*> tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		_map[_pNewAtom] = tmpMap;
	}
}

bool IcTable::editBond(Atom * _pAtom1, Atom * _pAtom2, double _newValue) {
	if (bondMap.find(_pAtom1) != bondMap.end() && bondMap[_pAtom1].find(_pAtom2) != bondMap[_pAtom1].end()) {
		for (vector<double*>::iterator k = bondMap[_pAtom1][_pAtom2].begin(); k!=bondMap[_pAtom1][_pAtom2].end(); k++) {
			*(*k) = _newValue;
			//cout << "UUU XXX Edited bond for " << *_pAtom1 << ":" << *_pAtom2 << " with value " << _newValue << endl;
		}
		return true;
	}
	return false;
}

bool IcTable::editAngle(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, double _newValue) {
	if (angleMap.find(_pAtom1) != angleMap.end() && angleMap[_pAtom1].find(_pAtom2) != angleMap[_pAtom1].end() && angleMap[_pAtom1][_pAtom2].find(_pAtom3) != angleMap[_pAtom1][_pAtom2].end()) {
		for (vector<double*>::iterator k = angleMap[_pAtom1][_pAtom2][_pAtom3].begin(); k!=angleMap[_pAtom1][_pAtom2][_pAtom3].end(); k++) {
			*(*k) = _newValue / 180.0 * M_PI;
			//cout << "UUU XXX Edited dihedral for " << *_pAtom1 << ":" << *_pAtom2 << ":" << *_pAtom3 << " with value " << _newValue << endl;
		}
		return true;
	}
	return false;
}

bool IcTable::editDihedral(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _newValue) {
	if (dihedralMap.find(_pAtom1) != dihedralMap.end() && dihedralMap[_pAtom1].find(_pAtom2) != dihedralMap[_pAtom1].end() && dihedralMap[_pAtom1][_pAtom2].find(_pAtom3) != dihedralMap[_pAtom1][_pAtom2].end() && dihedralMap[_pAtom1][_pAtom2][_pAtom3].find(_pAtom4) != dihedralMap[_pAtom1][_pAtom2][_pAtom3].end()) {
		for (vector<double*>::iterator k = dihedralMap[_pAtom1][_pAtom2][_pAtom3][_pAtom4].begin(); k!=dihedralMap[_pAtom1][_pAtom2][_pAtom3][_pAtom4].end(); k++) {
			*(*k) = _newValue / 180.0 * M_PI;
			//cout << "UUU XXX Edited dihedral for " << *_pAtom1 << ":" << *_pAtom2 << ":" << *_pAtom3 << ":" << *_pAtom4 << " with value " << _newValue << endl;
		}
		return true;
	} else if (dihedralMap.find(_pAtom1) != dihedralMap.end() && dihedralMap[_pAtom1].find(_pAtom3) != dihedralMap[_pAtom1].end() && dihedralMap[_pAtom1][_pAtom3].find(_pAtom2) != dihedralMap[_pAtom1][_pAtom3].end() && dihedralMap[_pAtom1][_pAtom3][_pAtom2].find(_pAtom4) != dihedralMap[_pAtom1][_pAtom3][_pAtom2].end()) {
		// inversed order of atoms 2-3, edit -value
		for (vector<double*>::iterator k = dihedralMap[_pAtom1][_pAtom3][_pAtom2][_pAtom4].begin(); k!=dihedralMap[_pAtom1][_pAtom3][_pAtom2][_pAtom4].end(); k++) {
			*(*k) = -_newValue / 180.0 * M_PI;
			//cout << "UUU XXX Edited dihedral for " << *_pAtom1 << ":" << *_pAtom2 << ":" << *_pAtom3 << ":" << *_pAtom4 << " (inverted)" << endl;
		}
		return true;
	}
	return false;
}

