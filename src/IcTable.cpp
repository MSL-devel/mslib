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

#include "IcTable.h"

using namespace MSL;
using namespace std;


IcTable::IcTable() {
}

IcTable::IcTable(const IcTable & _ic) {
	//copy(_ic);
}

IcTable::~IcTable() {
}

void IcTable::deletePointers() {
    for(vector<IcEntry *>::iterator it = this->begin(); it != this->end(); ++it) {
        delete *it;
    }
}

void IcTable::mapValues(IcEntry * _ic) {
	vector<double> & vals = _ic->getValues();
	Atom * pA1 = _ic->getAtom1();
	Atom * pA2 = _ic->getAtom2();
	Atom * pA3 = _ic->getAtom3();
	Atom * pA4 = _ic->getAtom4();
	if (_ic->isImproper()) {
		if (pA1 != NULL && pA3 != NULL) {
			bondMap[pA1][pA3].push_back(&(vals[0]));
			bondMap[pA3][pA1].push_back(&(vals[0]));
		}
		if (pA1 != NULL && pA2 != NULL && pA3 != NULL) {
			angleMap[pA1][pA3][pA2].push_back(&(vals[1]));
			angleMap[pA2][pA3][pA1].push_back(&(vals[1]));
		}
	} else {
		if (pA1 != NULL && pA2 != NULL) {
			bondMap[pA1][pA2].push_back(&(vals[0]));
			bondMap[pA2][pA1].push_back(&(vals[0]));
		}
		if (pA1 != NULL && pA2 != NULL && pA3 != NULL) {
			angleMap[pA1][pA2][pA3].push_back(&(vals[1]));
			angleMap[pA3][pA2][pA1].push_back(&(vals[1]));
		}
	}
	if (pA1 != NULL && pA2 != NULL && pA3 != NULL && pA4 != NULL) {
		dihedralMap[pA1][pA2][pA3][pA4].push_back(&(vals[2]));
		dihedralMap[pA4][pA3][pA2][pA1].push_back(&(vals[2]));
	}
	if (pA2 != NULL && pA3 != NULL && pA4 != NULL) {
		angleMap[pA2][pA3][pA4].push_back(&(vals[3]));
		angleMap[pA4][pA3][pA2].push_back(&(vals[3]));
	}
	if (pA3 != NULL && pA4 != NULL) {
		bondMap[pA3][pA4].push_back(&(vals[4]));
		bondMap[pA4][pA3].push_back(&(vals[4]));
	}
}

void IcTable::removeAtom(Atom * pAtom) {
	bool found = false;
	for (vector<IcEntry*>::iterator k=begin(); k!=end(); k++) {
		if ((*k)->removeAtom(pAtom)) {
			found = true;
		}
	}
	if (found) {
		updateMap(pAtom, NULL);
	}
}

void IcTable::updateMap(Atom * _pOldAtom, Atom * _pNewAtom) {
	replaceInMap(bondMap, _pOldAtom, _pNewAtom);
	replaceInMap(angleMap, _pOldAtom, _pNewAtom);
	replaceInMap(dihedralMap, _pOldAtom, _pNewAtom);
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom) {
	if (_pOldAtom == NULL) {
		return;
	}
	map<Atom*, map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		// erase and create a new entry with the new atom
		tmpMap = k->second;
		_map.erase(k);
		if (_pNewAtom != NULL) {
			_map[_pNewAtom] = tmpMap;
		}
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, map<Atom*, vector<double*> > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	if (_pOldAtom == NULL) {
		return;
	}

	map<Atom*, map<Atom*, map<Atom*, vector<double*> > > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, map<Atom*, vector<double*> > > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		if (_pNewAtom != NULL) {
			_map[_pNewAtom] = tmpMap;
		}
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, map<Atom*, vector<double*> > > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	if (_pOldAtom == NULL) {
		return;
	}

	map<Atom*, map<Atom*, vector<double*> > >::iterator k = _map.find(_pOldAtom);
	map<Atom*, vector<double*> > tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		if (_pNewAtom != NULL) {
			_map[_pNewAtom] = tmpMap;
		}
	}

	// call the subroutine that will take care of second level replacements
	for (k=_map.begin(); k!=_map.end(); k++) {
		replaceInMap(k->second, _pOldAtom, _pNewAtom);
	}
}

void IcTable::replaceInMap(map<Atom*, vector<double*> > & _map, Atom * _pOldAtom, Atom * _pNewAtom){
	if (_pOldAtom == NULL) {
		return;
	}

	map<Atom*, vector<double*> >::iterator k = _map.find(_pOldAtom);
	vector<double*> tmpMap;

	// find elements of the map that have the pointer to the old atom
	if (k != _map.end()) {
		tmpMap = k->second;
		_map.erase(k);
		if (_pNewAtom != NULL) {
			_map[_pNewAtom] = tmpMap;
		}
	}
}

bool IcTable::editBond(Atom * _pAtom1, Atom * _pAtom2, double _newValue) {
	if (_pAtom1 == NULL || _pAtom2 == NULL) {
		return false;
	}
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
	if (_pAtom1 == NULL || _pAtom2 == NULL || _pAtom3 == NULL) {
		return false;
	}
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
	if (_pAtom1 == NULL || _pAtom2 == NULL || _pAtom3 == NULL || _pAtom4 == NULL) {
		return false;
	}
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

bool IcTable::seed() {
	/*********************************************************
	 *  Auto seeding, finds the first IC that seems proper for
	 *  seeding each chain.  The current mechanism for finding
	 *  seeds for different chains works well as long as the IC
	 *  don't mix chains, which it shouldn't, if it happens it
	 *  would probably mnot work
	 *********************************************************/
	double d1For = 0.0;
	double d1Rev = 0.0;
	double d2 = 0.0;
	double aFor = 0.0;
	double aRev = 0.0;
	bool out = false;
	map<string, bool> chainSeeded; // to keep track of what chains we already seeded
	for (IcTable::iterator k=this->begin(); k!=this->end(); k++) {
		bool foundForward = false;
		bool foundReverse = false;
		Atom * IC1 = (*k)->getAtom1();
		Atom * IC2 = (*k)->getAtom2();
		Atom * IC3 = (*k)->getAtom3();
		Atom * IC4 = (*k)->getAtom4();
		if ((IC1 != NULL && chainSeeded.find(IC1->getChainId()) != chainSeeded.end()) || (IC2 != NULL && chainSeeded.find(IC2->getChainId()) != chainSeeded.end()) || (IC3 != NULL && chainSeeded.find(IC3->getChainId()) != chainSeeded.end()) || (IC4 != NULL && chainSeeded.find(IC4->getChainId()) != chainSeeded.end())) {
			// chain already seeded
			continue;
		}
		out = false;
		if (IC2 == NULL || IC3 == NULL) {
			// need the two middle atoms for this to work
			continue;
		}
		if (IC1 != NULL) {
			// first atom is good
			d1For = (*k)->getDistance1();
			aFor  = (*k)->getAngle1Radians();
			foundForward = true;
		}
		if (IC4 != NULL) {
			// last atom is good
			d1Rev = (*k)->getDistance2();
			aRev  = (*k)->getAngle2Radians();
			foundReverse = true;
		}
		if (foundForward || foundReverse) {
			bool foundD2 = false;
			// we need the distance betwee IC2 and IC3
			for (IcTable::iterator l=this->begin(); l!=this->end(); l++) {
				if (l == k) {
					continue;
				}
				if ((*l)->areD1Atoms(IC2, IC3)) {
					d2 = (*l)->getDistance1();
					foundD2 = true;
					break;
				} else if ((*l)->areD2Atoms(IC2, IC3)) {
					d2 = (*l)->getDistance2();
					foundD2 = true;
					break;
				}
			}
			if (foundD2) {
				if (foundForward) {
					if ((*k)->isImproper()) {
						CartesianGeometry::seedRadians(IC1->getCoor(), IC3->getCoor(), IC2->getCoor(), d1For, d2, aFor);
					} else {
						CartesianGeometry::seedRadians(IC1->getCoor(), IC2->getCoor(), IC3->getCoor(), d1For, d2, aFor);
					}
					IC1->setHasCoordinates();
					IC2->setHasCoordinates();
					IC3->setHasCoordinates();
					chainSeeded[IC1->getChainId()] = true;
					chainSeeded[IC2->getChainId()] = true;
					chainSeeded[IC3->getChainId()] = true;
					out = true;
				} else {
					CartesianGeometry::seedRadians(IC4->getCoor(), IC3->getCoor(), IC2->getCoor(), d1Rev, d2, aRev);
					IC4->setHasCoordinates();
					IC3->setHasCoordinates();
					IC2->setHasCoordinates();
					chainSeeded[IC4->getChainId()] = true;
					chainSeeded[IC3->getChainId()] = true;
					chainSeeded[IC2->getChainId()] = true;
					out = true;
				}
			}
		}
	}

	return out;
}

bool IcTable::seed(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3) {

	if (_pAtom1 == NULL || _pAtom2 == NULL || _pAtom3 == NULL) {
		return false;
	}

	double d1 = 0.0;
	double d2 = 0.0;
	double a = 0.0;
	bool foundD1 = false;
	bool foundD2 = false;
	bool foundA = false;
	for (IcTable::iterator k=this->begin(); k!=this->end(); k++) {
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
		CartesianGeometry::seedRadians(_pAtom1->getCoor(), _pAtom2->getCoor(), _pAtom3->getCoor(), d1, d2, a);
		_pAtom1->setHasCoordinates();
		_pAtom2->setHasCoordinates();
		_pAtom3->setHasCoordinates();
	//	cout << "UUU seeded " << _pAtom1->getCoor() << _pAtom2->getCoor() << _pAtom3->getCoor() << " " << d1 << " " << d2 << " " << a << endl;
		return true;
	} else {
		return false;
	}
}

