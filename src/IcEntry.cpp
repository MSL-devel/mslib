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

#include "IcEntry.h"
#include "IcTable.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("IcEntry");

IcEntry::IcEntry() {
	setup(NULL, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, false);
}

IcEntry::IcEntry(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {

	setup(&_atom1, &_atom2, &_atom3, &_atom4, _d1, _a1, _dihe, _a2, _d2, _improperFlag);
}

IcEntry::IcEntry(const IcEntry & _ic) {
	setup(NULL, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, false);
	copy(_ic);
}

IcEntry::~IcEntry() {
	/**********************************************
	 * make sure that the atoms will remove
	 * their pointers to this IcEntry before it
	 * is deleted from memory
	 **********************************************/
	if (pAtom1 != NULL) {
		pAtom1->removeIcEntry(this);
	}
	if (pAtom4 != NULL) {
		pAtom4->removeIcEntry(this);
	}
}

void IcEntry::setup(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag) {
	pParentTable = NULL;

	pAtom1 = _pAtom1;
	pAtom2 = _pAtom2;
	pAtom3 = _pAtom3;
	pAtom4 = _pAtom4;

	vals = vector<double>(5, 0.0);
	vals[0]   = _d1;
	vals[1]   = _a1 * M_PI / 180.0;  // store as radians
	vals[2] = _dihe * M_PI / 180.0;
	vals[3]   = _a2 * M_PI / 180.0;
	vals[4]   = _d2;

	improperFlag = _improperFlag;

	if (pAtom1 != NULL) {
		pAtom1->addIcEntry(this);
	}
	if (pAtom4 != NULL) {
		pAtom4->addIcEntry(this);
	}
}

void IcEntry::copy(const IcEntry & _ic) {
	pAtom1 = _ic.pAtom1;
	pAtom2 = _ic.pAtom2;
	pAtom3 = _ic.pAtom3;
	pAtom4 = _ic.pAtom4;

	vals[0]   = _ic.vals[0];
	vals[1]   = _ic.vals[1];
	vals[2] = _ic.vals[2];
	vals[3]   = _ic.vals[3];
	vals[4]   = _ic.vals[4];

	improperFlag = _ic.improperFlag;

	storedValues = _ic.storedValues;

}

void IcEntry::operator=(const IcEntry & _ic) {
	copy(_ic);
}

/*
bool IcEntry::build(Atom * _pAtom, map<Atom*, bool> & _exclude, bool _onlyActive) {
	/ ********************************************************
	 *   The ic can build either atom 1 or 4: the pointer tells
	 *   the function what atom we want to build
	 *
	 *   Used by the atom when it calls one of it own IcEntry,
	 *   since it does know if it is atom1 or atom4
	 ******************************************************** /
	if (pAtom1 == _pAtom) {
	        //MSLOUT.stream() << "Asked to build atom1, call build1"<<endl;
		return build1(_exclude, _onlyActive);
	} else if (pAtom4 == _pAtom) {
	        //MSLOUT.stream() << "Asked to build atom4, call build4"<<endl;
		return build4(_exclude, _onlyActive);
	} else {
		// the pointer is not pointing to atom 1 or 4
		return false;
	}
}

bool IcEntry::build1(map<Atom*, bool> & _exclude, bool _onlyActive) {
	if (pAtom1 == NULL) {
		return false;
	} else if (pAtom1->hasCoor()) {
		// already built
		return true;
	} else {
		// check that is a valid IC entry
		if (pAtom2 == NULL || pAtom3 == NULL || pAtom4 == NULL || vals[0] == 0.0 || vals[1] == 0.0) {
			return false;
		}
		if (_onlyActive && (!pAtom2->getActive() || !pAtom3->getActive() || !pAtom4->getActive())) {
			// some atoms are in inactive identities
			return false;
		}
		/ ********************************************************
		 *   Here the function calls a recursive build.  If 
		 *   pAtom2 3 and 4 already have coordinates, the buildFromIc()
		 *   returns true, else it will try to build them
		 *
		 *   For atom1 the build function is called differently if the
		 *   IC entry is an improper dihedral
		 *
		 *              1                  1
		 *               \                  \
		 *                2--3               3
		 *                    \             / \
		 *                     4           2   4
		 *               normal           improper
		 *
		 ******************************************************** /
		//_exclude.push_back(this);
		_exclude[pAtom1] = true;
	        //MSLOUT.stream() << "build1, try to build: "<<pAtom2->getName()<<" then "<<pAtom3->getName()<<" "<<pAtom4->getName()<<endl;
		if (pAtom2->buildFromIc(_exclude, _onlyActive) && pAtom3->buildFromIc(_exclude, _onlyActive) && pAtom4->buildFromIc(_exclude, _onlyActive)) {
			// atoms 2 3 and 4 have coordinates, build atom 1
			if (improperFlag) {
			 // MSLOUT.stream() << "build1, build atom1 improper dihedral : "<<pAtom1->getName()<<" with "<<-vals[2]<<endl;
				// improper dihedral, pass atoms as 3, 2, 4 and invert the sign of the dihedral
				pAtom1->setCoor(CartesianGeometry::buildRadians(pAtom3->getCoor(), pAtom2->getCoor(), pAtom4->getCoor(), vals[0], vals[1], -vals[2]));
			} else {

			  //MSLOUT.stream() << "build1, build atom1 normal dihedral : "<<pAtom1->getName()<<endl;
				// normal dihedral
				pAtom1->setCoor(CartesianGeometry::buildRadians(pAtom2->getCoor(), pAtom3->getCoor(), pAtom4->getCoor(), vals[0], vals[1], vals[2]));
			}
			return true;
		} else {
			cout << "      UUUI Unable to build1 " << pAtom4->getAtomId() << endl;
			return false;
		}
	}
}

bool IcEntry::build4(map<Atom*, bool> & _exclude, bool _onlyActive) {
	if (pAtom4 == NULL) {
		return false;
	} else if (pAtom4->hasCoor()) {
		// already built
		return true;
	} else {
		// check that is a valid IC entry
		if (pAtom3 == NULL || pAtom2 == NULL || pAtom1 == NULL || vals[4] == 0.0 || vals[3] == 0.0) {
			return false;
		}
		if (_onlyActive && (!pAtom3->getActive() || !pAtom2->getActive() || !pAtom1->getActive())) {
			// some atoms are in inactive identities
			return false;
		}
		/ ********************************************************
		 *   Here the function calls a recursive build.  If 
		 *   pAtom2 3 and 4 already have coordinates, the buildFromIc()
		 *   returns true, else it will try to build them
		 ******************************************************** /
		//_exclude.push_back(this);
		_exclude[pAtom4] = true;
		if (pAtom3->buildFromIc(_exclude, _onlyActive) && pAtom2->buildFromIc(_exclude, _onlyActive) && pAtom1->buildFromIc(_exclude, _onlyActive)) {
			pAtom4->setCoor(CartesianGeometry::buildRadians(pAtom3->getCoor(), pAtom2->getCoor(), pAtom1->getCoor(), vals[4], vals[3], vals[2]));
			return true;
		} else {
			cout << "      UUUI Unable to build4 " << pAtom4->getAtomId() << endl;
			return false;
		}
	}
}
*/

bool IcEntry::build(Atom * _pAtom, map<IcEntry*, bool> & _exclude, bool _onlyActive) {
	/********************************************************
	 *   The ic can build either atom 1 or 4: the pointer tells
	 *   the function what atom we want to build
	 *
	 *   Used by the atom when it calls one of it own IcEntry,
	 *   since it does know if it is atom1 or atom4
	 ********************************************************/
	if (pAtom1 == _pAtom) {
	        //MSLOUT.stream() << "Asked to build atom1, call build1"<<endl;
		return build1(_exclude, _onlyActive);
	} else if (pAtom4 == _pAtom) {
	        //MSLOUT.stream() << "Asked to build atom4, call build4"<<endl;
		return build4(_exclude, _onlyActive);
	} else {
		// the pointer is not pointing to atom 1 or 4
		return false;
	}
}

bool IcEntry::build1(map<IcEntry*, bool> & _exclude, bool _onlyActive) {
	if (pAtom1 == NULL) {
		return false;
	} else if (pAtom1->hasCoor()) {
		// already built
		return true;
	} else {
		// check that is a valid IC entry
		if (pAtom2 == NULL || pAtom3 == NULL || pAtom4 == NULL || vals[0] == 0.0 || vals[1] == 0.0) {
			return false;
		}
		if (_onlyActive && (!pAtom2->getActive() || !pAtom3->getActive() || !pAtom4->getActive())) {
			// some atoms are in inactive identities
			return false;
		}
		/********************************************************
		 *   Here the function calls a recursive build.  If 
		 *   pAtom2 3 and 4 already have coordinates, the buildFromIc()
		 *   returns true, else it will try to build them
		 *
		 *   For atom1 the build function is called differently if the
		 *   IC entry is an improper dihedral
		 *
		 *              1                  1
		 *               \                  \
		 *                2--3               3
		 *                    \             / \
		 *                     4           2   4
		 *               normal           improper
		 *
		 ********************************************************/
		//_exclude.push_back(this);
		///_exclude[pAtom1] = true;
		_exclude[this] = true;
	        //MSLOUT.stream() << "build1, try to build: "<<pAtom2->getName()<<" then "<<pAtom3->getName()<<" "<<pAtom4->getName()<<endl;
		if (pAtom2->buildFromIc(_exclude, _onlyActive) && pAtom3->buildFromIc(_exclude, _onlyActive) && pAtom4->buildFromIc(_exclude, _onlyActive)) {
			// atoms 2 3 and 4 have coordinates, build atom 1
			if (improperFlag) {
			 // MSLOUT.stream() << "build1, build atom1 improper dihedral : "<<pAtom1->getName()<<" with "<<-vals[2]<<endl;
				// improper dihedral, pass atoms as 3, 2, 4 and invert the sign of the dihedral
				pAtom1->setCoor(CartesianGeometry::buildRadians(pAtom3->getCoor(), pAtom2->getCoor(), pAtom4->getCoor(), vals[0], vals[1], -vals[2]));
			} else {

			  //MSLOUT.stream() << "build1, build atom1 normal dihedral : "<<pAtom1->getName()<<endl;
				// normal dihedral
				pAtom1->setCoor(CartesianGeometry::buildRadians(pAtom2->getCoor(), pAtom3->getCoor(), pAtom4->getCoor(), vals[0], vals[1], vals[2]));
			}
			return true;
		} else {
			//cout << "      UUUI Unable to build1 " << pAtom4->getAtomId() << endl;
			return false;
		}
	}
}

bool IcEntry::build4(map<IcEntry*, bool> & _exclude, bool _onlyActive) {
	if (pAtom4 == NULL) {
		return false;
	} else if (pAtom4->hasCoor()) {
		// already built
		return true;
	} else {
		// check that is a valid IC entry
		if (pAtom3 == NULL || pAtom2 == NULL || pAtom1 == NULL || vals[4] == 0.0 || vals[3] == 0.0) {
			return false;
		}
		if (_onlyActive && (!pAtom3->getActive() || !pAtom2->getActive() || !pAtom1->getActive())) {
			// some atoms are in inactive identities
			return false;
		}
		/********************************************************
		 *   Here the function calls a recursive build.  If 
		 *   pAtom2 3 and 4 already have coordinates, the buildFromIc()
		 *   returns true, else it will try to build them
		 ********************************************************/
		//_exclude.push_back(this);
		//_exclude[pAtom4] = true;
		_exclude[this] = true;
		if (pAtom3->buildFromIc(_exclude, _onlyActive) && pAtom2->buildFromIc(_exclude, _onlyActive) && pAtom1->buildFromIc(_exclude, _onlyActive)) {
			pAtom4->setCoor(CartesianGeometry::buildRadians(pAtom3->getCoor(), pAtom2->getCoor(), pAtom1->getCoor(), vals[4], vals[3], vals[2]));
			return true;
		} else {
			//cout << "      UUUI Unable to build4 " << pAtom4->getAtomId() << endl;
			return false;
		}
	}
}

string IcEntry::toString() const {
	string imp = " ";
	string A1 = "";
	string A2 = "";
	string A3 = "";
	string A4 = "";
	char c [1000];
	if (pAtom1!=NULL) {
		sprintf(c, "%1s %4d%1s %-4s", pAtom1->getChainId().c_str(), pAtom1->getResidueNumber(), pAtom1->getResidueIcode().c_str(), pAtom1->getName().c_str());
		A1 = c;
	}
	if (pAtom2!=NULL) {
		sprintf(c, "%1s %4d%1s %-4s", pAtom2->getChainId().c_str(), pAtom2->getResidueNumber(), pAtom2->getResidueIcode().c_str(), pAtom2->getName().c_str());
		A2 = c;
	}
	if (pAtom3!=NULL) {
		sprintf(c, "%1s %4d%1s %-4s", pAtom3->getChainId().c_str(), pAtom3->getResidueNumber(), pAtom3->getResidueIcode().c_str(), pAtom3->getName().c_str());
		A3 = c;
	}
	if (pAtom4!=NULL) {
		sprintf(c, "%1s %4d%1s %-4s", pAtom4->getChainId().c_str(), pAtom4->getResidueNumber(), pAtom4->getResidueIcode().c_str(), pAtom4->getName().c_str());
		A4 = c;
	}
	if (improperFlag) {
		imp = "*";
	}
	sprintf(c, "[%12s] [%12s] [%12s] [%12s] %8.3f %8.3f %8.3f %8.3f %8.3f [%1s]", A1.c_str(), A2.c_str(), A3.c_str(), A4.c_str(), vals[0], vals[1] * 180.0 / M_PI, vals[2] * 180.0 / M_PI, vals[3] * 180.0 / M_PI, vals[4], imp.c_str());
	return (string)c;
}

void IcEntry::fillFromCoor() {
//	vals[0] = 0.0;
//	vals[1] = 0.0;
//	vals[2] = 0.0;
//	vals[3] = 0.0;
//	vals[4] = 0.0;
	if (pAtom1 != NULL && pAtom1->hasCoor()) {
		if (improperFlag) {
			if (pAtom3 != NULL && pAtom3->hasCoor()) {
				vals[0] = pAtom1->distance(*pAtom3);
				if (pAtom2 != NULL && pAtom2->hasCoor()) {
					vals[1] = pAtom1->angleRadians(*pAtom3, *pAtom2);
				}
			}
		} else {
			if (pAtom2 != NULL && pAtom2->hasCoor()) {
				vals[0] = pAtom1->distance(*pAtom2);
				if (pAtom3 != NULL && pAtom3->hasCoor()) {
					vals[1] = pAtom1->angleRadians(*pAtom2, *pAtom3);
				}
			}
		}
		if (pAtom4 != NULL && pAtom4->hasCoor()) {
			vals[2] = pAtom1->dihedralRadians(*pAtom2, *pAtom3, *pAtom4);
		}
	} 
	if (pAtom4 != NULL && pAtom4->hasCoor()) {
		if (pAtom3 != NULL && pAtom3->hasCoor()) {
			vals[4] = pAtom4->distance(*pAtom3);
			if (pAtom2 != NULL && pAtom2->hasCoor()) {
				vals[3] = pAtom4->angleRadians(*pAtom3, *pAtom2);
			}
		}
	}
}

void IcEntry::setAtom1(Atom & _atom) {
	if (pAtom1 != NULL) {
		// remove from the old atom 1 the pointer to the IC
		pAtom1->removeIcEntry(this);
	}
	updateParentMap(pAtom1, &_atom);
	pAtom1 = &_atom;
}
void IcEntry::setAtom2(Atom & _atom) {
	if (pAtom2 != NULL) {
		pAtom2->removeIcEntry(this);
	}
	updateParentMap(pAtom2, &_atom);
	pAtom2 = &_atom;
}
void IcEntry::setAtom3(Atom & _atom) {
	if (pAtom3 != NULL) {
		pAtom3->removeIcEntry(this);
	}
	updateParentMap(pAtom3, &_atom);
	pAtom3 = &_atom;
}
void IcEntry::setAtom4(Atom & _atom) {
	if (pAtom4 != NULL) {
		pAtom4->removeIcEntry(this);
	}
	updateParentMap(pAtom4, &_atom);
	pAtom4 = &_atom;
}

void IcEntry::updateParentMap(Atom * _pOldAtom, Atom * _pNewAtom) {
	if (pParentTable != NULL) {
		pParentTable->updateMap(_pOldAtom, _pNewAtom);
	}
}

bool IcEntry::removeAtom(Atom * _pAtom) {
	bool found = false;
	if (_pAtom != NULL) {
		if (pAtom1 == _pAtom) {
			pAtom1 = NULL;
			found = true;
		}
		if (pAtom2 == _pAtom) {
			pAtom2 = NULL;
			found = true;
		}
		if (pAtom3 == _pAtom) {
			pAtom3 = NULL;
			found = true;
		}
		if (pAtom4 == _pAtom) {
			pAtom4 = NULL;
			found = true;
		}
	//	if (found) {
	//		updateParentMap(_pAtom, NULL);
	//	}
	}
	return found;
}
