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

#include "CharmmTopologyResidue.h"
#include <stdio.h>

using namespace MSL;
using namespace std;


CharmmTopologyResidue::CharmmTopologyResidue() {
	setup("", false, 0.0, "", "");
}

CharmmTopologyResidue::CharmmTopologyResidue(string _name, bool _isPatch, double _charge, string _firstPatch, string _lastPatch) {
	setup(_name, _isPatch, _charge, _firstPatch, _lastPatch);
}

CharmmTopologyResidue::CharmmTopologyResidue(const CharmmTopologyResidue & _res) {
	copy(_res);
}

CharmmTopologyResidue::~CharmmTopologyResidue() {
	deletePointers();
}

void CharmmTopologyResidue::operator=(const CharmmTopologyResidue & _res) {
	copy(_res);
}

void CharmmTopologyResidue::deletePointers() {
	for (vector<TopolAtom*>::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
}

void CharmmTopologyResidue::reset() {
	deletePointers();
	atoms.clear();
	IcTable.clear();
	bonds.clear();
	angles.clear();
	dihedrals.clear();
	impropers.clear();
	donors.clear();
	acceptors.clear();
	deletes.clear();
}

void CharmmTopologyResidue::setup(string _name, bool _isPatch, double _charge, string _firstPatch, string _lastPatch) {
	name = _name;
	isPatch = _isPatch;
	charge = _charge;
	firstPatch = _firstPatch;
	lastPatch = _lastPatch;
}

void CharmmTopologyResidue::copy(const CharmmTopologyResidue & _res) {
	reset();
	name = _res.name;
	isPatch = _res.isPatch;
	charge = _res.charge;
	firstPatch = _res.firstPatch;
	lastPatch = _res.lastPatch;

	for (vector<TopolAtom*>::const_iterator k=_res.atoms.begin(); k!=_res.atoms.end(); k++) {
		TopolAtom * tmp = new TopolAtom;
		tmp->name = (*k)->name;
		tmp->type = (*k)->type;
		tmp->partialCharge = (*k)->partialCharge;
		tmp->element = (*k)->element;
		tmp->group = (*k)->group;
		atoms.push_back(tmp);
	}
	IcTable = _res.IcTable;
	bonds = _res.bonds;
	angles = _res.angles;
	dihedrals = _res.dihedrals;
	impropers = _res.impropers;
	donors = _res.donors;
	acceptors = _res.acceptors;
	deletes = _res.deletes;

}

void CharmmTopologyResidue::unlinkAtoms() {
	// IC TABLE
	for (vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {
		k->pAtoms.clear();
	}

	// BONDS
	for (vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {
		k->pAtoms.clear();
	}

	// ANGLES
	for (vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {
		k->pAtoms.clear();
	}

	// DIHEDRALS
	for (vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {
		k->pAtoms.clear();
	}

	// IMPROPERS
	for (vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {
		k->pAtoms.clear();
	}

	// DONORS
	for (vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {
		k->pAtoms.clear();
	}

	// ACCEPTORS
	for (vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {
		k->pAtoms.clear();
	}
}

void CharmmTopologyResidue::linkAtoms() {
	// IC TABLE
	for (vector<IcLine>::iterator k=IcTable.begin(); k!=IcTable.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// BONDS
	for (vector<Bond>::iterator k=bonds.begin(); k!=bonds.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// ANGLES
	for (vector<Angle>::iterator k=angles.begin(); k!=angles.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// DIHEDRALS
	for (vector<Dihedral>::iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// IMPROPERS
	for (vector<Improper>::iterator k=impropers.begin(); k!=impropers.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// DONORS
	for (vector<Donor>::iterator k=donors.begin(); k!=donors.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}

	// ACCEPTORS
	for (vector<Acceptor>::iterator k=acceptors.begin(); k!=acceptors.end(); k++) {
		k->pAtoms.clear();
		for (vector<string>::iterator l=k->atoms.begin(); l!=k->atoms.end(); l++) {
			CharmmTopologyResidue * pRes = NULL;
			if (l->substr(0,1) == "-") {
				pRes = pPrevTopologyRes;
				// call previous
			} else if (l->substr(0,1) == "+") {
				// call next
				pRes = pNextTopologyRes;
			}
			if (pRes == NULL) {
				k->pAtoms.push_back(NULL);
			} else {
				map<string, TopolAtom*>::iterator found = pRes->atomMap.find(*l);
				if (found != atomMap.end()) {
					k->pAtoms.push_back(found->second);
				} else {
					k->pAtoms.push_back(NULL);
				}
			}
		}
	}
	
}

bool CharmmTopologyResidue::applyPatch(const CharmmTopologyResidue & _patch) {
	if (!_patch.isPatch) {
		return false;
	}

	// take care of the deletions
	for (vector<Delete>::const_iterator k=_patch.deletes.begin(); k!=_patch.deletes.end(); k++) {
		if (k->type == "ATOM") {
			//cout << "UUU deleting atom " << k->atoms[0] << endl;
			deleteTopolAtom(k->atoms[0]);
		} else if (k->type == "BOND") {
			//cout << "UUU deleting bonds " << k->atoms[0] << " " << k->atoms[1] << endl;
			unsigned int index = 0;
			while (k->atoms.size() > index+1) {
				deleteBond(k->atoms[index], k->atoms[index+1]);
				index += 2;
			}
		} else if (k->type == "ANGL" || k->type == "THET") {
			//cout << "UUU deleting angles " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << endl;
			unsigned int index = 0;
			while (k->atoms.size() > index+2) {
				deleteAngle(k->atoms[index], k->atoms[index+1], k->atoms[index+2]);
				index += 3;
			}
		} else if (k->type == "DIHE" || k->type == "PHI") {
			//cout << "UUU deleting dihedral " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << " " << k->atoms[3] << endl;
			unsigned int index = 0;
			while (k->atoms.size() > index+3) {
				deleteDihedral(k->atoms[index], k->atoms[index+1], k->atoms[index+2], k->atoms[index+3]);
				index += 4;
			}
		} else if (k->type == "IMPR" || k->type == "IMPH") {
			//cout << "UUU deleting improper " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << " " << k->atoms[3] << endl;
			unsigned int index = 0;
			while (k->atoms.size() > index+3) {
				deleteImproper(k->atoms[index], k->atoms[index+1], k->atoms[index+2], k->atoms[index+3]);
				index += 4;
			}
		} else if (k->type == "DONO") {
			if (k->atoms.size() == 4) {
				//cout << "UUU deleting donor " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << " " << k->atoms[3] << endl;
				// HERE VARIABLE NUMBER 1 to 4!!!!
				deleteDonor(k->atoms[0], k->atoms[1], k->atoms[2], k->atoms[3]);
			} else if (k->atoms.size() == 3) {
				//cout << "UUU deleting donor " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << endl;
				// HERE VARIABLE NUMBER 1 to 4!!!!
				deleteDonor(k->atoms[0], k->atoms[1], k->atoms[2]);
			} else if (k->atoms.size() == 2) {
				cout << "UUU deleting donor " << k->atoms[0] << " " << k->atoms[1] << endl;
				// HERE VARIABLE NUMBER 1 to 4!!!!
				deleteDonor(k->atoms[0], k->atoms[1]);
			} else {
				cout << "UUU deleting donor " << k->atoms[0] << endl;
				// HERE VARIABLE NUMBER 1 to 3!!!!
				deleteDonor(k->atoms[0]);
			}
		} else if (k->type == "ACCE") {
			if (k->atoms.size() == 3) {
				//cout << "UUU deleting acceptor " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << endl;
				deleteAcceptor(k->atoms[0], k->atoms[1], k->atoms[2]);
			} else if (k->atoms.size() == 2) {
				//cout << "UUU deleting acceptor " << k->atoms[0] << " " << k->atoms[1] << endl;
				deleteAcceptor(k->atoms[0], k->atoms[1]);
			} else {
				deleteAcceptor(k->atoms[0]);
			}
		} else if (k->type == "IC") {
			//cout << "UUU deleting ic " << k->atoms[0] << " " << k->atoms[1] << " " << k->atoms[2] << " " << k->atoms[3] << endl;
			deleteIcLine(k->atoms[0], k->atoms[1], k->atoms[2], k->atoms[3]);
		}
	}

	// check the new atoms, if they are already present in an existing group then
	// assign those to the same group.  All atoms 
	map<int, int> groupMap;
	vector<int> newGroups(_patch.atoms.size(), -1);
	map<unsigned int, unsigned int> patchToResAtomMap;
	for (vector<TopolAtom*>::const_iterator k=_patch.atoms.begin(); k!=_patch.atoms.end(); k++) {
		for (vector<TopolAtom*>::iterator l=atoms.begin(); l!=atoms.end(); l++) {
			if ((*k)->name == (*l)->name) {
				groupMap[(*k)->group] = (*l)->group;
				newGroups[k-_patch.atoms.begin()] = (*l)->group;
				patchToResAtomMap[k-_patch.atoms.begin()] = l-atoms.begin();
			}
		}
	}
	// look for highest group
	int maxGroup = 0;
	for (vector<TopolAtom*>::iterator l=atoms.begin(); l!=atoms.end(); l++) {
		if (l==atoms.begin() || maxGroup < (*l)->group) {
			maxGroup = (*l)->group;
		}
	}
	groupMap[-1] = maxGroup;
	int newGroupNum = maxGroup + 1;
	for (vector<TopolAtom*>::const_iterator k=_patch.atoms.begin(); k!=_patch.atoms.end(); k++) {
		if (newGroups[k-_patch.atoms.begin()] > -1) {
			// this atom was found in the residue, the group is assigned
			continue;
		} else if ((*k)->group == -1 && newGroups[k-_patch.atoms.begin()] == -1) {
			// this atom was not in a group in the patch (-1) and 
			// the atom was not found in the residue:
			// assign it to the last group
			newGroups[k-_patch.atoms.begin()] = maxGroup;
		} else if ((*k)->group > -1 && newGroups[k-_patch.atoms.begin()] == -1 && groupMap.find((*k)->group) == groupMap.end()) {
			// this atom was not found and it belong to a group in the patch, and
			// no other atom in the group was found: add a new group
			newGroups[k-_patch.atoms.begin()] = newGroupNum;
			groupMap[(*k)->group] = newGroupNum;
			newGroupNum++;
		} else if (groupMap.find((*k)->group) != groupMap.end()) {
			// if the group is mapped assign it to it
			newGroups[k-_patch.atoms.begin()] = groupMap[(*k)->group];
		} else {
			// if everything fails assign to last group
			newGroups[k-_patch.atoms.begin()] = maxGroup;
		}
	}

	for (vector<TopolAtom*>::const_iterator k=_patch.atoms.begin(); k!=_patch.atoms.end(); k++) {
		map<unsigned int, unsigned int>::iterator found = patchToResAtomMap.find(k-_patch.atoms.begin());
		if (found != patchToResAtomMap.end()) {
			// exsting atom: edit its type and charge
			//cout << "UUU atom " << (*k)->name << " existed in residue, updating it" << endl;
			(atoms[found->second])->type = (*k)->type;
			(atoms[found->second])->partialCharge = (*k)->partialCharge;
		} else {
			//cout << "UUU atom " << (*k)->name << " did not exist in residue, adding it to group " << newGroups[k-_patch.atoms.begin()] << endl;
			addTopolAtom((*k)->name, (*k)->type, (*k)->partialCharge, (*k)->element, newGroups[k-_patch.atoms.begin()]);
		}
	}

	for (vector<IcLine>::const_iterator k=_patch.IcTable.begin(); k!=_patch.IcTable.end(); k++) {
		IcTable.push_back(*k);
	}

	for (vector<Bond>::const_iterator k=_patch.bonds.begin(); k!=_patch.bonds.end(); k++) {
		bonds.push_back(*k);
	}

	for (vector<Angle>::const_iterator k=_patch.angles.begin(); k!=_patch.angles.end(); k++) {
		angles.push_back(*k);
	}

	for (vector<Dihedral>::const_iterator k=_patch.dihedrals.begin(); k!=_patch.dihedrals.end(); k++) {
		dihedrals.push_back(*k);
	}

	for (vector<Improper>::const_iterator k=_patch.impropers.begin(); k!=_patch.impropers.end(); k++) {
		impropers.push_back(*k);
	}

	for (vector<Donor>::const_iterator k=_patch.donors.begin(); k!=_patch.donors.end(); k++) {
		donors.push_back(*k);
	}

	for (vector<Acceptor>::const_iterator k=_patch.acceptors.begin(); k!=_patch.acceptors.end(); k++) {
		acceptors.push_back(*k);
	}
	sortAtomsByGroup();

	return true;
}

string CharmmTopologyResidue::toString() const {
	char line [1000];
	stringstream ss;
	sprintf(line, "RESI %-4s       %5.2f", name.c_str(), charge);
	ss << line << endl;
	int group = -1;
	for (vector<TopolAtom*>::const_iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (group != (*k)->group) {
			ss << "GROUP" << endl;
			group = (*k)->group;
		}
		cout << "UUUR >" << (*k)->name.c_str() << "< >" << (*k)->type.c_str() << "< >" << (*k)->partialCharge << "<" << endl;
		sprintf(line, "ATOM %-4s %-4s  %5.2f", (*k)->name.c_str(), (*k)->type.c_str(), (*k)->partialCharge);
		ss << line << endl;
	}
	for (vector<Bond>::const_iterator k=bonds.begin(); k!=bonds.end(); k++) {
		if (k->type == 1) {
			ss << "BOND " << k->atoms[0] << " " <<  k->atoms[1] << endl;
		} else if (k->type == 2) {
			ss << "DOUB " << k->atoms[0] << " " <<  k->atoms[1] << endl;
		} else if (k->type == 3) {
			ss << "TRIP " << k->atoms[0] << " " <<  k->atoms[1] << endl;
		}
	}
	for (vector<Angle>::const_iterator k=angles.begin(); k!=angles.end(); k++) {
		ss << "ANGL " << k->atoms[0] << " " <<  k->atoms[1] << " " <<  k->atoms[2] << endl;
	}
	for (vector<Dihedral>::const_iterator k=dihedrals.begin(); k!=dihedrals.end(); k++) {
		ss << "DIHE " << k->atoms[0] << " " <<  k->atoms[1] << " " <<  k->atoms[2] << " " <<  k->atoms[3] << endl;
	}
	for (vector<Improper>::const_iterator k=impropers.begin(); k!=impropers.end(); k++) {
		ss << "IMPR " << k->atoms[0] << " " <<  k->atoms[1] << " " <<  k->atoms[2] << " " <<  k->atoms[3] << endl;
	}
	for (vector<IcLine>::const_iterator k=IcTable.begin(); k!=IcTable.end(); k++) {
		if (k->improperFlag) {
			sprintf(line, "IC %-4s %-4s *%-4s %-4s  %8.3f %8.3f %8.3f %8.3f %8.3f", k->atoms[0].c_str(), k->atoms[1].c_str(), k->atoms[2].c_str(), k->atoms[3].c_str(), k->values[0], k->values[1], k->values[2], k->values[3], k->values[4]);
		} else {
			sprintf(line, "IC %-4s %-4s  %-4s %-4s  %8.3f %8.3f %8.3f %8.3f %8.3f", k->atoms[0].c_str(), k->atoms[1].c_str(), k->atoms[2].c_str(), k->atoms[3].c_str(), k->values[0], k->values[1], k->values[2], k->values[3], k->values[4]);
		}
		ss << line << endl;
	}
	ss << "PATCHING FIRST " << firstPatch << " LAST " << lastPatch << endl;
	return ss.str();

}

void CharmmTopologyResidue::sortAtomsByGroup() {
	vector<TopolAtom*> newOrder;
	for (vector<TopolAtom*>::iterator k=atoms.begin(); k!=atoms.end(); k++) {
//		cout << (*k)->name << endl;
		bool found = false;
		for (vector<TopolAtom*>::iterator l=newOrder.begin(); l!=newOrder.end(); l++) {
			if ((*k)->group < (*l)->group) {
				newOrder.insert(l, k, k+1);
				found = true;
				break;
			}
		}
		if (!found) {
			newOrder.insert(newOrder.end(), k, k+1);
		}
	}
	atoms = newOrder;
//	for (vector<TopolAtom*>::iterator k=atoms.begin(); k!=atoms.end(); k++) {
//		cout << (*k)->name << endl;
//	}
}




vector<string> CharmmTopologyResidue::getAllTopoAtomNames(){
	vector<string> result;
	for (uint i = 0; i < atoms.size();i++){
		result.push_back(atoms[i]->name);
	}

	return result;
}
