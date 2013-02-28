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

#include "Residue.h"
#include "AtomGroup.h"

using namespace MSL;
using namespace std;


AtomGroup::AtomGroup() {
	pParentResidue = NULL;
	groupNumber = 0;
	cachedCenter = CartesianPoint(0.0, 0.0, 0.0);
	stamp = 1941992349; // random number
}

AtomGroup::AtomGroup(Residue * _pParentResidue) {
	pParentResidue = _pParentResidue;
	groupNumber = 0;
}

AtomGroup::~AtomGroup() {
}

/*
CartesianPoint AtomGroup::getGeometricCenter(unsigned int _stamp) {
	/ ************************************************************
	 * A trick for speed, to prevent to recalculate the same center over
	 * and over if the function is called on all atoms.
	 *
	 * A stamp is given and if the stamp is identical to the current stamp,
	 * the precalculated center from the previous call is given,
	 * otherwise it is calculated and the result is cached
	 ************************************************************ /
	if (stamp != _stamp) {
		stamp = _stamp;
		CartesianPoint out;

		for (unsigned int i=0; i<size(); i++) {
			out += (*this)[i]->getCoor();
		}
		cachedCenter = out/size();
	}
	return cachedCenter;
}
*/

void AtomGroup::setResidueName(string _resname) {
	if (pParentResidue != NULL) {
		pParentResidue->setResidueName(_resname);
	} else {
		residueName = _resname;
	}
}

string AtomGroup::getResidueName() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getResidueName();
	} else {
		return residueName;
	}
}

void AtomGroup::setResidueNumber(int _resnum) {
	if (pParentResidue != NULL) {
		pParentResidue->setResidueNumber(_resnum);
	} else {
		residueNumber = _resnum;
	}
}

int AtomGroup::getResidueNumber() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getResidueNumber();
	} else {
		return residueNumber;
	}
}

void AtomGroup::setResidueIcode(string _icode) {
	if (pParentResidue != NULL) {
		pParentResidue->setResidueIcode(_icode);
	} else {
		residueIcode = _icode;
	}
}

string AtomGroup::getResidueIcode() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getResidueIcode();
	} else {
		return residueIcode;
	}
}

void AtomGroup::setChainId(string _chainId) {
	if (pParentResidue != NULL) {
		pParentResidue->setChainId(_chainId);
	} else {
		chainId = _chainId;
	}
}

string AtomGroup::getChainId() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getChainId();
	} else {
		return chainId;
	}
}

// the name space defines if the atoms and residues are named with
// PDB or charmm names (or whatever else)
void AtomGroup::setNameSpace(string _nameSpace) {
	if (pParentResidue != NULL) {
		pParentResidue->setNameSpace(_nameSpace);
	} else {
		nameSpace = _nameSpace;
	}
}

string AtomGroup::getNameSpace() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getNameSpace();
	} else {
		return nameSpace;
	}
}

bool AtomGroup::getActive() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getActive();
	} else {
		return true;
	}
}
bool AtomGroup::getHidden() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getHidden();
	} else {
		return false;
	}
}

void AtomGroup::updateResidueMap(Atom * _atom) {
	if (pParentResidue != NULL) {
		pParentResidue->updateAtomMap(_atom);
	}
}

unsigned int AtomGroup::getGroupNumber() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getGroupNumber(this);
	} else {
		return groupNumber;
	}
}

Position * AtomGroup::getParentPosition() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getParentPosition();
	} else {
		return NULL;
	}
}

Chain * AtomGroup::getParentChain() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getParentChain();
	} else {
		return NULL;
	}
}

System * AtomGroup::getParentSystem() const {
	if (pParentResidue != NULL) {
		return pParentResidue->getParentSystem();
	} else {
		return NULL;
	}
}

unsigned int AtomGroup::getIdentityIndex() {
	if (pParentResidue != NULL) {
		return pParentResidue->getIdentityIndex();
	} else {
		return 0;
	}
}

bool AtomGroup::isPositionNterminal() const {
	if (pParentResidue != NULL) {
		return pParentResidue->isPositionNterminal();
	}
	return false;
}
bool AtomGroup::isPositionCterminal() const {
	if (pParentResidue != NULL) {
		return pParentResidue->isPositionCterminal();
	}
	return false;
}

