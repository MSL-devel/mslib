/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#include "Position.h"
#include "Chain.h"

using namespace MSL;
using namespace std;


Position::Position() {
	setup(1, "", "A");
}
Position::Position(const std::string _positionId){
  string chainId;
  int resNum;
  string icode;
  MslTools::parsePositionId(_positionId,chainId,resNum,icode);
  setup(resNum,icode,chainId);
}
Position::Position(int _resNum, string _icode) {
	setup(_resNum, _icode, "A");
}

Position::Position(const AtomPointerVector & _atoms, string _resName, int _resNum, string _icode) {
	addIdentity(Residue(_atoms, _resName, residueNumber, residueIcode));
	setup(_resNum, _icode, "A");
}

Position::Position(const Residue & _residue, int _resNum, string _icode) {
	addIdentity(_residue);
	setup(_resNum, _icode, "A");
}

Position::Position(const Position & _position) {
	pParentChain = NULL;
	copy(_position);
}

Position::~Position() {
	deletePointers();
	linkedPositions.clear();
}

void Position::operator=(const Position & _position) {
	copy(_position);
	updateChainsActiveAtomList();
}

void Position::setup(int _resNum, string _insertionCode, string _chainId) {
	pParentChain = NULL;
	//addIdentity(_residue);
	currentIdentityIterator = identities.begin();
	//setActiveAtomsVector();
	residueNumber = _resNum;
	residueIcode = _insertionCode;
	chainId = _chainId;
	index = 0;
	foundIdentity = identityMap.end();
	positionType = UNLINKED;
}

void Position::copy(const Position & _position) {
	deletePointers();
	residueNumber = _position.residueNumber;
	residueIcode = _position.residueIcode;
	updateChainMap();
	chainId = _position.getChainId();
	for (vector<Residue*>::const_iterator k=_position.identities.begin(); k!=_position.identities.end(); k++) {
		addIdentity(**k);
	}
	currentIdentityIterator = identities.begin() + (_position.currentIdentityIterator - _position.identities.begin());
	setActiveAtomsVector();
	updateAllAtomsList();
	index = 0;
	foundIdentity = identityMap.end();
	identityReverseLookup.clear();

	// Correct behavior ?
	positionType = UNLINKED;
	linkedPositions.clear();
}

void Position::deletePointers() {
	for (vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
		delete *k;
		*k = NULL;
	}
	identities.clear();
	identityIndex.clear();
	currentIdentityIterator = identities.begin();
	identityMap.clear();
	activeAtoms.clear();
	activeAndInactiveAtoms.clear();
	identityReverseLookup.clear();
	updateChainsActiveAtomList();
}

void Position::addIdentity(vector<string> _atomNames,string _residueName){

	Residue *foo = new Residue();
	foo->setChainId(getChainId());
	foo->setResidueNumber(getResidueNumber());
	foo->setResidueIcode(getResidueIcode());
	foo->setResidueName(_residueName);


	for (uint i = 0; i < _atomNames.size();i++){
		foo->addAtom(_atomNames[i]);
	}

	addIdentity(*foo);

	delete(foo);
}

void Position::addIdentity(AtomPointerVector _atoms, string _name) {

	addIdentity(Residue(_atoms, _name, residueNumber, residueIcode));

}

void Position::addIdentity(const Residue & _residue) {

	string name = _residue.getResidueName();
	foundIdentity=identityMap.find(name);

	if (foundIdentity==identityMap.end()) {
		/******************************************
		 *  This identity DOES NOT exist: 
		 *   - add it
		 ******************************************/
		unsigned int iteratorPosition = 0;
		if (identities.size() > 0) {
			// take note of where the iterator to the current identity was
			// (necessary in case the internal size has changed)
			iteratorPosition = currentIdentityIterator - identities.begin();
		}

		identities.push_back(new Residue(_residue));
		identities.back()->setParentPosition(this);
		identityMap[name] = identities.back();
		identityIndex[identities.back()] = identities.size() - 1;

		// update the reverse lookup map
		identityReverseLookup.clear();
		for (map<string, Residue*>::iterator k=identityMap.begin(); k!=identityMap.end(); k++) {
			identityReverseLookup[k->second] = k;
		}


		// add the new atoms to the atom list
		activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), identities.back()->getAtomPointers().begin(), identities.back()->getAtomPointers().end());

		// restore the iterator 
		currentIdentityIterator = identities.begin() + iteratorPosition;
		if (identities.size() == 1) {
			// if it is the first identity update the current atom vector
			activeAtoms = (*currentIdentityIterator)->getAtomPointers();
			updateChainsActiveAtomList();
		}

	} else {
		/******************************************
		 *  This identity DOES ALREADY exist: 
		 *   - do nothing and return
		 ******************************************/
		 return;
	}
}		

bool Position::removeIdentity(string _name) {

	foundIdentity=identityMap.find(_name);

	if (foundIdentity!=identityMap.end()) {
		for (vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
			if ((*foundIdentity).second == *k) {
				// deallocate from memory and remove from list
				delete *k;
				identities.erase(k);
				// erase from the map
				identityMap.erase(foundIdentity);
				foundIdentity = identityMap.end();
				// update the reverse lookup map
				identityReverseLookup.clear();
				for (map<string, Residue*>::iterator k=identityMap.begin(); k!=identityMap.end(); k++) {
					identityReverseLookup[k->second] = k;
				}
				identityIndex.clear();
				for (vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
					identityIndex[*k] = k - identities.begin();
				}
				setActiveAtomsVector();
				updateAllAtomsList();
				updateChainsActiveAtomList();
				return true;
			}
		}
	}
	return false;
}

void Position::removeAllIdentities() {
	deletePointers();
	addIdentity(Residue("DUM", 1, ""));
	currentIdentityIterator = identities.begin();
	setActiveAtomsVector();
	updateAllAtomsList();
	updateChainsActiveAtomList();
}


void Position::addAtoms(const AtomPointerVector & _atoms) {
	/***********************************************
	 *  This function splits the AtomPointerVector in a number
	 *  of AtomPointerVector objects by residue name and then calls
	 *  the addAtoms(const AtomPointerVector & _atoms) function
	 *  of the residues to take care of the rest
	 ***********************************************/
	map<string, AtomPointerVector> dividedInIdentities;

	// store the order of the identities so that it will be preserved
	vector<string> identityOrder;
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		if (dividedInIdentities.find((*k)->getResidueName()) == dividedInIdentities.end()) {
			identityOrder.push_back((*k)->getResidueName());
		}
		dividedInIdentities[(*k)->getResidueName()].push_back(*k);
	}
	//cout << "UUU Position divided atom vector in " << dividedInIdentities.size() << " identities" << endl;

	bool callChainUpdate_flag = false;

	for (unsigned int i=0; i<identityOrder.size(); i++) {
		foundIdentity=identityMap.find(identityOrder[i]);
		if (foundIdentity!=identityMap.end()) {
			/***********************************************
			 *  A residue with the residue name already EXISTS, add
			 *  the atoms to it
			 ***********************************************/
			(*foundIdentity).second->addAtoms(dividedInIdentities[identityOrder[i]]);
			if (foundIdentity->second == *currentIdentityIterator) {
				// we are modifying the current identity
				// we need to update the activeAtoms list
				activeAtoms = (*currentIdentityIterator)->getAtomPointers();
				callChainUpdate_flag = true;
			}
		} else {
			/***********************************************
			 *  A residue with the residue name DOES NOT EXIST, 
			 *  create a new residue first and add
			 *  the atoms to it
			 *
			 *  k = iterator, pointer to element of map
			 *  *k = element of map
			 *  (*k).second  = an AtomPointerVector
			 *  *((*k).second.begin()) = Atom *
			 *
			 ***********************************************/
			unsigned int iteratorPosition = 0;
			if (identities.size() > 0) {
				// take note of where the iterator to the current identity was
				// (necessary in case the internal size has changed)
				iteratorPosition = currentIdentityIterator - identities.begin();
			}
			Atom * tmpAtom = *(dividedInIdentities[identityOrder[i]].begin());
			identities.push_back(new Residue(tmpAtom->getResidueName(), residueNumber, residueIcode));
			identities.back()->setParentPosition(this);
			identityMap[tmpAtom->getResidueName()] = identities.back();
			identities.back()->addAtoms(dividedInIdentities[identityOrder[i]]);
			// restore the iterator 
			currentIdentityIterator = identities.begin() + iteratorPosition;
			if (identities.size() == 1) {
				// if it is the first identity update the current atom vector
				activeAtoms = (*currentIdentityIterator)->getAtomPointers();
				callChainUpdate_flag = true;
				//identityIndex[identities.back()] = identities.size() - 1;
			}
			identityIndex[identities.back()] = identities.size() - 1;
		}
	}

	identityReverseLookup.clear();
	for (map<string, Residue*>::iterator k=identityMap.begin(); k!=identityMap.end(); k++) {
		identityReverseLookup[k->second] = k;
	}
	if (callChainUpdate_flag) {
		updateChainsActiveAtomList();
	}
	updateAllAtomsList();
}


void Position::setChainId(string _chainId) {
	if (pParentChain != NULL) {
		pParentChain->setChainId(_chainId);
	} else {
		chainId = _chainId;
	}
}

string Position::getChainId() const {
	if (pParentChain != NULL) {
		return pParentChain->getChainId();
	} else {
		return chainId;
	}
}

void Position::updateResidueMap(Residue * _residue) {
	// if a residue changes its name it calls this function
	// to update the map
	for (map<string, Residue*>::iterator k=identityMap.begin(); k!=identityMap.end(); k++) {
		if (k->second == _residue) {
			if (k->first != _residue->getResidueName()) {
				identityMap.erase(k);
				identityMap[_residue->getResidueName()] = _residue;
				break;
			}
		}
	}
	identityReverseLookup.clear();
	for (map<string, Residue*>::iterator k=identityMap.begin(); k!=identityMap.end(); k++) {
		identityReverseLookup[k->second] = k;
	}

}

void Position::updateChainMap(){
	/******************************************************
	 * Update the map in the chain that maps the
	 * resnum and icode to the position address
	 ******************************************************/
	if (pParentChain != NULL) {
		pParentChain->updatePositionMap(this);
	}
}

void Position::updateChainsActiveAtomList() {
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
	if (pParentChain != NULL) {
		pParentChain->updateIndexing();
	}
}

void Position::updateChainsAllAtomList() {
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
	if (pParentChain != NULL) {
		pParentChain->updateAllAtomIndexing();
	}
}


// the name space defines if the atoms and residues are named with
// PDB or charmm names (or whatever else)
void Position::setNameSpace(string _nameSpace) {
	if (pParentChain != NULL) {
		pParentChain->setNameSpace(_nameSpace);
	} else {
		nameSpace = _nameSpace;
	}
}

string Position::getNameSpace() const {
	if (pParentChain != NULL) {
		return pParentChain->getNameSpace();
	} else {
		return nameSpace;
	}
}

int Position::getIndex() const {
       cerr << "DEPRECATED Function Position::getIndex() has been replaced by Position::getIndexInChain() or Position::getIndexInSystem(), for now this will return getIndexInSystem()\n";

       return getIndexInSystem();
}
int Position::getIndexInChain() const {
    	if (pParentChain != NULL) {
		return pParentChain->getPositionIndex(this);
	}

	return 0;
}

int Position::getIndexInSystem() const {
	if (pParentChain != NULL) {
		return pParentChain->getPositionIndexInSystem(this);
	}
	return 0;
} 
int Position::getReverseIndexInChain() const {
	if (pParentChain != NULL) {
		return (getIndexInChain() - pParentChain->positionSize());
	}

	return 0;
}

int Position::getReverseIndexInSystem() const{
       if (pParentChain != NULL) {
	 return (pParentChain->getReversePositionIndexInSystem(this));
       }  	

       return 0;
}

System * Position::getParentSystem() const {
	if (pParentChain != NULL) {
		return pParentChain->getParentSystem();
	} else {
		return NULL;
	}
}

bool Position::copyCoordinatesOfAtoms(vector<string> _sourcePosNames, bool _atomsWithoutCoorOnly, vector<string> _targetPosNames, string _sourceIdentity, string _targetIdentity) {

	/*****************************************************
	 *  Copy coordinates from atoms of one identity to
	 *  those of other identity
	 *
	 *  If the target position names (_targetPosNames) aren't
	 *  given, assumes that the same name are used (i.e. N -> N)
	 *
	 *  If the source identity is not given, assumes the 
	 *  current identity
	 *
	 *  If the target identity is not given, assumes all
	 *  other identities
	 *
	 *  If no atoms are given, assumes all atoms with the
	 *  same name
	 *
	 *  If _atomsWithoutCoorOnly is true, it only copies
	 *  onto atoms that do not have active coordinates
	 *
	 *  For example, to copy the N HN CA C O atoms
	 *  from the current identity to all other identities
	 *  use
	 *     vector<string> names;
	 *     names.push_back("N");
	 *     names.push_back("HN");
	 *     names.push_back("CA");
	 *     names.push_back("C");
	 *     names.push_back("O");
	 *     copyCoordinatesOfAtoms(names);
	 ******************************************************/

	map<string, Residue*>::iterator sourceIdentityIt;
	if (_sourceIdentity == "") {
		if (identities.size() == 1) {
			// just one identity, nothing to do
			return true;
		} else {
			_sourceIdentity = (*currentIdentityIterator)->getResidueName();
		}
	}

	if (!identityExists(_sourceIdentity)) {
		return false;
	} else {
		sourceIdentityIt = foundIdentity;
	}

	if (_sourcePosNames.size() == 0) {
		// default the list to all atoms of the residue
		AtomPointerVector atoms = sourceIdentityIt->second->getAtomPointers();
		for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
			_sourcePosNames.push_back((*k)->getName());
		}
	}
	
	if (_targetPosNames.size() == 0) {
		_targetPosNames = _sourcePosNames;
	} else if (_targetPosNames.size() != _sourcePosNames.size()) {
		return false;
	}
	
	for (map<string, Residue*>::iterator targetIdentityIt=identityMap.begin(); targetIdentityIt!=identityMap.end(); targetIdentityIt++) {
		if (targetIdentityIt == sourceIdentityIt) {
			continue;
		} else if (_targetIdentity != "" && _targetIdentity != targetIdentityIt->second->getResidueName()) {
			continue;
		}

		for (unsigned int i=0; i<_sourcePosNames.size(); i++) {
			// for each atom to be copied
			if (sourceIdentityIt->second->atomExists(_sourcePosNames[i])) {
				// get the atom in the source
				Atom * pAtomSource = &(sourceIdentityIt->second->getLastFoundAtom());
				if (targetIdentityIt->second->atomExists(_targetPosNames[i])) {
					// get the atom in the target and copy the coordinates
					Atom * pAtomTarget = &(targetIdentityIt->second->getLastFoundAtom());
					if (!_atomsWithoutCoorOnly || !pAtomTarget->hasCoor()) {
						pAtomTarget->setCoor(pAtomSource->getCoor());
					}
				}
			}
		}
	}
	return true;
}



bool Position::setActiveIdentity(unsigned int _i, bool _applyToLinked) {

	if (_applyToLinked && positionType == SLAVE) {

	       // Make sure we have a master position
	       if (linkedPositions.size() == 0){
		    std::cerr << "ERROR 2342 Position::setActiveIdentity(), asking to set active identity to a SLAVE position that has no master! "<<(*this).toString()<<std::endl;
		    exit(2342);
	       }

		// put the master in charge to change them all
		return linkedPositions[0]->setActiveIdentity(_i, true);
	}
	if (_i >= identities.size()) {
		return false;
	}
	if (currentIdentityIterator != identities.begin() + _i) {
		currentIdentityIterator = identities.begin() + _i;
		setActiveAtomsVector();
	}
	// apply to all linked positions
	bool OK = true;
	if (_applyToLinked && positionType == MASTER) {
		for (unsigned int i=0; i<linkedPositions.size(); i++) {
			if (!linkedPositions[i]->setActiveIdentity(_i, false)) {
				OK = false;
			}
		}
	}
	return OK;
}
bool Position::setActiveIdentity(std::string _resName, bool _applyToLinked) {


       if (_applyToLinked && positionType == SLAVE) {

	     // Make sure we have a master position
	     if (linkedPositions.size() == 0){
	       std::cerr << "ERROR 2343 Position::setActiveIdentity(), asking to set active identity to a SLAVE position that has no master! "<<(*this).toString()<<std::endl;
	       exit(2343);
	     }

	     // put the master in charge to change them all
	     return linkedPositions[0]->setActiveIdentity(_resName, true);
	}


	bool OK = false;
	for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
		if ((*k)->getResidueName() == _resName) {
			currentIdentityIterator = k;
			setActiveAtomsVector();

			// dwkulp commented out 11/23/11 . 
			// First, it should only be executed with "_applyToLinked" is set
			// Second, we are doing this outside the loop already
			// Third, if _applyToLinked = false we can get to this loop as a SLAVE position, which I think is unintended.

			/*
			  // apply to all linked positions
			  for (unsigned int i=0; i<linkedPositions.size(); i++) {
			    linkedPositions[i]->setActiveIdentity(_resName);
			  }
			*/
			OK = true;
			break;
		} 
	}

	// apply to all linked positions
	if (OK && _applyToLinked && positionType == MASTER) {
		for (unsigned int i=0; i<linkedPositions.size(); i++) {
			if (!linkedPositions[i]->setActiveIdentity(_resName, false)) {
				OK = false;
			}
		}
	}
	return OK;
}


void Position::setActiveRotamer(unsigned int _index, bool _applyToLinked) {
	if (_applyToLinked && positionType == SLAVE) {
		// put the master in charge to change them all
		linkedPositions[0]->setActiveRotamer(_index, true);
		return;
	}
	unsigned int tot = 0;
	unsigned int prevTot = 0;
	//std::cout << "Num identities: "<<identities.size()<< " on pos: "<<(*this).toString()<<std::endl;
	for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
		prevTot = tot;
		//tot += (*k)->getNumberOfAltConformations();
		tot += (*k)->getNumberOfRotamers();
		//std::cout << "Number rotamers: "<<(*k)->getNumberOfRotamers()<<" "<<tot<<" "<<_index<<std::endl;
		if (tot > _index) {
			setActiveIdentity(k-identities.begin());
			(*k)->setActiveConformation(_index - prevTot);
			//std::cout << "Set ("<<getPositionId()<<") conformation to : "<<getCurrentIdentity().getResidueName()<<" rotamer: "<<_index-prevTot<<" absolute: "<<_index<<std::endl;
			break;
		}
	}
	if (_applyToLinked && positionType == MASTER) {
		// apply to all linked SLAVE positions
	  // std::cout << "This position has "<<linkedPositions.size()<<" linked positions"<<std::endl;
		for (unsigned int i=0; i<linkedPositions.size(); i++) {
		  //std::cout << "Set to linked position: "<<linkedPositions[i]->toString()<<std::endl;
			linkedPositions[i]->setActiveRotamer(_index, false);
		}
	}
}

void Position::setActiveRotamer(std::string _identity, unsigned int _n, bool _applyToLinked) {
	if (_applyToLinked && positionType == SLAVE) {
		// put the master in charge to change them all
		linkedPositions[0]->setActiveRotamer(_identity, _n, true);
		return;
	}
	if (identityExists(_identity)) {
		foundIdentity->second->setActiveConformation(_n);
		setActiveIdentity(identityIndex[foundIdentity->second]);
	}
	if (_applyToLinked && positionType == MASTER) {
		// apply to all linked SLAVE positions
		for (unsigned int i=0; i<linkedPositions.size(); i++) {
			linkedPositions[i]->setActiveRotamer(_identity, _n, false);
		}
	}
}
