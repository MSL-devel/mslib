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


#ifndef POSITION_H
#define POSITION_H

#include <string>
#include <map>

#include "Residue.h"


namespace MSL { 
class Chain;
class System;

class Position {
	public:
		Position();
		Position(int _resNum, std::string _icode="");
		Position(const std::string _positionId);
		Position(const AtomPointerVector & _atoms, std::string _resName, int _resNum, std::string _icode="");
		Position(const Residue & _residue, int _resNum, std::string _icode="");
		Position(const Position & _position);
		~Position();

		void operator=(const Position & _position); // assignment

		std::string getPositionId(unsigned int _skip=0) const; // return "A,37", "A,37A" or "37" is skip is set to 1
	
		std::string getResidueName() const;
		void setResidueName(std::string _resname);

		void setResidueNumber(int _resnum);
		int getResidueNumber() const;

		void setResidueIcode(std::string _icode);
		std::string getResidueIcode() const;

		void setChainId(std::string _chainId);
		std::string getChainId() const;

		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		void setParentChain(Chain * _chain);
		Chain * getParentChain() const;
		System * getParentSystem() const;

		/* ADD RESIDUE IDENTITIES */
		void addIdentity(const Residue & _residue);
		void addIdentity(AtomPointerVector _atoms, std::string _name);
		void addIdentity(std::vector<std::string> _atomNames, std::string _residueName);

		void addAtoms(const AtomPointerVector & _atoms);

		bool removeIdentity(std::string _resName);
		void removeAllIdentities();

		/* UPDATES REQUESTED BY IDENTITIES */
		void updateResidueMap(Residue * _residue);


		unsigned int getIdentityIndex(Residue * _pRes);

		/***************************************************
		 *  Alternate identities
		 *
		 *  The current identity is obtained using an
		 *  iterator (currentIdentityIterator ) that points to
		 *  the address of the current identity (NOTE:
		 *  the iterator is actually a pointer to a pointer
		 *  of Residue)
		 *
		 *  Using an interator for the current idendity,
		 *  instead of a Residue pointer, allow to
		 *  keep track of the current identity without
		 *  having another variable (but makes understanding
		 *  the code a bit harder).
		 *
		 *  The current cartesian point is
		 *  *(*currentIdentityIterator)
		 *
		 ***************************************************/
	//	unsigned int size() const; // number of identities -- removed, substituted by identitySize()
		unsigned int identitySize() const; // number of identities
		unsigned int residueSize() const; // number of identities
		unsigned int atomSize() const; // number of atoms in the current identity
		unsigned int allAtomSize() const; // number of atoms including the inactive identities
		bool setActiveIdentity(unsigned int _i);
		bool setActiveIdentity(std::string _resName);
		int getActiveIdentity() const;
		size_t getNumberOfIdentities() const;

		Residue & operator()(unsigned int _index); // (n) returns the n-th identity
		Residue & operator()(std::string _residueId); // "ILE"
		Atom & operator[](unsigned int _index); // [n] returns the n-th atom of the active identity
		Atom & operator[](std::string _atomId); // "ILE,CA" or just "CA"
		Residue & getIdentity(unsigned int _index); 
		Residue & getIdentity(std::string _name); /* HOW DO WE HANDLE THE ERROR IF _name DOES NOT EXIST? */
		Residue & getResidue(unsigned int _index); // alias for getIdentity
		Residue & getResidue(std::string _name); // alias for getIdentity
		Residue & getCurrentIdentity();
		Residue & getLastFoundIdentity();
		Residue & getLastFoundResidue();
		AtomPointerVector & getAtomPointers(); // only active
		AtomPointerVector & getAllAtomPointers(); // all atoms, including the inactive
		Atom & getAtom(std::string _atomId); // get an atom from the active identity ("CA") or any identity "LEU,CA"
		Atom & getAtom(unsigned int _index); // get an atom from the active identity
		Atom & getAtom(std::string _identity, std::string _name);

		unsigned int getTotalNumberOfRotamers() const;  // this returns the sum of the alt confs for all identities
		void setActiveRotamer(unsigned int _index);  // this sets the position to the identity and conformation given by the index of all alt conf at all positions
		void setActiveRotamer(std::string _identity, unsigned int _n);  // this sets the position to the identity and conformation given by the index of all alt conf at all positions

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		// each position stores its absolute index in the System's std::vector<Position*>
		void setIndex(unsigned int _index);
		int getIndex() const; // cerr DEPRECATED message and returns index from System
		int getIndexInChain() const; // returns index from Chain
		int getIndexInSystem() const; // returns index from System
		int getReverseIndexInChain() const; // returns index from Chain
		int getReverseIndexInSystem() const; // returns index from System

		bool identityExists(std::string _identityId);
		bool residueExists(std::string _identityId); // alias for identityExists
		bool atomExists(std::string _atomId);// takes "CA" or "ILE,CA" ("A,37,ILE,CA" also works)
		bool atomExists(std::string _identity, std::string _name);// check in a specific identity
		//DEPRECATED exists functions
	//	bool exists(std::string _name);// check the existance of atom names in the current identity
	//	bool exists(std::string _name, std::string _identity);// check in a specific identity
		Atom & getLastFoundAtom();

		bool copyCoordinatesOfAtoms(std::vector<std::string> _sourcePosNames=std::vector<std::string>(), std::vector<std::string> _targePosNames=std::vector<std::string>(), std::string _sourceIdentity="", std::string _targetIdentity="");

		friend std::ostream & operator<<(std::ostream &_os, const Position & _pos)  {_os << _pos.toString(); return _os;};
		std::string toString() const;

		void renumberNoUpdate(int _residueNumber); // special function called by the chain to renumber all residues

		double getSasa() const;

		enum LinkedPositionType { UNLINKED=0, MASTER=1, SLAVE=2 };

		std::vector<Position *>& getLinkedPositions(); // gets list of index for each linked position in this system.
		void addLinkedPosition(Position &_pos);
		int getLinkedPositionType();
		void setLinkedPositionType(int _linkPositionType);

	private:

		void deletePointers();
		void setup(int _resNum, std::string _insertionCode, std::string _chainId);
		void copy(const Position & _position);
		void setActiveAtomsVector();
		void updateChainMap();
		void updateChainsActiveAtomList();
		void updateChainsAllAtomList();
		void updateAllAtomsList();


		unsigned int index; // index of the position in the System's std::vector<Position*>

		Chain * pParentChain;
		
		int residueNumber;
		std::string residueIcode;
		std::string chainId;
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		std::vector<Residue*> identities;
		std::vector<Residue*>::iterator currentIdentityIterator;
		std::map<std::string, Residue*> identityMap;
		std::map<Residue*, std::map<std::string, Residue*>::iterator > identityReverseLookup;
		std::map<Residue*, unsigned int> identityIndex;
		//  reverse std::map??
		AtomPointerVector activeAtoms;
		AtomPointerVector activeAndInactiveAtoms;
		std::map<std::string, Residue*>::iterator foundIdentity;
		
		std::vector<Position *> linkedPositions;
		int positionType;

		
};

// INLINED FUNCTIONS
inline std::string Position::getResidueName() const {return (*currentIdentityIterator)->getResidueName();};
inline void Position::setResidueName(std::string _resName) {return (*currentIdentityIterator)->setResidueName(_resName);};
inline void Position::setResidueNumber(int _resnum) {residueNumber = _resnum; updateChainMap();};
inline int Position::getResidueNumber() const {return residueNumber;};
inline void Position::setResidueIcode(std::string _icode) {residueIcode = _icode; updateChainMap();};
inline std::string Position::getResidueIcode() const {return residueIcode;};
inline void Position::setParentChain(Chain * _chain) {pParentChain = _chain;};
inline Chain * Position::getParentChain() const {return pParentChain;};
//inline unsigned int Position::size() const {std::cerr << "WARNING: using deprecated Position::size() function.  Use identitySize() instead" << std::endl; return identities.size();}; // number of identities
inline unsigned int Position::identitySize() const {return identities.size();}; // number of identities
inline unsigned int Position::residueSize() const {return identities.size();}; // number of identities
inline unsigned int Position::atomSize() const { return activeAtoms.size(); }; // number of atoms in the current identity
inline unsigned int Position::allAtomSize() const { return activeAndInactiveAtoms.size(); }; // number of atoms in the current identity
inline bool Position::setActiveIdentity(unsigned int _i) {
	if (_i >= identities.size()) {
		return false;
	}
	if (currentIdentityIterator != identities.begin() + _i) {
		currentIdentityIterator = identities.begin() + _i;
		setActiveAtomsVector();
	}
	return true;
}
inline bool Position::setActiveIdentity(std::string _resName) { for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) { if ((*k)->getResidueName() == _resName) { currentIdentityIterator = k; setActiveAtomsVector(); return true; } } return false; }
inline int Position::getActiveIdentity() const {return currentIdentityIterator - identities.begin();};
inline size_t Position::getNumberOfIdentities() const {return identities.size();};
inline Residue & Position::operator()(unsigned int _index) {return *identities[_index];}; // (n) returns the n-th identity
inline Residue & Position::operator()(std::string _residueId) {return getIdentity(_residueId);}
inline Atom & Position::operator[](unsigned int _index) {return getAtom(_index);}; // [n] returns the n-th atom of the active identity
inline Atom & Position::operator[](std::string _atomId) {
	return getAtom(_atomId);
} // "CA" or "ILE,CA"
inline Residue & Position::getIdentity(unsigned int _index) {return *identities[_index];}
inline Residue & Position::getIdentity(std::string _identityId) {
	//return *identityMap[_name];
	if (identityExists(_identityId)) {
		return *(foundIdentity->second);
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 53809: identity " << _identityId << " does not exist in residue at inline inline Residue & Position::getIdentity(string _identityId)" << std::endl;
		exit(53809);
	}
}
inline Residue & Position::getResidue(unsigned int _index) {return getIdentity(_index);} // alias for getIdentity
inline Residue & Position::getResidue(std::string _identityId) {return getIdentity(_identityId);}; // alias for getIdentity
inline Residue & Position::getCurrentIdentity() {return *(*currentIdentityIterator);}
inline Residue & Position::getLastFoundIdentity() {return *(foundIdentity->second);}
inline Residue & Position::getLastFoundResidue() {return getLastFoundIdentity();}
inline AtomPointerVector & Position::getAtomPointers() {return activeAtoms;}
inline AtomPointerVector & Position::getAllAtomPointers() {return activeAndInactiveAtoms;}
inline Atom & Position::getAtom(std::string _atomId) {
//	return (*currentIdentityIterator)->getAtom(_name);
	if (atomExists(_atomId)) {
		return foundIdentity->second->getLastFoundAtom();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 53812: atom " << _atomId << " does not exist in residue at inline Atom & Position::getAtom(string _atomId)" << std::endl;
		exit(53812);
	}
}
inline Atom & Position::getAtom(std::string _identity, std::string _name) {
//	return (*currentIdentityIterator)->getAtom(_name);
	if (atomExists(_identity, _name)) {
		return foundIdentity->second->getLastFoundAtom();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 53817: atom " << _identity << "," << _name << " does not exist in residue at inline Atom & Position::getAtom(string _identity, string _name)" << std::endl;
		exit(53817);
	}
}
inline Atom & Position::getAtom(unsigned int _index) {return (*(*currentIdentityIterator))[_index];}; // [n] returns the n-th atom of the active identity
inline void Position::setIndex(unsigned int _index) {index = _index;}
inline bool Position::identityExists(std::string _identityId) {
	if (_identityId == "") {
		// nothing given, use the current residue name
		_identityId = (*currentIdentityIterator)->getResidueName();
	}
	foundIdentity=identityMap.find(_identityId);
	if (foundIdentity != identityMap.end()) {
		return true;
	} else {
		// try to parse an identityId
		std::string chainid;
		int resnum;
		std::string icode;
		std::string identity;
		bool OK = MslTools::parseIdentityId(_identityId, chainid, resnum, icode, identity, 2);
		if (!OK) {
			// was the residue identity not specified ("A,37")
			OK = MslTools::parsePositionId(_identityId, chainid, resnum, icode, 1);
			if (OK) {
				// "A,37" use default residue name
				identity = (*currentIdentityIterator)->getResidueName();
			}
		}
		if(OK) {
			foundIdentity = identityMap.find(identity);
			if (foundIdentity != identityMap.end()) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
}
inline bool Position::residueExists(std::string _identityId) {return identityExists(_identityId);} // alias for identityExists
inline bool Position::atomExists(std::string _atomId) {
	// this accepts either "CA" or "ILE,CA", or even "A,37,ILE,CA" (the chain and resnum are ignored
	std::string chainid;
	int resnum;
	std::string icode;
	std::string identity;
	std::string atomName;
	bool OK = MslTools::parseAtomOfIdentityId(_atomId, chainid, resnum, icode, identity, atomName, 3);
	// if "CA" was given, identity will be = "" and the next function will return the atom for the current
	// identity
	if(OK) {
		return atomExists(identity, atomName);
	}
	return false;
}
inline bool Position::atomExists(std::string _identity, std::string _name) {
	if (_identity == "") {
		// identity not given, use the current
		_identity = (*currentIdentityIterator)->getResidueName();
	}

	foundIdentity = identityMap.find(_identity);
	if (foundIdentity != identityMap.end()) {
		return foundIdentity->second->atomExists(_name);
	}
	return false;
}
//inline bool Position::exists(std::string _name) {std::cerr << "DEPRECATED: Position::exists(string), use Position::atomExist(string)" << std::endl; return atomExists(_name);}
//inline bool Position::exists(std::string _name, std::string _identity) {std::cerr << "DEPRECATED: Position::exists(string), use Position::atomExist(string)" << std::endl; return atomExists(_identity, _name);}
inline Atom & Position::getLastFoundAtom() {return foundIdentity->second->getLastFoundAtom();}
//inline Atom & Position::getLastFoundAtom() {return (*currentIdentityIterator)->getLastFoundAtom();}
inline void Position::setActiveAtomsVector() {
	if (identities.size() > 0) {
		activeAtoms = (*currentIdentityIterator)->getAtomPointers();
		updateChainsActiveAtomList();
	}
}
inline void Position::wipeAllCoordinates() {for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {(*k)->wipeAllCoordinates();}}
inline void Position::updateAllAtomsList() {
	activeAndInactiveAtoms.clear();
	for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
		activeAndInactiveAtoms.insert(activeAndInactiveAtoms.end(), (*k)->getAtomPointers().begin(), (*k)->getAtomPointers().end());
	}
	updateChainsAllAtomList();
}
inline unsigned int Position::getTotalNumberOfRotamers() const {
	unsigned int out = 0;
	for (std::vector<Residue*>::const_iterator k=identities.begin(); k!=identities.end(); k++) {
		out += (*k)->getNumberOfAltConformations();
	}
	return out;
}
inline void Position::setActiveRotamer(unsigned int _index) {
	unsigned int tot = 0;
	unsigned int prevTot = 0;
	for (std::vector<Residue*>::iterator k=identities.begin(); k!=identities.end(); k++) {
		prevTot = tot;
		tot += (*k)->getNumberOfAltConformations();
		if (tot > _index) {
			setActiveIdentity(k-identities.begin());
			(*k)->setActiveConformation(_index - prevTot);
			return;
		}
	}
}

inline void Position::setActiveRotamer(std::string _identity, unsigned int _n) {
	if (identityExists(_identity)) {
		foundIdentity->second->setActiveConformation(_n);
		setActiveIdentity(identityIndex[foundIdentity->second]);
	}
 }

inline unsigned int Position::getIdentityIndex(Residue * _pRes) {return identityIndex[_pRes]; }
inline std::string Position::toString() const {
	std::stringstream ss;
	ss << getPositionId() <<  " [";
	//ss << getChainId() << " " << getResidueNumber() << getResidueIcode() << " [";
	unsigned int active = getActiveIdentity();
	for (std::vector<Residue*>::const_iterator k=identities.begin(); k!=identities.end(); k++) {
		ss << (*k)->getResidueName();
		if (identities.size() > 1 && k-identities.begin() == active) {
			ss << "*";
		}
		if (k == identities.end() - 1) {
			ss << "]";
		} else {
			ss << ", ";
		}
	}
	return ss.str();
}
inline void Position::renumberNoUpdate(int _resnum) {
	// this function is different than the setResidueNumber because it does not
	// cause updates.  It is called by the chain to renumber all the residue
	// numbers and do a final update
	residueNumber = _resnum;
}
inline double Position::getSasa() const {
	return (*currentIdentityIterator)->getSasa();
}



inline std::vector<Position *>& Position::getLinkedPositions() {
	return linkedPositions;
	
}
inline void Position::addLinkedPosition(Position &_pos){
	linkedPositions.push_back(&_pos);
}
inline int Position::getLinkedPositionType() { return positionType; }
inline void Position::setLinkedPositionType(int _lpt){ positionType = _lpt;}

inline std::string Position::getPositionId(unsigned int _skip) const {
	return MslTools::getPositionId(getChainId(), getResidueNumber(), getResidueIcode(), _skip);
}

}

#endif
