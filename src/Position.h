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

		std::string getRotamerId(unsigned int _skip=0) const;

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

		bool removeIdentity(std::string _resName, bool _allowEmpty=false);
		//void removeAllIdentities();

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
		bool setActiveIdentity(unsigned int _i, bool _applyToLinked=true);
		bool setActiveIdentity(std::string _resName, bool _applyToLinked=true);
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

		void setRotamerSamplingLevel(std::string _label);
		bool defineRotamerSamplingLevels(std::map<std::string,std::map<std::string,unsigned int> >& _levels);

		unsigned int getTotalNumberOfRotamers() const;  // this returns the sum of the alt confs for all identities
		unsigned int getTotalNumberOfRotamers(unsigned int _index) const;  // this returns the number of the alt confs the i-th identity
		unsigned int getTotalNumberOfRotamers(std::string _identityId);  // this returns the number of the alt confs a given identity, i.e. "A,37,ILE"
		void setActiveRotamer(unsigned int _index, bool _applyToLinked=true);  // this sets the position to the identity and conformation given by the index of all alt conf at all positions
		void setActiveRotamer(std::string _identity, unsigned int _n, bool _applyToLinked=true);  // this sets the position to the identity and conformation given by the index of all alt conf at all positions


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

		bool copyCoordinatesOfAtoms(std::vector<std::string> _sourcePosNames=std::vector<std::string>(), bool _atomsWithoutCoorOnly=true, std::vector<std::string> _targePosNames=std::vector<std::string>(), std::string _sourceIdentity="", std::string _targetIdentity="");

		friend std::ostream & operator<<(std::ostream &_os, const Position & _pos)  {_os << _pos.toString(); return _os;};
		std::string toString() const;

		void renumberNoUpdate(int _residueNumber); // special function called by the chain to renumber all residues

		double getSasa() const;

		enum LinkedPositionType { UNLINKED=0, MASTER=1, SLAVE=2 };

		std::vector<Position *>& getLinkedPositions(); // gets list of index for each linked position in this system.
		bool addLinkedPosition(Position &_pos); // this makes the position a MASTER
		unsigned int getLinkedPositionType() const;
	//	void setLinkedPositionType(unsigned int _linkPositionType);

		/***************************************************
		 *  Saving coordinates to buffers:
		 *
		 *  coordinates can be saved to named buffers (string _coordName),
		 *  and copied back from them
		 *
		 *  The difference between save coordinates to a buffer, and 
		 *  having multiple alternate coor is that the saved coord 
		 *  are simply a buffer that can be restored
		 *
		 *  Coor can be saved to buffer with two different commands:
		 *    saveCoor:
		 *      - saveCoor saves ONLY the current coor
		 *      - when restored with applySavedCoor, a buffer created with
		 *        saveCoor will replace the CURRENT coorinate only
		 *    saveAltCoor:
		 *      - saveAltCoor saves ALL alternative coordinates and
		 *        also remembers what was the current coordinate
		 *      - when restored with the same applySavedCoor, a buffer
		 *        created with saveAltCoor will wipe off all alternative
		 *        cordinates and recreate the situation that was present
		 *        when the buffer was saved
		 *
		 *  More details in Atom.h
		 ***************************************************/
		void saveCoor(std::string _coordName);
		void saveAltCoor(std::string _coordName);
		bool applySavedCoor(std::string _coordName);
		void clearSavedCoor(std::string _coordName="");		

		/***************************************************
		 *  Ask the chain if this position is N- or C-terminal
		 ***************************************************/
		bool isPositionNterminal() const;
		bool isPositionCterminal() const;

		/***************************************************
		 *  IDENTITIES CAN BE TEMPORARILY HIDDEN
		 *
		 *
		 *  If hidden the identity will be still stored but 
		 *  it is like it is not present.
		 *
		 *  All changes will be updated on all linked positions.
		 *
		 *  HOW TO OPERATE
		 *  For example, let's say there are 5 identities loaded.
		 *  The absolute index is 0 to 4 independently if a
		 *  identity is hidden. The relative index only considers 
		 *  unhidden identities. When everything is active
		 *  they are the same:
		 *    Relative index: 0 1 2 3 4
		 *    Absative index: 0 1 2 3 4
		 * 
		 *  Hide identity 1 with absolute index (the Position will
		 *  behave like it had only 4 identitys)
		 *    pos.hideIdentityAbsIndex(1);
		 *    Rel: 0 1 2 3 4  >>  0 - 1 2 3  (4)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *    pos.hideIdentityAbsIndex(3);
		 *    Rel: 0 - 1 2 3  >>  0 - 1 - 2  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *               ^              ^
		 *
		 *  Hide by relative index
		 *    pos.hideIdentityRelIndex(1);
		 *    Rel: 0 - 1 - 2  >>  0 - - - 1  (2)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *             ^              ^
		 *
		 *  Unhide a identity
		 *    pos.unhideIdentityAbsIndex(3);
		 *    Rel: 0 - - - 1  >>  0 - - 1 2  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *             ^              ^
		 *
		 *  Hide all identitys but one (absolute undex)
		 *    pos.hideAllIdentitiesButOneAbsIndex(3)
		 *    Rel: 0 - - 1 2  >>  - - - 0 -  (1)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *               ^              ^
		 *
		 *  Unhide a identity
		 *    pos.unhideIdentityAbsIndex(1);
		 *    Rel: - - - 0 -  >>  - 0 - 1 -  (2)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *  Hide all identitys but one (relative undex)
		 *    pos.hideAllIdentitiesButOneRelIndex(0)
		 *    Rel: - 0 - 1 -  >>  - 0 - - -  (1)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *           ^              ^
		 *
		 *  Hide all identitys except the first 3
		 *    pos.hideAllIdentitiesButFirstN(3);
		 *    Rel: - 0 - - -  >>  0 1 2 - -  (3)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *         ^ ^ ^          ^ ^ ^
		 *
		 *  Unhide all identitys
		 *    pos.unhideAllIdentities()
		 *    Rel: 0 1 2 - -  >>  0 1 2 3 4  (5)
		 *    Abs: 0 1 2 3 4  >>  0 1 2 3 4  (5)
		 *         ^ ^ ^ ^ ^      ^ ^ ^ ^ ^
		 *
		 ***************************************************/

		/****************************************************
		* TODO: Implement the following interfaces 
		* NOTE: 1) Be careful to update the identityIndex everytime an identity is hidden or unhidden
		* 	2) The identityMap will always contain a residue irrespective of if it is hidden or not
		* bool hideIdentityAbsIndex(unsigned int _absoluteIndex); // hide identity based on absolute index
		* bool hideAllIdentitiesButOneRelIndex(unsigned int _keepThisIndex); // turns all identitys off except one, expressed as relative index
		* bool hideAllIdentitiesButOneAbsIndex(unsigned int _keepThisIndex); // turns all identitys off except one, expressed as absolute index
		* bool hideAllIdentitiesButFirstN(unsigned int _numberToKeepAbsIndex); // turns all identitys off except the first N, expressed as absolute index
		 ***************************************************/

		bool unhideIdentity(std::string _resName);
		bool unhideAllIdentities(); 

		bool hideAllIdentitiesButOne(std::string _resName);
		bool hideIdentities(std::vector<std::string> _resNames);
		bool hideIdentity(std::string _resName);

		bool getHidden(const Residue* _pRes) const; // is the residue hidden?



	private:


		void deletePointers();
		void setup(int _resNum, std::string _insertionCode, std::string _chainId);
		void copy(const Position & _position);
		void setActiveAtomsVector();
		void updateChainMap();
		void updateChainsActiveAtomList();
		void updateChainsAllAtomList();
		void updateAllAtomsList();
		void updateIdentityIndex();

		bool hideAltIdentity(unsigned int _absoluteIndex, unsigned int _relativeIndex, unsigned int _indexInHidden);
		bool unhideAltIdentity(Residue* _pRes, unsigned int _idx, unsigned int _hiddenIdx);
		bool unhideIdentityAbsIndex(unsigned int _absoluteIndex); // unhide a specific identity based on absolute index
		bool hideIdentityRelIndex(unsigned int _relativeIndex); // hide identity based on relative index

		unsigned int index; // index of the position in the System's std::vector<Position*>

		Chain * pParentChain;
		
		int residueNumber;
		std::string residueIcode;
		std::string chainId;
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		std::vector<Residue*> identities;
		std::vector<Residue*> hiddenIdentities;
		std::vector<unsigned int> hiddenIdentityIndeces;

		std::vector<Residue*>::iterator currentIdentityIterator;
		std::map<std::string, Residue*> identityMap;
		std::map<Residue*, std::map<std::string, Residue*>::iterator > identityReverseLookup;
		std::map<Residue*, unsigned int> identityIndex; // this is the relative index
		//  reverse std::map??
		AtomPointerVector activeAtoms;
		AtomPointerVector activeAndInactiveAtoms;
		std::map<std::string, Residue*>::iterator foundIdentity;

		std::vector<Position *> linkedPositions;
		unsigned int positionType;

		
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
	if (foundIdentity != identityMap.end() && identityIndex.find(foundIdentity->second) != identityIndex.end()) {
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
			if (foundIdentity != identityMap.end() && identityIndex.find(foundIdentity->second) != identityIndex.end() ) {
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
	if (foundIdentity != identityMap.end() && identityIndex.find(foundIdentity->second) != identityIndex.end()) {
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
		//out += (*k)->getNumberOfAltConformations();
		out += (*k)->getNumberOfRotamers();
	}
	return out;
}
inline unsigned int Position::getTotalNumberOfRotamers(unsigned int _index) const {
	//return identities[_index]->getNumberOfAltConformations();
	return identities[_index]->getNumberOfRotamers();
}
inline unsigned int Position::getTotalNumberOfRotamers(std::string _identityId) {
	if (identityExists(_identityId)) {
		//return foundIdentity->second->getNumberOfAltConformations();
		return foundIdentity->second->getNumberOfRotamers();
	}
	return 0;
}

inline unsigned int Position::getIdentityIndex(Residue * _pRes) {
	if(identityIndex.find(_pRes) != identityIndex.end()) {
		return identityIndex[_pRes];
	 } else {
		 return -1;
	 }
}
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
inline bool Position::addLinkedPosition(Position &_pos){
	/************************************************************
	 *     TODO:
	 *
	 *     Here we need to check that the linked positions have
	 *     the same number of identities and each identity the
	 *     same number of rotamers
	 ************************************************************/
	 if (positionType == SLAVE) {
		 std::cerr << "WARNING 39143: cannot add linked positions to a SLAVE position " << _pos << std::endl;
		 return false;
	 }
	positionType = MASTER;
	linkedPositions.push_back(&_pos);
	_pos.positionType = SLAVE;
	_pos.linkedPositions = linkedPositions;
	_pos.linkedPositions.insert(_pos.linkedPositions.begin(), this); // the first is the master
	for (unsigned int i=0; i<linkedPositions.size(); i++) {
		// add the position to all other slaves
		linkedPositions[i]->linkedPositions.push_back(&_pos);
	}
	return true;
}
inline unsigned int Position::getLinkedPositionType() const { return positionType; }
//inline void Position::setLinkedPositionType(unsigned int _lpt){ positionType = _lpt;}

inline std::string Position::getPositionId(unsigned int _skip) const {
	return MslTools::getPositionId(getChainId(), getResidueNumber(), getResidueIcode(), _skip);
}
inline std::string Position::getRotamerId(unsigned int _skip) const {
	// Gets the identity from the current identity
	// Gets the conformation from the 0th atom of the current identity
	//   Alternatively one could check all the atoms current conformation to make sure they are all the same...

	if ((*currentIdentityIterator)->size() == 0){
		return "";
	}else {
		return MslTools::getRotamerId(getChainId(),getResidueNumber(),getResidueIcode(),getResidueName(),(*currentIdentityIterator)->getAtom(0).getActiveConformation());
	}
}
inline void Position::saveCoor(std::string _coordName) {activeAndInactiveAtoms.saveCoor(_coordName);}
inline void Position::saveAltCoor(std::string _coordName) {activeAndInactiveAtoms.saveAltCoor(_coordName);}
inline bool Position::applySavedCoor(std::string _coordName) {return activeAndInactiveAtoms.applySavedCoor(_coordName);}
inline void Position::clearSavedCoor(std::string _coordName) {activeAndInactiveAtoms.clearSavedCoor(_coordName);}
inline void Position::setRotamerSamplingLevel(std::string _label) {
	// adopt a certain sampling level for all residues in the position
	for (std::vector<Residue*>::const_iterator k=identities.begin(); k!=identities.end(); k++) {
		(*k)->setRotamerSamplingLevel(_label);
	}

}
inline bool Position::defineRotamerSamplingLevels(std::map<std::string,std::map<std::string,unsigned int> >& _levels) {
	// define the number of rotamers at certain sampling levels for all residues in the position
	for(std::map<std::string,std::map<std::string,unsigned int> >::iterator lev = _levels.begin(); lev != _levels.end(); lev++) {
		for (std::vector<Residue*>::const_iterator k=identities.begin(); k!=identities.end(); k++) {
			std::string resName = (*k)->getResidueName();
			if(lev->second.find(resName) != lev->second.end()) {
				(*k)->defineRotamerSamplingLevel(lev->first,lev->second[resName]);
			}
		}
	}
	return true;
}
}

#endif
