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


#ifndef CHAIN_H
#define CHAIN_H

#include <map>
#include <string>
#include <vector>
#include "MslExceptions.h"

#include "MslTools.h"
#include "Position.h"


namespace MSL { 
class System;


class Chain {
	public:
		Chain();
		Chain(std::string _chainId);
		Chain(const std::vector<Residue> & _residues, std::string _chainId);
		Chain(const AtomPointerVector & _atomss, std::string _chainId);
		Chain(const Chain & _chain);
		~Chain();
		
		void operator=(const Chain & _chain); // assignment
		
		void setChainId(std::string _chainId);
		std::string getChainId() const;

		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;


		void addResidue(AtomPointerVector _atoms, std::string _name);
		void addResidue(AtomPointerVector _atoms, std::string _name, unsigned int _resNum, std::string _iCode="");
		void addResidue(const Residue & _residue);
		void addResidue(const Residue & _residue, unsigned int _resNum, std::string _iCode="");
		void addResidues(const std::vector<Residue> & _residues);

		bool removeResidue(int _resNum, std::string _iCode="");
		void removeAllResidues();

		void addAtoms(const AtomPointerVector & _atoms, bool _keepOrder=false);

		bool addIdentityToPosition(AtomPointerVector _atoms, std::string _name, unsigned int _resNum, std::string _iCode="");
		bool addIdentityToPosition(const Residue & _residue, unsigned int _resNum, std::string _iCode="");

		void setParentSystem(System * _system);
		System * getParentSystem() const;


	//	unsigned int size() const; // number of positions -- removed, substituted by positionSize()
		unsigned int positionSize() const; // number of positions
		std::vector<Position*> & getPositions();
		Position & getPosition(unsigned int _index);
//		Position & getPositionByIndex(unsigned int _index);
		//Position & getPositionByIndexTmp(unsigned int _index);
		Position & getPosition(std::string _positionId);
		Position & getPosition(int _resNum, std::string _iCode);

		int getPositionIndex(const Position * _pPos) const;  // CHANGED, returns index in Chain
		int getPositionIndexInSystem(const Position * _pPos) const;  // returns index in System

		int getReversePositionIndex(const Position * _pPos) const;  // CHANGED, returns index-size() in Chain
		int getReversePositionIndexInSystem(const Position * _pPos) const;  // returns index-size()- in System

	//	Residue & operator()(std::string _resNumAndIcode); // returns the active residue at position _resNumAndIcode (enter as "75" or "75A")
//		Residue & operator()(int _resNum); // same as above, but it takes and int and assumes the insertion code to be blank (i.e. 75)
		Position & operator()(std::string _positionId); // returns the active residue at position (enter as "75" or "75A")
		Position & operator()(unsigned int _index); // return the i-th position in the molecule (not the residue number)
		//Residue & getResidueByIndex(unsigned int _index);
		//Residue & getResidueByIndexTmp(unsigned int _index);
		Residue & getIdentity(unsigned int _index);
		Residue & getIdentity(std::string _identityId);
		Residue & getIdentity(int _resNum, std::string _iCode, std::string _identity="");
		// wrapper functions, same of getIdentity
		Residue & getResidue(unsigned int _index);
		Residue & getResidue(std::string _identityId);
		Residue & getResidue(int _resNum, std::string _iCode, std::string _identity="");
		AtomPointerVector & getAtomPointers();
		AtomPointerVector & getAllAtomPointers();
		unsigned int atomSize();
		unsigned int allAtomSize();
		Atom & operator[](unsigned int _index);
		Atom & operator[](std::string _atomId); // "A,37,CA" or "A,37,ILE,CA"
		Atom & getAtom(unsigned int _index);
		Atom & getAtom(int _resnum, std::string _icode, std::string _name);
		Atom & getAtom(int _resnum, std::string _icode, std::string _identity, std::string _name);
		Atom & getAtom(std::string _atomId);
		
		// check the existance of chains, residues, atoms
		bool positionExists(int _resNumm, std::string _icode); // residue, by int
		bool positionExists(std::string _positionId); // residue by string (possibly with insertion code
		bool identityExists(std::string _identityId); // "37,ILE" "37A,ILE"
		bool identityExists(int _resnum, std::string _icode, std::string _identity="");
		bool residueExists(std::string _identityId); // "37,ILE" "37A,ILE"
		bool residueExists(int _resnum, std::string _icode, std::string _identity="");
		bool atomExists(std::string _atomId); // atom
		bool atomExists(int _resNum, std::string _icode, std::string _name); // atom
		bool atomExists(int _resNum, std::string _icode, std::string _identity, std::string _name); // atom
		// DEPRECATED EXISTS FUNCTIONS
	//	bool exists(int _resNum); // position, by int
	//	bool exists(std::string _resNumAndIcode); // position by string (possibly with insertion code
	//	bool exists(int _resNum, std::string _name); // atom
	//	bool exists(std::string _resNumAndIcode, std::string _name);
	//	bool exists(int _resNum, std::string _name, std::string _identity); // atom w/ specified identity
	//	bool exists(std::string _resNumAndIcode, std::string _name, std::string _identity);

		Position & getLastFoundPosition();
		Residue & getLastFoundIdentity();
		Residue & getLastFoundResidue();
		Atom & getLastFoundAtom();

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates
		bool defineRotamerSamplingLevels(std::map<std::string,std::map<std::string,unsigned int> >& _levels);

		/* UPDATES REQUESTED BY POSITIONS */
		void updatePositionMap(Position * _position);
		void updateIndexing();
		void updateAllAtomIndexing();
	//	void swapInActiveList(Position * _position, AtomPointerVector & _atoms);
	
		/* RENUMBER THE WHOLE CHAIN */
		void renumberChain(int _start);


		std::string toString() const;
		friend std::ostream & operator<<(std::ostream &_os, const Chain & _chain)  {_os << _chain.toString(); return _os;};
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
		 *  Ask if a position is N- or C-terminal
		 ***************************************************/
		bool isPositionNterminal(const Position * _pPos) const;
		bool isPositionCterminal(const Position * _pPos) const;

	private:

		void deletePointers();
		void setup(std::string _chainId);
		void copy(const Chain & _chain);
		void updateSystemMap();
		void updateSystemActiveAtomList();
		void updateSystemAllAtomList();

		std::string chainId;
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		std::vector<Position*> positions;
		/*********************************************
		 *  The following is a 2D std::map [int residue number][std::string insertion code] 
		 *********************************************/
		std::map<int, std::map<std::string, Position*> > positionMap;

		/*********************************************
		 *  We keep a list of active atoms
		 *
		 *  We also keep an index of where the residue
		 *  atoms start and end so that we can quickly
		 *  swap them when a position changes identity
		 *********************************************/
		AtomPointerVector activeAtoms;
		AtomPointerVector activeAndInactiveAtoms;
		bool noUpdateIndex_flag;
		//struct ResidueAtoms {
		//	unsigned int start;
		//	unsigned int size;
		//};
		//std::map<Position *, ResidueAtoms> residueLookupMap;

		System * pParentSystem;

		std::map<std::string, Position*>::iterator foundPosition;

		
};

// INLINED FUNCTIONS
inline void Chain::setChainId(std::string _chainId) {chainId = _chainId; updateSystemMap();}
inline std::string Chain::getChainId() const {return chainId;}
inline void Chain::setParentSystem(System * _system) {pParentSystem = _system;}
inline System * Chain::getParentSystem() const {return pParentSystem;}
//inline unsigned int Chain::size() const {std::cerr << "WARNING: using deprecated Chain::size() function.  Use positionSize() instead" << std::endl; return positions.size();}
inline unsigned int Chain::positionSize() const {return positions.size();}
inline std::vector<Position*> & Chain::getPositions() {return positions;};
inline Position & Chain::getPosition(unsigned int _index) {return *(positions[_index]);}
//inline Position & Chain::getPositionByIndex(unsigned int _index) {std::cerr << "DEPRECATED Position & Chain::getPositionByIndex(unsigned int _index)" << std::endl; return *(positions[_index]);}
//inline Position & Chain::getPositionByIndexTmp(unsigned int _index) {std::cerr << "DEPRECATED Position & Chain::getPositionByIndex(unsigned int _index)" << std::endl; return *(positions[_index]);}
inline Position & Chain::getPosition(std::string _positionId) {
	if (positionExists(_positionId)) {
		return getLastFoundPosition();
	} else {
	        throw MslNotFoundException(MslTools::stringf("Chain::getPosition(): %s",_positionId.c_str()));
		// we should add try... catch support here
		//std::cerr << "ERROR 63812: position " << _positionId << " does not exist in chain at inline inline Position & Chain::getPosition(string _positionId)" << std::endl;
		//exit(63812);
	}
}
inline Position & Chain::getPosition(int _resNum, std::string _iCode) {return *(positionMap[_resNum][_iCode]);}
//inline Residue & Chain::operator()(std::string _resNumAndIcode) {int resNum=0; std::string iCode; MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode); return positionMap[resNum][iCode]->getCurrentIdentity();}
//inline Residue & Chain::operator()(int _resNum) {return positionMap[_resNum][""]->getCurrentIdentity();}
//inline Residue & Chain::getResidueByIndex(unsigned int _index) {return (positions[_index])->getCurrentIdentity();}
inline Residue & Chain::getIdentity(unsigned int _index) {return (positions[_index])->getCurrentIdentity();}
inline Residue & Chain::getIdentity(std::string _identityId) {
	if (identityExists(_identityId)) {
		return getLastFoundIdentity();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 63817: identity " << _identityId << " does not exist in chain at inline Residue & Chain::getIdentity(string _identityId)" << std::endl;
		exit(63817);
	}
}
inline Residue & Chain::getIdentity(int _resNum, std::string _iCode, std::string _identity) {
	if (identityExists(_resNum, _iCode, _identity)) {
		return getLastFoundIdentity();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 63822: identity " << _resNum << _iCode << "," << _identity << " does not exist in chain at inline Residue & Chain::getIdentity(int _resNum, string _iCode, string _identity)" << std::endl;
		exit(63822);
	}
}
inline Residue & Chain::getResidue(unsigned int _index) {return getIdentity(_index);}
inline Residue & Chain::getResidue(std::string _identityId) {return getIdentity(_identityId);}
inline Residue & Chain::getResidue(int _resNum, std::string _iCode, std::string _identity) {return getIdentity(_resNum, _iCode, _identity);}
inline AtomPointerVector & Chain::getAtomPointers() {return activeAtoms;}
inline AtomPointerVector & Chain::getAllAtomPointers() {return activeAndInactiveAtoms;}
inline unsigned int Chain::atomSize() {return activeAtoms.size();}
inline unsigned int Chain::allAtomSize() {return activeAndInactiveAtoms.size();}
inline Position & Chain::operator()(std::string _positionId) {return getPosition(_positionId);}
inline Atom & Chain::operator[](unsigned int _index) {return getAtom(_index);}
inline Atom & Chain::operator[](std::string _atomId) {return getAtom(_atomId);}
inline Position & Chain::operator()(unsigned int _index) {return getPosition(_index);}
inline Atom & Chain::getAtom(unsigned int _index) {return *(activeAtoms[_index]);}
//inline bool Chain::exists(int _resNum) {return positionMap.find(_resNum) != positionMap.end();}
//inline bool Chain::exists(int _resNum) {foundPosition = positionMap.find(_resNum); return foundPosition != positionMap.end();}
inline Position & Chain::getLastFoundPosition() {return *(foundPosition->second);}
inline Residue & Chain::getLastFoundIdentity() {return foundPosition->second->getLastFoundIdentity();}
inline Residue & Chain::getLastFoundResidue() {return getLastFoundIdentity();}
inline Atom & Chain::getLastFoundAtom() {return foundPosition->second->getLastFoundIdentity().getLastFoundAtom();}
inline void Chain::wipeAllCoordinates() {for (std::vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {(*k)->wipeAllCoordinates();}}
inline bool Chain::positionExists(std::string _positionId) {
	// this accepts either "37" or "37A", or even "A,37" (the chain id is ignored
	std::string chainid;
	int resnum;
	std::string icode;
	bool OK = MslTools::parsePositionId(_positionId, chainid, resnum, icode, 1);
	if (OK) {
		return positionExists(resnum, icode);
	}
	return false;
}
inline bool Chain::positionExists(int _resNum, std::string _icode) {
	std::map<int, std::map<std::string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find(_icode);
		return foundPosition != found->second.end();
	}
	return false;
}
inline bool Chain::identityExists(std::string _identityId) {
	// this accepts either "37,LEU" or "37A,LEU", or even "A,37,LEU" (the chain id is ignored
	std::string chainid;
	int resnum;
	std::string icode;
	std::string identity;
	bool OK = MslTools::parseIdentityId(_identityId, chainid, resnum, icode, identity, 1);
	if (!OK) {
		// was the residue identity not specified ("A,37")
		OK = MslTools::parsePositionId(_identityId, chainid, resnum, icode, 1);
		if (OK) {
			// "A,37" use default residue name
			identity = "";
		}
	}
	if (OK) {
		return identityExists(resnum, icode, identity);
	}
	return false;
}
inline bool Chain::identityExists(int _resNum, std::string _icode, std::string _identity) {
	std::map<int, std::map<std::string, Position*> >::iterator found=positionMap.find(_resNum);
	if (found != positionMap.end()) {
		foundPosition=found->second.find(_icode);
		if (foundPosition != found->second.end()) {
			return foundPosition->second->identityExists(_identity);
		}
	}
	return false;
}
inline bool Chain::atomExists(std::string _atomId) {
	// this accepts either "CA" or "ILE,CA", or even "A,37,ILE,CA" (the chain and resnum are ignored
	std::string chainid;
	int resnum;
	std::string icode;
	std::string atomName;
	bool OK = MslTools::parseAtomId(_atomId, chainid, resnum, icode, atomName, 1);
	if (OK) {
		return atomExists(resnum, icode, atomName);
	} else {
		std::string identity;
		OK = MslTools::parseAtomOfIdentityId(_atomId, chainid, resnum, icode, identity, atomName, 1);
		if (OK) {
			return atomExists(resnum, icode, identity, atomName);
		}
	}
	return false;
}
inline bool Chain::atomExists(int _resNum, std::string _icode, std::string _name) {
	if (positionExists(_resNum, _icode)) {
		return foundPosition->second->atomExists(_name);
	}
	return false;
}
inline bool Chain::atomExists(int _resNum, std::string _icode, std::string _identity, std::string _name) {
	if (identityExists(_resNum, _icode, _identity)) {
		return foundPosition->second->getLastFoundIdentity().atomExists(_name);
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 43849: atom " << _resNum << _icode << "," << _identity << "," << _name << " does not exist in inline bool Chain::atomExists(int _resNum, std::string _icode, std::string _identity, std::string _name)" << std::endl;
		exit(43853);
	}
}
inline Atom & Chain::getAtom(int _resnum, std::string _icode, std::string _name) {
	if (atomExists(_resnum, _icode, _name)) {
		return getLastFoundAtom();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 43853: atom " << _resnum << _icode << "," << _name << " does not exist in inline Atom & getAtom(int _resnum, string _icode, string _name)" << std::endl;
		exit(43853);
	}
}
inline Atom & Chain::getAtom(int _resnum, std::string _icode, std::string _identity, std::string _name) {
	if (atomExists(_resnum, _icode, _identity, _name)) {
		return getLastFoundAtom();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 43858: atom " << _resnum << _icode << "," << _identity << "," << _name << " does not exist in inline Atom & getAtom(int _resnum, string _icode, string _identity, string _name)" << std::endl;
		exit(43858);
	}
}
inline Atom & Chain::getAtom(std::string _atomId) {
	if (atomExists(_atomId)) {
		return getLastFoundAtom();
	} else {
		// we should add try... catch support here
		std::cerr << "ERROR 43863: atom " << _atomId << " does not exist in inline Atom & getAtom(string _atomId)" << std::endl;
		exit(43863);
	}
}
inline void Chain::saveCoor(std::string _coordName) {activeAndInactiveAtoms.saveCoor(_coordName);}
inline void Chain::saveAltCoor(std::string _coordName) {activeAndInactiveAtoms.saveAltCoor(_coordName);}
inline bool Chain::applySavedCoor(std::string _coordName) {return activeAndInactiveAtoms.applySavedCoor(_coordName);}
inline void Chain::clearSavedCoor(std::string _coordName) {activeAndInactiveAtoms.clearSavedCoor(_coordName);}

inline bool Chain::defineRotamerSamplingLevels(std::map<std::string,std::map<std::string,unsigned int> >& _levels) {
	bool success = true;
	for(std::vector<Position*>::iterator pos = positions.begin(); pos != positions.end(); pos++) {
		success = success && (*pos)->defineRotamerSamplingLevels(_levels);
	}
	return success;
}
inline bool Chain::isPositionNterminal(const Position * _pPos) const {
	if (positions.size() > 0 && _pPos == positions[0]) {
		return true;
	}
	return false;
}
inline bool Chain::isPositionCterminal(const Position * _pPos) const {
	if (positions.size() > 0 && _pPos == positions.back()) {
		return true;
	}
	return false;
}

}
#endif
