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

#ifndef CHAIN_H
#define CHAIN_H

#include <map>
#include <string>
#include <vector>

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

		void addAtoms(const AtomPointerVector & _atoms);

		bool addIdentityToPosition(AtomPointerVector _atoms, std::string _name, unsigned int _resNum, std::string _iCode="");
		bool addIdentityToPosition(const Residue & _residue, unsigned int _resNum, std::string _iCode="");

		void setParentSystem(System * _system);
		System * getParentSystem() const;


		unsigned int size() const; // number of positions
		std::vector<Position*> & getPositions();
		//Position & getPosition(size_t _n);
		Position & getPositionByIndex(size_t _n);
		Position & getPosition(std::string _resNumAndIcode);
		Position & getPosition(int _resNum, std::string _iCode="");
		unsigned int getPositionIndex(const Position * _pPos) const;	

		Residue & operator()(std::string _resNumAndIcode); // returns the active residue at position _resNumAndIcode (enter as "75" or "75A")
		Residue & operator()(int _resNum); // same as above, but it takes and int and assumes the insertion code to be blank (i.e. 75)
		Residue & getResidueByIndex(size_t _n);
		Residue & getResidue(int _resNum, std::string _iCode="");
		AtomPointerVector & getAtoms();
		AtomPointerVector & getAllAtoms();
		unsigned int atomSize();
		unsigned int allAtomSize();
		Atom & operator[](size_t _n);
		Atom & getAtom(size_t _n);
		
		// check the existance of chains, residues, atoms
		bool exists(int _resNum); // residue, by int
		bool exists(std::string _resNumAndIcode); // residue by std::string (possibly with insertion code
		bool exists(int _resNum, std::string _name); // atom
		bool exists(std::string _resNumAndIcode, std::string _name);
		bool exists(int _resNum, std::string _name, std::string _identity); // atom w/ specified identity
		bool exists(std::string _resNumAndIcode, std::string _name, std::string _identity);

		Position & getLastFoundPosition();
		Residue & getLastFoundResidue();
		Atom & getLastFoundAtom();

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		/* UPDATES REQUESTED BY POSITIONS */
		void updatePositionMap(Position * _position);
		void updateIndexing();
		void updateAllAtomIndexing();
	//	void swapInActiveList(Position * _position, AtomPointerVector & _atoms);
	
		/* RENUMBER THE WHOLE CHAIN */
		void renumberChain(int _start);

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
inline unsigned int Chain::size() const {return positions.size();}
inline std::vector<Position*> & Chain::getPositions() {return positions;};
//inline Position & Chain::getPosition(size_t _n) {return *(positions[_n]);}
inline Position & Chain::getPositionByIndex(size_t _n) {return *(positions[_n]);}
inline Position & Chain::getPosition(std::string _resNumAndIcode) {int resNum=0; std::string iCode; MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode); return *(positionMap[resNum][iCode]);}
inline Position & Chain::getPosition(int _resNum, std::string _iCode) {return *(positionMap[_resNum][_iCode]);}
inline Residue & Chain::operator()(std::string _resNumAndIcode) {int resNum=0; std::string iCode; MslTools::splitIntAndString(_resNumAndIcode, resNum, iCode); return positionMap[resNum][iCode]->getCurrentIdentity();}
inline Residue & Chain::operator()(int _resNum) {return positionMap[_resNum][""]->getCurrentIdentity();}
inline Residue & Chain::getResidueByIndex(size_t _n) {return (positions[_n])->getCurrentIdentity();}
inline Residue & Chain::getResidue(int _resNum, std::string _iCode) {return positionMap[_resNum][_iCode]->getCurrentIdentity();}
inline AtomPointerVector & Chain::getAtoms() {return activeAtoms;}
inline AtomPointerVector & Chain::getAllAtoms() {return activeAndInactiveAtoms;}
inline unsigned int Chain::atomSize() {return activeAtoms.size();}
inline unsigned int Chain::allAtomSize() {return activeAndInactiveAtoms.size();}
inline Atom & Chain::operator[](size_t _n) {return *(activeAtoms[_n]);}
inline Atom & Chain::getAtom(size_t _n) {return *(activeAtoms[_n]);}
//inline bool Chain::exists(int _resNum) {return positionMap.find(_resNum) != positionMap.end();}
//inline bool Chain::exists(int _resNum) {foundPosition = positionMap.find(_resNum); return foundPosition != positionMap.end();}
inline Position & Chain::getLastFoundPosition() {return *(foundPosition->second);}
inline Residue & Chain::getLastFoundResidue() {return foundPosition->second->getCurrentIdentity();}
inline Atom & Chain::getLastFoundAtom() {return foundPosition->second->getLastFoundAtom();}
inline void Chain::wipeAllCoordinates() {for (std::vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {(*k)->wipeAllCoordinates();}}
inline void Chain::renumberChain(int _start) {
	for (unsigned int i=0; i<positions.size(); i++) {
		positions[i]->renumberNoUpdate(_start + i);
	}
	positionMap.clear();
	for (std::vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {
		positionMap[positions.back()->getResidueNumber()][positions.back()->getResidueIcode()] = *k;
	}
}

}

#endif
