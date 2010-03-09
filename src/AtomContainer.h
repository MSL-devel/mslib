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


#ifndef ATOMCONTAINER_H
#define ATOMCONTAINER_H

#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBReader.h"
#include "PDBWriter.h"


namespace MSL { 
class Residue;

/**************************************************
 *   TO DO
 *
 *   Reorganize the exist and getAtom(string) and the internal
 *   atom map to use the AtomId and AtomOfIdentityId, with
 *   some support for accepting both (with an internal check)
 **************************************************/

class AtomContainer {
	public:
		AtomContainer();
		AtomContainer(const AtomPointerVector & _atoms);
		AtomContainer(const AtomContainer & _AC);
		~AtomContainer();
		void operator=(const AtomContainer & _AC); // assignment
		
		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		/* ADD ATOMS TO THE END */
		void addAtom(const Atom & _atom);
		void addAtom(std::string _atomId, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0));
		void addAtom(std::string _atomId, double _x, double _y, double _z);
		void addAtoms(const AtomPointerVector & _atoms);

		/* INSERT ATOMS AT A CERTAIN POSITION */
		void insertAtom(const Atom & _atom, unsigned int _skipPositions);
		void insertAtom(std::string _atomId, const CartesianPoint & _coor, unsigned int _skipPositions);
		void insertAtom(std::string _atomId, double _x, double _y, double _z, unsigned int _skipPositions);
		void insertAtoms(const AtomPointerVector & _atoms, unsigned int _skipPositions);

		/* REMOVE ATOMS */
		bool removeAtom(std::string _atomId);
		void removeAtom(unsigned int _n);
		void removeAllAtoms();

		
		unsigned int size() const;
		Atom & operator[](size_t _n);
		Atom & operator()(std::string _atomId); // use argument as ("A,7,CA");
		Atom & getAtom(size_t _n);
		Atom & getAtom(std::string _atomId);
		AtomPointerVector getResidue(std::string _positionId) const;
		AtomPointerVector & getAtoms();

		bool atomExists(std::string _atomId);
		Atom & getLastFoundAtom();
		void updateAtomMap(Atom * _atom);

		/* I/O */
		bool readPdb(std::string _filename); // add atoms or alt coor
		bool writePdb(std::string _filename);

		// print the atom container using the AtomPointerVector toString
		std::string toString() const;
		friend std::ostream & operator<<(std::ostream &_os, const AtomContainer & _atomContainer)  {_os << _atomContainer.toString(); return _os;};

	private:
		void setup();
		void copy(const AtomContainer & _AC);
		void deletePointers();		
		//std::string getMapKey(std::string _chainId, int _resNum, std::string _iCode, std::string _name);
		void reset();

		AtomPointerVector atoms;
		std::map<std::string, Atom*> atomMap;
		std::map<std::string, Atom*> atomMapWithIdentities;
		std::map<std::string, Atom*>::iterator found; // without Identity
		std::map<std::string, Atom*>::iterator found2; // with Identity
		
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		PDBReader * pdbReader;
		PDBWriter * pdbWriter;
		
};
// inlined functions
inline unsigned int AtomContainer::size() const {return atoms.size();}
inline Atom & AtomContainer::operator[](size_t _n) {return *(atoms[_n]);}
inline Atom & AtomContainer::getAtom(size_t _n) {return *atoms[_n];}
inline AtomPointerVector & AtomContainer::getAtoms() {return atoms;}
inline Atom & AtomContainer::getLastFoundAtom() { return *(found->second);}

inline bool AtomContainer::readPdb(std::string _filename) {reset(); if (!pdbReader->open(_filename) || !pdbReader->read()) return false; addAtoms(pdbReader->getAtoms()); return true;}
inline bool AtomContainer::writePdb(std::string _filename) {if (!pdbWriter->open(_filename)) return false; bool result = pdbWriter->write(atoms); pdbWriter->close();return result;}
inline std::string AtomContainer::toString() const {return atoms.toString();}

}

#endif
