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

#ifndef ATOMCONTAINER_H
#define ATOMCONTAINER_H

#include "AtomVector.h"
#include "MslTools.h"
#include "PDBReader.h"
#include "PDBWriter.h"

using namespace std;

class Residue;

class AtomContainer {
	public:
		AtomContainer();
		AtomContainer(const AtomVector & _atoms);
		AtomContainer(const AtomContainer & _AC);
		~AtomContainer();
		void operator=(const AtomContainer & _AC); // assignment
		
		void setNameSpace(string _nameSpace);
		string getNameSpace() const;

		/* ADD ATOMS TO THE END */
		void addAtom(const Atom & _atom);
		void addAtom(string _name, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0));
		void addAtom(string _name, double _x, double _y, double _z);
		void addAtoms(const AtomVector & _atoms);

		/* INSERT ATOMS AT A CERTAIN POSITION */
		void insertAtom(const Atom & _atom, unsigned int _skipPositions);
		void insertAtom(string _name, const CartesianPoint & _coor, unsigned int _skipPositions);
		void insertAtom(string _name, double _x, double _y, double _z, unsigned int _skipPositions);
		void insertAtoms(const AtomVector & _atoms, unsigned int _skipPositions);

		/* REMOVE ATOMS */
		bool removeAtom(string _name);
		void removeAtom(unsigned int _n);
		void removeAllAtoms();

		
		unsigned int size() const;
		Atom & operator[](size_t _n);
		Atom & operator()(string _chain_resnum_name); // use argument as ("A 7 CA");
		Atom & getAtom(size_t _n);
		Atom & getAtom(string _chain_resnum_name);
		AtomVector & getAtoms();

		bool exists(string _chain_resnum_name);
		Atom & getLastFoundAtom();
		void updateAtomMap(Atom * _atom);

		/* I/O */
		bool readPdb(string _filename); // add atoms or alt coor
		bool writePdb(string _filename);

		// print the atom container using the AtomVector toString
		string toString() const;
		friend ostream & operator<<(ostream &_os, const AtomContainer & _atomContainer)  {_os << _atomContainer.toString(); return _os;};

	private:
		void setup();
		void copy(const AtomContainer & _AC);
		void deletePointers();		
		string getMapKey(string _chainId, int _resNum, string _iCode, string _name);
		void reset();

		AtomVector atoms;
		map<string, Atom*> atomMap;
		map<string, Atom*>::iterator found;
		
		string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		PDBReader * pdbReader;
		PDBWriter * pdbWriter;
		
};
// inlined functions
inline unsigned int AtomContainer::size() const {return atoms.size();}
inline Atom & AtomContainer::operator[](size_t _n) {return *(atoms[_n]);}
inline Atom & AtomContainer::getAtom(size_t _n) {return *atoms[_n];}
inline AtomVector & AtomContainer::getAtoms() {return atoms;}
inline Atom & AtomContainer::getLastFoundAtom() { return *(found->second);}

inline bool AtomContainer::readPdb(string _filename) {reset(); if (!pdbReader->open(_filename) || !pdbReader->read()) return false; addAtoms(pdbReader->getAtoms()); return true;}
inline bool AtomContainer::writePdb(string _filename) {if (!pdbWriter->open(_filename)) return false; bool result = pdbWriter->write(atoms); pdbWriter->close();return result;}
inline string AtomContainer::toString() const {return atoms.toString();}

#endif
