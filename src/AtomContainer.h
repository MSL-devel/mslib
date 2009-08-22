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
		void addAtoms(const AtomVector & _atoms);

		/* INSERT ATOMS AT A CERTAIN POSITION */
		void insertAtom(const Atom & _atom, unsigned int _skipPositions);
		void insertAtom(string _name, const CartesianPoint & _coor, unsigned int _skipPositions);
		void insertAtoms(const AtomVector & _atoms, unsigned int _skipPositions);

		/* REMOVE ATOMS */
		bool removeAtom(string _name);
		void removeAtom(unsigned int _n);
		void removeAllAtoms();

		
		unsigned int size() const;
		Atom & operator[](size_t _n);
		Atom & operator()(string _chain_resnum_name); // use argument as ("A, 7, CA");
		Atom & getAtom(size_t _n);
		Atom & getAtom(string _chain_resnum_name);
		AtomVector & getAtoms();

		bool exists(string _chain_resnum_name);
		Atom & getFoundAtom();
		void updateAtomMap(Atom * _atom);

	private:
		void setup();
		void copy(const AtomContainer & _AC);
		void deletePointers();		
		string getMapKey(string _chainId, int _resNum, string _iCode, string _name);

		AtomVector atoms;
		map<string, Atom*> atomMap;
		map<string, Atom*>::iterator found;
		
		string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd
		
};
// inlined functions
inline unsigned int AtomContainer::size() const {return atoms.size();}
inline Atom & AtomContainer::operator[](size_t _n) {return *(atoms[_n]);}
inline Atom & AtomContainer::getAtom(size_t _n) {return *atoms[_n];}
inline AtomVector & AtomContainer::getAtoms() {return atoms;}
inline Atom & AtomContainer::getFoundAtom() { return *(found->second);}


#endif
