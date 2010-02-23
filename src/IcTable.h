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

#ifndef ICTABLE_H
#define ICTABLE_H

#include <map>

#include "IcEntry.h"


namespace MSL { 
class IcTable : public std::vector<IcEntry*> {
	/****************************************************
	 *  This class builds atom a1 based on a2 a3 a4 and
	 *  the a1-a2 distance, a1-a2-a3 angle and a1-a2-a3-a4
	 *  dihedral
	 *
	 *  a1-a2-a3-a4
	 ****************************************************/
	public:
		IcTable();
		IcTable(const IcTable & _ic);
		~IcTable();

		/********************************************************
		 *  Set the internal coordinates from existing coordinates
		 ********************************************************/
		void fillFromCoor();
		//void buildAll(); // build all possible coordinates from the IC table
		void printTable() const;

		/********************************************************
		 *  Save and restore IC entries to buffers
		 ********************************************************/
		void saveToBuffer(std::string _name);
		bool restoreFromBuffer(std::string _name);
		void clearAllBuffers();
		void updateMap(Atom * _pOldAtom, Atom * _pNewAtom);

		bool editBond(Atom * _pAtom1, Atom * _pAtom2, double _newValue);
		bool editAngle(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, double _newValue);
		bool editDihedral(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _newValue);


		void push_back(IcEntry * _ic);

	private:
		void setup();
		void copy(const IcTable & _ic);
		void deletePointers();

		void mapValues(IcEntry * _ic);

		void replaceInMap(std::map<Atom*, std::map<Atom*, std::map<Atom*, std::map<Atom*, std::vector<double*> > > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom);
		void replaceInMap(std::map<Atom*, std::map<Atom*, std::map<Atom*, std::vector<double*> > > > & _map, Atom * _pOldAtom, Atom * _pNewAtom);
		void replaceInMap(std::map<Atom*, std::map<Atom*, std::vector<double*> > > & _map, Atom * _pOldAtom, Atom * _pNewAtom);
		void replaceInMap(std::map<Atom*, std::vector<double*> > & _map, Atom * _pOldAtom, Atom * _pNewAtom);

		std::map<Atom*, std::map<Atom*, std::vector<double*> > > bondMap;
		std::map<Atom*, std::map<Atom*, std::map<Atom*, std::vector<double*> > > > angleMap;
		std::map<Atom*, std::map<Atom*, std::map<Atom*, std::map<Atom*, std::vector<double*> > > > > dihedralMap;

};

inline void IcTable::push_back(IcEntry * _ic) {
	_ic->setParentIcTable(this);
	
	std::vector<IcEntry*>::push_back(_ic);
	mapValues(_ic);
}

}

#endif

