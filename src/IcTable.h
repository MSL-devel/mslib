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
         *
         * NOTE:  The IcTable DOES NOT own the IcEntry pointers.
         *        It WILL NOT automatically delete the
         *        IcEntry pointers upon destruction.
         *        Responsibility for memory management is
         *        left to whomever created the IcEntry pointer.
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

		bool seed(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3);
		bool seed();


		void push_back(IcEntry * _ic);
                void deletePointers();

		void removeAtom(Atom * pAtom);

	private:
		void setup();
		void copy(const IcTable & _ic);
		

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
inline void IcTable::printTable() const {
	for (IcTable::const_iterator k=this->begin(); k!=this->end(); k++) {
		std::cout << **k << std::endl;
	}
}

} // end namespace MSL

#endif

