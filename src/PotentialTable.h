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

#ifndef POTENTIALTABLE_H
#define POTENTIALTABLE_H



/**
 * This class contains all the information about a given potential table.
 * A potential table is basically a map, with a generic string key and 
 * a double value result. 
 * 
 * 
 * 
 */

#include <string>
#include <map>
using namespace std;

class PotentialTable {
	public:
                	PotentialTable();


			void   addPotential(string _key, double _value);
	                double getPotential(string _key);


			void   setPotentialName(string _name);
			string getPotentialName();
	protected:

			string fileName;
			string potentialName;

			map<string,double> table;
};

/* INLINES */

inline PotentialTable::PotentialTable() { fileName = ""; potentialName = "";}
inline double PotentialTable::getPotential(string _key) { map<string,double>::iterator it = table.find(_key); if (it == table.end()) return 0.0; return it->second;}
inline void PotentialTable::addPotential(string _key, double _value) { table[_key] = _value; }
inline void   PotentialTable::setPotentialName(string _name) { potentialName = _name; }
inline string PotentialTable::getPotentialName() { return potentialName; }

#endif
