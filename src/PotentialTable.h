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

#ifndef POTENTIALTABLE_H
#define POTENTIALTABLE_H



/**
 * This class contains all the information about a given potential table.
 * A potential table is basically a std::map, with a generic std::string key and 
 * a double value result. 
 * 
 * 
 * 
 */

#include <string>
#include <map>

namespace MSL { 
class PotentialTable {
	public:
                	PotentialTable();


			void   addPotential(std::string _key, double _value);
	                double getPotential(std::string _key);


			void   setPotentialName(std::string _name);
			std::string getPotentialName();
	protected:

			std::string fileName;
			std::string potentialName;

			std::map<std::string,double> table;
};

/* INLINES */

inline PotentialTable::PotentialTable() { fileName = ""; potentialName = "";}
inline double PotentialTable::getPotential(std::string _key) { std::map<std::string,double>::iterator it = table.find(_key); if (it == table.end()) return 0.0; return it->second;}
inline void PotentialTable::addPotential(std::string _key, double _value) { table[_key] = _value; }
inline void   PotentialTable::setPotentialName(std::string _name) { potentialName = _name; }
inline std::string PotentialTable::getPotentialName() { return potentialName; }

}

#endif
