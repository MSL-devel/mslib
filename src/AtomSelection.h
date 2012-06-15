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


#ifndef ATOMSELECTION_H
#define ATOMSELECTION_H

// STL includes
#include <sstream>
//#include <tr1/unordered_map>

// MSL includes
#ifndef __TESTING__
#include "LogicalCondition.h"
#else
#include "LogicalParser.h"
#endif

#include "Hash.h"
#include "AtomPointerVector.h"

namespace MSL { 
class AtomSelection {

	public:
		AtomSelection();
		AtomSelection(AtomPointerVector &_data);
		~AtomSelection();


		AtomPointerVector& select(std::string _selectString,bool _selectAllAtoms=false);
		AtomPointerVector& inverseSelect(std::string _selectString,bool _selectAllAtoms=false);
		void clearStoredSelection(std::string _name);
 		void clearStoredSelections();


		AtomPointerVector& getSelection(std::string _selectName);
		bool selectionExists(std::string _selectName);
		unsigned int selectionSize(std::string _selectName);

		bool getDebugFlag();
		void setDebugFlag(bool _flag);

		inline unsigned int size(std::string _name);

	private:
#ifndef __TESTING__
		AtomPointerVector& logicalSelect(std::string _selectString, std::string _name, AtomPointerVector & _atoms, bool _selectAllAtoms);
		LogicalCondition cond;
#else
		LogicalParser lp;
#endif
		AtomPointerVector *data;
		Hash<std::string,AtomPointerVector>::Table storedSelections;

		bool debug;
};

// INLINES
inline bool AtomSelection::getDebugFlag(){ return debug; }
inline void AtomSelection::setDebugFlag(bool _flag){ debug = _flag;}
inline unsigned int AtomSelection::size(std::string _name) {
	std::string name = MslTools::toUpper(_name);
	if (storedSelections.find(name) != storedSelections.end()) { 
		return storedSelections[name].size();
	} else { 
		return 0;
	}
}


inline void AtomSelection::clearStoredSelection(std::string _name) {
	_name = MslTools::toUpper(_name);
	Hash<std::string,AtomPointerVector>::Table::iterator it = storedSelections.find(_name);
	if (it != storedSelections.end()){
		storedSelections.erase(it);
	}
	AtomPointerVector::iterator avIt;
	for (avIt = data->begin();avIt != data->end();avIt++){
		(*avIt)->clearFlag(_name);
	}
}
inline void AtomSelection::clearStoredSelections() { 
	storedSelections.clear(); 
	AtomPointerVector::iterator avIt;
	for (avIt = data->begin();avIt != data->end();avIt++){
		(*avIt)->clearAllFlags();
	}
	storedSelections["_EMPTY_"]; // create an empty selection
}
}

#endif

