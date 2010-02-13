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

#ifndef ATOMSELECTION_H
#define ATOMSELECTION_H

// STL includes
#include <sstream>
//#include <tr1/unordered_map>

// MSL includes
#include "LogicalParser.h"
#include "Hash.h"
#include "AtomPointerVector.h"



class AtomSelection {

	public:
		AtomSelection();
		AtomSelection(AtomPointerVector &_data);
		~AtomSelection();


		AtomPointerVector& select(string _selectString,bool _selectAllAtoms=false);
		
		
		void clearStoredSelection(string _name) {}
 		inline void clearStoredSelections() { storedSelections.clear(); }


		AtomPointerVector& getSelection(string _selectName);
		bool selectionExists(string _selectName);

		bool getDebugFlag();
		void setDebugFlag(bool _flag);

		inline unsigned int size(string _name);

	private:
		LogicalParser lp;
		AtomPointerVector *data;
		Hash<string,AtomPointerVector>::Table storedSelections;

		bool debug;
};

// INLINES
inline bool AtomSelection::getDebugFlag(){ return debug; }
inline void AtomSelection::setDebugFlag(bool _flag){ debug = _flag;}
inline unsigned int AtomSelection::size(string _name) {
	string name = MslTools::toUpper(_name);
	if (storedSelections.find(name) != storedSelections.end()) { 
		return storedSelections[name].size();
	} else { 
		return 0;
	}
}

#endif

