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

#ifndef ATOMDIHEDRALRELATIONSHIP_H
#define ATOMDIHEDRALRELATIONSHIP_H

#include "AtomGeometricRelationship.h"

class AtomDihedralRelationship : public AtomGeometricRelationship {
	public:
		AtomDihedralRelationship();
		AtomDihedralRelationship(const AtomPointerVector & _atoms);
		AtomDihedralRelationship(const AtomDihedralRelationship & _AAR);
		~AtomDihedralRelationship();
		void setDihedralSelectionType();
		void setImproperSelectionType();
		bool isDihedralSelectionType() const;

	private:
		void calcValue();
		void checkSelected(string _selection1, string _selection2);
		bool dihedralSelection; // if true select using the two central atoms (dihe), if false use the first atom (improper)

};

inline void AtomDihedralRelationship::checkSelected(string _selection1, string _selection2) {
	if (dihedralSelection) {
		// a dihe is selected if the two central atoms are in different selections
		if ( (atoms[1]->getSelectionFlag(_selection1) && atoms[2]->getSelectionFlag(_selection2)) || (atoms[1]->getSelectionFlag(_selection2) && atoms[2]->getSelectionFlag(_selection1)) ) {
			selected = true;
		} else {
			selected = false;
		}
	} else {
		// an improper is selected if the two central atoms are in different selections
		if ( atoms[0]->getSelectionFlag(_selection1) && atoms[0]->getSelectionFlag(_selection2)) {
			selected = true;
		} else {
			selected = false;
		}
	}
}

inline void AtomDihedralRelationship::calcValue() {
	value = atoms[0]->dihedralRadians(*(atoms[1]),*(atoms[2]),*(atoms[3]));
}
inline void AtomDihedralRelationship::setDihedralSelectionType() {dihedralSelection = true;}
inline void AtomDihedralRelationship::setImproperSelectionType() {dihedralSelection = false;}
inline bool AtomDihedralRelationship::isDihedralSelectionType() const {return dihedralSelection;}


#endif
