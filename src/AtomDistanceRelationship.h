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

#ifndef ATOMDISTANCERELATIONSHIP_H
#define ATOMDISTANCERELATIONSHIP_H

#include "AtomGeometricRelationship.h"

class AtomDistanceRelationship : public AtomGeometricRelationship {
	public:
		AtomDistanceRelationship();
		AtomDistanceRelationship(const AtomPointerVector & _atoms);
		AtomDistanceRelationship(const AtomDistanceRelationship & _ADR);
		~AtomDistanceRelationship();

	private:
		void calcValue();
		void checkSelected(string _selection1, string _selection2);

};

inline void AtomDistanceRelationship::checkSelected(string _selection1, string _selection2) {
	if ( (atoms[0]->getSelectionFlag(_selection1) && atoms[1]->getSelectionFlag(_selection2)) || (atoms[0]->getSelectionFlag(_selection2) && atoms[1]->getSelectionFlag(_selection1)) ) {
		selected = true;
	} else {
		selected = false;
	}
}

inline void AtomDistanceRelationship::calcValue() {
	value = atoms[0]->distance(*atoms[1]);
}


#endif
