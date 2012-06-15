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

#ifndef ATOMANGLERELATIONSHIP_H
#define ATOMANGLERELATIONSHIP_H

#include "AtomGeometricRelationship.h"

namespace MSL { 
class AtomAngleRelationship : public AtomGeometricRelationship {
	public:
		AtomAngleRelationship();
		AtomAngleRelationship(const AtomPointerVector & _atoms);
		AtomAngleRelationship(const AtomAngleRelationship & _AAR);
		~AtomAngleRelationship();

	private:
		void calcValue();
		void checkSelected(std::string _selection1, std::string _selection2);

};

inline void AtomAngleRelationship::checkSelected(std::string _selection1, std::string _selection2) {
	// an angle is selected if the central atoms is present in both selections
	if ( atoms[1]->getSelectionFlag(_selection1) && atoms[1]->getSelectionFlag(_selection2) )  {
		selected = true;
	} else {
		selected = false;
	}
}

inline void AtomAngleRelationship::calcValue() {
	value = atoms[0]->angleRadians(*(atoms[1]),*(atoms[2]));
}


}

#endif
