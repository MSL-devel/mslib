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

#ifndef ATOMGEOMETRICRELATIONSHIP_H
#define ATOMGEOMETRICRELATIONSHIP_H

#include <vector>
#include <string>

#include "AtomPointerVector.h"


/*******************************************************************
 *  This object is a template for 3 types of geometric relationships
 *  for atoms (distance, angle and dihedral angle)
 *
 *  The purpose of the object is to calculate the value of the 
 *  distance/angle for energy trms only once for a group of atoms and 
 *  store it if it is needed in the same cycle of energy calculation 
 *  without recomputing it (for example, the distance is first calculated 
 *  for the vdw term but not recomputed when the same std::pair of atoms
 *  is used by the electrostatic or solvation terms
 *
 *  The object knows if it is the same cycle (return stored) or new
 *  cycle (calculate) based on a "stamp", a value that is passed by
 *  the EnergySet that is unique to each cycle of energy calculation
 *******************************************************************/

namespace MSL { 
class AtomGeometricRelationship {
	public:
		virtual ~AtomGeometricRelationship();

		void setAtoms(AtomPointerVector _atoms);
		AtomPointerVector & getAtomPointers();
		double getValue(unsigned int _stamp=0);
		bool isSelected(std::string _selection1, std::string _selection2, unsigned int _stamp=0);

	protected:
		AtomGeometricRelationship();
	//	AtomGeometricRelationship(const AtomGeometricRelationship & _AGR);  // not needed
		virtual void calcValue()=0;
		virtual void checkSelected(std::string _selection1, std::string _selection2)=0;

		double value;
		bool selected;
		unsigned int stamp;
		
		AtomPointerVector atoms;

};

inline void AtomGeometricRelationship::setAtoms(AtomPointerVector _atoms) {atoms = _atoms;}
inline AtomPointerVector & AtomGeometricRelationship::getAtomPointers() {return atoms;}
inline double AtomGeometricRelationship::getValue(unsigned int _stamp) {
	if (_stamp == 0 || _stamp != stamp) {
		calcValue();
	}
	stamp = _stamp;
	return value;
}
inline bool AtomGeometricRelationship::isSelected(std::string _selection1, std::string _selection2, unsigned int _stamp) {
	if (_stamp == 0 || _stamp != stamp) {
		checkSelected(_selection1, _selection2);
	}
	stamp = _stamp;
	return selected;
}
//inline bool AtomGeometricRelationship::isSelected(std::string _selection1, std::string _selection2, unsigned int _stamp) {return false;}
//inline void AtomGeometricRelationship::calcValue() {}

}

#endif

