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

#ifndef THREEBODYINTERACTION_H
#define THREEBODYINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "Interaction.h"


namespace MSL { 
class ThreeBodyInteraction: public Interaction {

	/*******************************************************
	 *   Inherits from Interaction (a prototype object
	 *      AtomPointerVector pAtoms  (the atoms) 
	 *      std::vector<double> param  (the parameters
	 *******************************************************/

	public:

		virtual ~ThreeBodyInteraction();

		/* setting and getting the atoms */
		void setAtoms(std::vector<Atom*> _atoms);
		void setAtoms(Atom & _a1, Atom & _a2, Atom & _a3);

		double getAngle() const;
		double getAngleRadians() const;
		
		bool isSelected(std::string _selection1, std::string _selection2) const;
		bool isActive() const;

		virtual double getEnergy()=0;
		virtual double getEnergy(double _angle,std::vector<double> *_ad=NULL)=0;


		friend std::ostream & operator<<(std::ostream &_os, ThreeBodyInteraction & _term) {_os << _term.toString(); return _os;};
		virtual std::string toString() =0;

	protected:
		ThreeBodyInteraction();

};


inline void ThreeBodyInteraction::setAtoms(std::vector<Atom*> _atoms) { if (_atoms.size() != 3) {std::cerr << "ERROR 34812: invalid number of atoms in inline void ThreeBodyInteraction::setAtoms(std::vector<Atom*> _atoms)" << std::endl; exit(34812);} pAtoms = _atoms;}
inline void ThreeBodyInteraction::setAtoms(Atom & _a1, Atom & _a2, Atom & _a3) { pAtoms[0] = &_a1; pAtoms[1] = &_a2; pAtoms[2] = &_a3; }
inline double ThreeBodyInteraction::getAngle() const {return pAtoms[0]->angle(*pAtoms[1], *pAtoms[2]);}
inline double ThreeBodyInteraction::getAngleRadians() const {return pAtoms[0]->angleRadians(*pAtoms[1], *pAtoms[2]);}
inline bool ThreeBodyInteraction::isSelected(std::string _selection1, std::string _selection2) const {
	if (pAtoms[1]->getSelectionFlag(_selection1) && pAtoms[1]->getSelectionFlag(_selection2)) {
		return true;
	} else {
		return false;
	}
}
inline bool ThreeBodyInteraction::isActive() const {return pAtoms[0]->getActive() && pAtoms[1]->getActive() && pAtoms[2]->getActive();}


}

#endif

