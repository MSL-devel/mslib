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

#ifndef THREEBODYINTERACTION_H
#define THREEBODYINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "Interaction.h"

using namespace std;

class ThreeBodyInteraction: public Interaction {

	/*******************************************************
	 *   Inherits from Interaction (a prototype object
	 *      AtomVector pAtoms  (the atoms) 
	 *      vector<double> param  (the parameters
	 *******************************************************/

	public:

		~ThreeBodyInteraction();

		/* setting and getting the atoms */
		void setAtoms(vector<Atom*> _atoms);
		void setAtoms(Atom & _a1, Atom & _a2, Atom & _a3);

		double getAngle() const;
		double getAngleRadians() const;
		
		bool isSelected(string _selection1, string _selection2) const;
		bool isActive() const;

		virtual double getEnergy()=0;
		virtual double getEnergy(double _angle)=0;

		friend ostream & operator<<(ostream &_os, ThreeBodyInteraction & _term) {_os << _term.toString(); return _os;};
		virtual string toString() const=0;

	protected:
		ThreeBodyInteraction();

};


inline void ThreeBodyInteraction::setAtoms(vector<Atom*> _atoms) { if (_atoms.size() != 3) {cerr << "ERROR 34812: invalid number of atoms in inline void ThreeBodyInteraction::setAtoms(vector<Atom*> _atoms)" << endl; exit(34812);} pAtoms = _atoms;}
inline void ThreeBodyInteraction::setAtoms(Atom & _a1, Atom & _a2, Atom & _a3) { pAtoms[0] = &_a1; pAtoms[1] = &_a2; pAtoms[2] = &_a3; }
inline double ThreeBodyInteraction::getAngle() const {return pAtoms[0]->angle(*pAtoms[1], *pAtoms[2]);}
inline double ThreeBodyInteraction::getAngleRadians() const {return pAtoms[0]->angleRadians(*pAtoms[1], *pAtoms[2]);}
inline bool ThreeBodyInteraction::isSelected(string _selection1, string _selection2) const {
	if (pAtoms[1]->getSelectionFlag(_selection1) && pAtoms[1]->getSelectionFlag(_selection2)) {
		return true;
	} else {
		return false;
	}
}
inline bool ThreeBodyInteraction::isActive() const {return pAtoms[0]->getActive() && pAtoms[1]->getActive() && pAtoms[2]->getActive();}


#endif

