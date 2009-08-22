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

#ifndef CHARMMVDWINTERACTION_H
#define CHARMMVDWINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "CharmmEnergy.h"

using namespace std;

class CharmmVdwInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmVdwInteraction();
		CharmmVdwInteraction(Atom & _a1, Atom & _a2, double _rmin, double _Emin);

		// add an operator= 
		CharmmVdwInteraction(const CharmmVdwInteraction & _interaction);
		~CharmmVdwInteraction();

		/* setting and getting the parameters */
		void setParams(vector<double> _params);
		void setParams(double _rmin, double _Emin);
		double getRmin() const;
		double getEmin() const;
		
		double getEnergy();
		double getEnergy(double _distance);

		friend ostream & operator<<(ostream &_os, CharmmVdwInteraction & _term) {_os << _term.toString(); return _os;};
		string toString() const;

		//unsigned int getType() const;
		string getName() const;
		
	private:
		void setup(Atom * _a1, Atom * _a2, double _rmin, double _Emin);
		void copy(const CharmmVdwInteraction & _interaction);
		double distance;

		//static const unsigned int type = 0;
		static const string typeName;
		

};

inline void CharmmVdwInteraction::setParams(vector<double> _params) { if (_params.size() != 2) {cerr << "ERROR 46123: invalid number of parameters in inline void CharmmVdwInteraction::setParams(vector<double> _params)" << endl; exit(46123);} params = _params;}
inline void CharmmVdwInteraction::setParams(double _rmin, double _Emin) {params[0] = _rmin; params[1] = _Emin;}
inline double CharmmVdwInteraction::getRmin() const {return params[0];};
inline double CharmmVdwInteraction::getEmin() const {return params[1];};
inline double CharmmVdwInteraction::getEnergy() {
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}
inline double CharmmVdwInteraction::getEnergy(double _distance) {
	distance = _distance;
	energy = CharmmEnergy::instance()->LJ(_distance, params[0], params[1]);
	return energy;
}
inline string CharmmVdwInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM VDW %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], distance, energy); return (string)c; };
//inline unsigned int CharmmVdwInteraction::getType() const {return type;}
inline string CharmmVdwInteraction::getName() const {return typeName;}

#endif

