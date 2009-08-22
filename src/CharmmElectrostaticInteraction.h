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

#ifndef CHARMMELECTROSTATICINTERACTION_H
#define CHARMMELECTROSTATICINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "CharmmEnergy.h"

using namespace std;

class CharmmElectrostaticInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmElectrostaticInteraction();
		CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, bool _is14=false);

		// add an operator= 
		CharmmElectrostaticInteraction(const CharmmElectrostaticInteraction & _interaction);
		~CharmmElectrostaticInteraction();

		/* setting and getting the parameters */
		void setParams(vector<double> _params);
		double getDielectricConstant() const;
		double getElec14factor() const;
		
		double getEnergy();
		double getEnergy(double _distance);

		friend ostream & operator<<(ostream &_os, CharmmElectrostaticInteraction & _term) {_os << _term.toString(); return _os;};
		string toString() const;

		//unsigned int getType() const;
		string getName() const;

		void update();
		
	private:
		void setup(Atom * _a1, Atom * _a2, bool _is14);
		void copy(const CharmmElectrostaticInteraction & _interaction);
		//static const unsigned int type = 1;
		static const string typeName;

		double distance;

		bool is14;
		double Kq_q1_q1_rescal_over_diel;
		bool useRiel;
		

};

inline void CharmmElectrostaticInteraction::setParams(vector<double> _params) { if (_params.size() != 0) {cerr << "ERROR 41822: invalid number of parameters in inline void CharmmElectrostaticInteraction::setParams(vector<double> _params)" << endl; exit(41822);} params = _params;}
inline double CharmmElectrostaticInteraction::getDielectricConstant() const {return params[0];};
inline double CharmmElectrostaticInteraction::getElec14factor() const {return params[1];};
inline double CharmmElectrostaticInteraction::getEnergy() {
	return getEnergy(getDistance());
}
inline double CharmmElectrostaticInteraction::getEnergy(double _distance) {
	distance = _distance;
	if (useRiel) {
		energy = CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel/_distance);
		//return CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel);
	} else {
		// R-dependent dielectric, divideant by distance
		energy = CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel);
		//return CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel/_distance);
	}
	return energy;
};
inline string CharmmElectrostaticInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM ELEC %s %s %+6.3f %+6.3f %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), pAtoms[0]->getCharge(), pAtoms[1]->getCharge(), params[0], params[1], distance, energy); return (string)c; };
//inline unsigned int CharmmElectrostaticInteraction::getType() const {return type;}
inline string CharmmElectrostaticInteraction::getName() const {return typeName;}

#endif

