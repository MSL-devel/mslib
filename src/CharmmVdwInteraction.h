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
		
		double getEnergy(); // wrapper function
		double getEnergy(double _distance); // used with no cutoffs
		double getEnergy(double _distance, double _groupDistance);// used with cutoffs

		friend ostream & operator<<(ostream &_os, CharmmVdwInteraction & _term) {_os << _term.toString(); return _os;};
		string toString() const;

		//unsigned int getType() const;
		string getName() const;

		// use cutoffs for non bonded interactions
		void setUseNonBondCutoffs(bool _flag, double _ctonnb=0.0, double _ctofnb=0.0);
		bool getUseNonBondCutoffs() const;
		double getNonBondCutoffOn() const;
		double getNonBondCutoffOff() const;

/*
		// use cutoffs for non bonded interactions
		void setUseNonBondCutoffs(bool _flag);
		bool getUseNonBondCutoffs() const;

		// start point for the cutoff
		void setNonBondCutoffOn(double _cutoff);
		double getNonBondCutoffOn() const;

		// end point for the cutoff
		void setNonBondCutoffOff(double _cutoff);
		double getNonBondCutoffOff() const;
*/
		
	private:
		void setup(Atom * _a1, Atom * _a2, double _rmin, double _Emin);
		void copy(const CharmmVdwInteraction & _interaction);
		double distance;

		//static const unsigned int type = 0;
		static const string typeName;
		
		bool useNonBondCutoffs;
		double nonBondCutoffOn;
		double nonBondCutoffOff;

};

inline void CharmmVdwInteraction::setParams(vector<double> _params) { if (_params.size() != 2) {cerr << "ERROR 46123: invalid number of parameters in inline void CharmmVdwInteraction::setParams(vector<double> _params)" << endl; exit(46123);} params = _params;}
inline void CharmmVdwInteraction::setParams(double _rmin, double _Emin) {params[0] = _rmin; params[1] = _Emin;}
inline double CharmmVdwInteraction::getRmin() const {return params[0];};
inline double CharmmVdwInteraction::getEmin() const {return params[1];};
inline double CharmmVdwInteraction::getEnergy() {
	if (useNonBondCutoffs) {
		// with cutoffs
		return getEnergy(pAtoms[0]->distance(*pAtoms[1]), pAtoms[0]->groupDistance(*pAtoms[1]));
	} else {
		// no cutoffs
		return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
	}
}
inline double CharmmVdwInteraction::getEnergy(double _distance) {
	// called if there are no cutoffs
	distance = _distance;
	energy = CharmmEnergy::instance()->LJ(_distance, params[0], params[1]);
	return energy;
}
inline double CharmmVdwInteraction::getEnergy(double _distance, double _groupDistance) {
	// called if there are cutoffs
	distance = _distance;
	double factor = 1.0;
	if (_groupDistance  > nonBondCutoffOff) {
		// out of cutofnb, return 0
		energy = 0.0;
		return energy;
	} else if (_groupDistance > nonBondCutoffOn) {
		// between cutofnb and cutonnb, calculate the switching factor based on the distance
		// between the geometric centers of the atom groups that the two atoms belong to
		factor = CharmmEnergy::instance()->switchingFunction(_groupDistance, nonBondCutoffOn, nonBondCutoffOff);
	}
	energy = CharmmEnergy::instance()->LJ(_distance, params[0], params[1]) * factor;
	return energy;
}
inline string CharmmVdwInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM VDW %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], distance, energy); return (string)c; };
//inline unsigned int CharmmVdwInteraction::getType() const {return type;}
inline string CharmmVdwInteraction::getName() const {return typeName;}
inline void CharmmVdwInteraction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {useNonBondCutoffs = _flag; nonBondCutoffOn = _ctonnb; nonBondCutoffOff = _ctofnb;}
inline bool CharmmVdwInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline double CharmmVdwInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline double CharmmVdwInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}
/*
inline void CharmmVdwInteraction::setUseNonBondCutoffs(bool _flag) {useNonBondCutoffs = _flag;}
inline bool CharmmVdwInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline void CharmmVdwInteraction::setNonBondCutoffOn(double _cutoff) {nonBondCutoffOn = _cutoff;}
inline double CharmmVdwInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline void CharmmVdwInteraction::setNonBondCutoffOff(double _cutoff) {nonBondCutoffOff = _cutoff;}
inline double CharmmVdwInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}
*/
#endif

