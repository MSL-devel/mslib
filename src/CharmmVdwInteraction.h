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


namespace MSL { 
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
		void setParams(std::vector<double> _params);
		void setParams(double _rmin, double _Emin);
		double getRmin() const;
		double getEmin() const;
		
		double getEnergy(); // wrapper function
		double getEnergy(double _distance, std::vector<double> *_dd=NULL); // used with no cutoffs
		double getEnergy(double _distance, double _groupDistance, std::vector<double> *_dd=NULL);// used with cutoffs

		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, double rmin, double Emin);

		friend std::ostream & operator<<(std::ostream &_os, CharmmVdwInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::string getName() const;
		void setName(std::string _name);

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

		//static const unsigned int type = 0;
		std::string typeName;
		
		bool useNonBondCutoffs;
		double nonBondCutoffOn;
		double nonBondCutoffOff;

};

inline void CharmmVdwInteraction::setParams(std::vector<double> _params) { if (_params.size() != 2) {std::cerr << "ERROR 46123: invalid number of parameters in inline void CharmmVdwInteraction::setParams(std::vector<double> _params)" << std::endl; exit(46123);} params = _params;}
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
inline double CharmmVdwInteraction::getEnergy(double _distance, std::vector<double> *dd) {
	// called if there are no cutoffs
	return CharmmEnergy::instance()->LJ(_distance, params[0], params[1],dd);
}
inline double CharmmVdwInteraction::getEnergy(double _distance, double _groupDistance, std::vector<double> *dd) {

	return CharmmEnergy::instance()->LJSwitched(_distance, params[0], params[1],_groupDistance,nonBondCutoffOn,nonBondCutoffOff,dd);
}
inline std::string CharmmVdwInteraction::toString() { 
	char c [1000]; 
	sprintf(c, "CHARMM VDW %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1],pAtoms[0]->distance(*pAtoms[1]), getEnergy()); 
	return (std::string)c; 
}
//inline unsigned int CharmmVdwInteraction::getType() const {return type;}
inline std::string CharmmVdwInteraction::getName() const {return typeName;}
inline void CharmmVdwInteraction::setName(std::string _name) {typeName = _name;}
inline void CharmmVdwInteraction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {useNonBondCutoffs = _flag; nonBondCutoffOn = _ctonnb; nonBondCutoffOff = _ctofnb;}
inline bool CharmmVdwInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline double CharmmVdwInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline double CharmmVdwInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}

inline std::vector<double> CharmmVdwInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],*pAtoms[1],params[0],params[1]);
}

/*
inline void CharmmVdwInteraction::setUseNonBondCutoffs(bool _flag) {useNonBondCutoffs = _flag;}
inline bool CharmmVdwInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline void CharmmVdwInteraction::setNonBondCutoffOn(double _cutoff) {nonBondCutoffOn = _cutoff;}
inline double CharmmVdwInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline void CharmmVdwInteraction::setNonBondCutoffOff(double _cutoff) {nonBondCutoffOff = _cutoff;}
inline double CharmmVdwInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}
*/
}

#endif

