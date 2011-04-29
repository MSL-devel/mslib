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


namespace MSL { 
class CharmmElectrostaticInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmElectrostaticInteraction();
		//CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, bool _is14=false);
		CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, double _dielectricConstant=1.0, double _14rescaling=1.0, bool _useRdielectric=false);

		// add an operator= 
		CharmmElectrostaticInteraction(const CharmmElectrostaticInteraction & _interaction);
		~CharmmElectrostaticInteraction();

		/* setting and getting the parameters */
		void setParams(std::vector<double> _params);
		double getDielectricConstant() const;
		double getElec14factor() const;
		
		double getEnergy(); // wrapper function
		double getEnergy(double &_distance, std::vector<double> *_dd=NULL); // used with no cutoffs

		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, bool _is14=false);
		double getEnergy(double _distance, double _groupDistance);// used with cutoffs

		friend std::ostream & operator<<(std::ostream &_os, CharmmElectrostaticInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() const;

		//unsigned int getType() const;
		std::string getName() const;

		// use cutoffs for non bonded interactions
		void setUseNonBondCutoffs(bool _flag, double _ctonnb=0.0, double _ctofnb=0.0);
		bool getUseNonBondCutoffs() const;
		double getNonBondCutoffOn() const;
		double getNonBondCutoffOff() const;

		
	private:
		void setup(Atom * _pA1, Atom * _pA2, double _dielectricConstant, double _14rescaling, bool _useRdielectric);
		void copy(const CharmmElectrostaticInteraction & _interaction);
		//static const unsigned int type = 1;
		static const std::string typeName;
		void update();

		double distance;

		bool is14;
		double Kq_q1_q1_rescal_over_diel;
		bool useRiel;
		
		bool useNonBondCutoffs;
		double nonBondCutoffOn;
		double nonBondCutoffOff;

};

inline void CharmmElectrostaticInteraction::setParams(std::vector<double> _params) { if (_params.size() != 0) {std::cerr << "ERROR 41822: invalid number of parameters in inline void CharmmElectrostaticInteraction::setParams(std::vector<double> _params)" << std::endl; exit(41822);} params = _params;}
inline double CharmmElectrostaticInteraction::getDielectricConstant() const {return params[0];};
inline double CharmmElectrostaticInteraction::getElec14factor() const {return params[1];};
inline double CharmmElectrostaticInteraction::getEnergy() {
	if (useNonBondCutoffs) {
		// with cutoffs
	        distance = pAtoms[0]->distance(*pAtoms[1]), pAtoms[0]->groupDistance(*pAtoms[1]);
	} else {
		// no cutoffs
		distance = pAtoms[0]->distance(*pAtoms[1]);
	}
	return getEnergy(distance);
}
 inline double CharmmElectrostaticInteraction::getEnergy(double &_distance, std::vector<double> *_dd) {
	distance = _distance;
	if (useRiel) {
		energy = CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel/_distance);
	} else {
		// R-dependent dielectric, divideant by distance
		energy = CharmmEnergy::instance()->coulombEnerPrecomputed(_distance, Kq_q1_q1_rescal_over_diel);
	}

	if (_dd != NULL){
		CharmmEnergy::instance()->coulombEnerGrad(*_dd, _distance, Kq_q1_q1_rescal_over_diel,useRiel);
	} 

	return energy;
}
inline double CharmmElectrostaticInteraction::getEnergy(double _distance, double _groupDistance) {

	double factor = 0.0;
	if (useRiel) {
	  factor = Kq_q1_q1_rescal_over_diel/_distance;
	} else {
		// R-dependent dielectric, divideant by distance
	  factor = Kq_q1_q1_rescal_over_diel;
	}

	energy = CharmmEnergy::instance()->coulombEnerPrecomputedSwitched(_distance,factor,_groupDistance,nonBondCutoffOn,nonBondCutoffOff);

	return energy;
}
inline std::string CharmmElectrostaticInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM ELEC %s %s %+6.3f %+6.3f %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), pAtoms[0]->getCharge(), pAtoms[1]->getCharge(), params[0], params[1], distance, energy); return (std::string)c; };
//inline unsigned int CharmmElectrostaticInteraction::getType() const {return type;}
inline std::string CharmmElectrostaticInteraction::getName() const {return typeName;}
inline void CharmmElectrostaticInteraction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {useNonBondCutoffs = _flag; nonBondCutoffOn = _ctonnb; nonBondCutoffOff = _ctofnb;}
inline bool CharmmElectrostaticInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline double CharmmElectrostaticInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline double CharmmElectrostaticInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}
inline std::vector<double> CharmmElectrostaticInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],*pAtoms[1],is14);
}

}

#endif

