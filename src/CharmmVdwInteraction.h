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
		double getEnergy(std::vector<double> *_dd); // used by minimizer - does not apply switching function even if cutoffs are in place 
		double getEnergy(double _distance, std::vector<double> *_dd=NULL); // used with no cutoffs
		double getEnergy(double _distance, double _groupDistance);// used with cutoffs

		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, double rmin, double Emin);

		friend std::ostream & operator<<(std::ostream &_os, CharmmVdwInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::string getName() const;

		// use cutoffs for non bonded interactions
		void setUseNonBondCutoffs(bool _flag, double _ctonnb=0.0, double _ctofnb=0.0);
		bool getUseNonBondCutoffs() const;
		double getNonBondCutoffOn() const;
		double getNonBondCutoffOff() const;

		std::pair<double,std::vector<double> > partialDerivative();
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
		static const std::string typeName;
		
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
inline double CharmmVdwInteraction::getEnergy(std::vector<double> *_dd) {
	if(_dd) {
		// get the gradient
		double distance = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),_dd);
		/*
		if(useNonBondCutoffs) {
			std::cerr << "WARNING 56783: Gradient not implemented for CharmmVdwInteraction with cutoffs" << std::endl;
			double groupDistance = pAtoms[0]->groupDistance(*pAtoms[1]);
			return getEnergy(distance,groupDistance,_dd); 
			return getEnergy(distance,_dd);
		} else {
			return getEnergy(distance,_dd);
		}
		*/
		return getEnergy(distance,_dd);
	}
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}
inline double CharmmVdwInteraction::getEnergy(double _distance, std::vector<double> *_dd) {
	// called if there are no cutoffs
	return CharmmEnergy::instance()->LJ(_distance, params[0], params[1],_dd);
}
inline double CharmmVdwInteraction::getEnergy(double _distance, double _groupDistance) {

	return CharmmEnergy::instance()->LJSwitched(_distance, params[0], params[1],_groupDistance,nonBondCutoffOn,nonBondCutoffOff);
}
inline std::string CharmmVdwInteraction::toString() { 
	char c [1000]; 
	sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1],pAtoms[0]->distance(*pAtoms[1]), getEnergy()); 
	return (std::string)c; 
}
//inline unsigned int CharmmVdwInteraction::getType() const {return type;}
inline std::string CharmmVdwInteraction::getName() const {return typeName;}
inline void CharmmVdwInteraction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {useNonBondCutoffs = _flag; nonBondCutoffOn = _ctonnb; nonBondCutoffOff = _ctofnb;}
inline bool CharmmVdwInteraction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
inline double CharmmVdwInteraction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
inline double CharmmVdwInteraction::getNonBondCutoffOff() const {return nonBondCutoffOff;}

inline std::vector<double> CharmmVdwInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],*pAtoms[1],params[0],params[1]);
}

inline std::pair<double,std::vector<double> > CharmmVdwInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
	return partials;
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

