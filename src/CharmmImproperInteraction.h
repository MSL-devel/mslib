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

#ifndef CHARMMIMPROPERINTERACTION_H
#define CHARMMIMPROPERINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "FourBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
class CharmmImproperInteraction: public FourBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmImproperInteraction();
		CharmmImproperInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, double _Kpsi, double _Psi0Radians); //Psi0 must be in Radians

		// should implement an operator= as well 
		CharmmImproperInteraction(const CharmmImproperInteraction & _interaction);
		~CharmmImproperInteraction();

		/* setting and getting the parameters */
		void setParams(std::vector<double> _params);
		void setParams(double _Kpsi, double _Psi0Radians);
		double getMinAngle() const;
		double getConstant() const;
		
		double getEnergy();
		double getEnergy(double &_angleDegrees,std::vector<double> *_ad=NULL);
		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, Atom& a2, Atom& a3, Atom& a4, double Kpsi, double Psi0Radians);


		friend std::ostream & operator<<(std::ostream &_os, CharmmImproperInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() const;

		//unsigned int getType() const;
		std::string getName() const;
		
		bool isSelected(std::string _selection1, std::string _selection2) const;

	private:
		void setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, double _Kpsi, double _Psi0Radians);
		void copy(const CharmmImproperInteraction & _interaction);
		double angle;

		//static const unsigned int type = 6;
		static const std::string typeName;
		

};

inline void CharmmImproperInteraction::setParams(std::vector<double> _params) { if (_params.size() != 2) {std::cerr << "ERROR 58128: invalid number of parameters in inline void CharmmImproperInteraction::setParams(std::vector<double> _params)" << std::endl; exit(58128);} params = _params;}
inline void CharmmImproperInteraction::setParams(double _Kpsi, double _Psi0Radians) {params[0] = _Kpsi; params[1] = _Psi0Radians;}
inline double CharmmImproperInteraction::getMinAngle() const {return params[1];};
inline double CharmmImproperInteraction::getConstant() const {return params[0];};
inline double CharmmImproperInteraction::getEnergy() {
	angle = pAtoms[0]->dihedralRadians(*pAtoms[1], *pAtoms[2], *pAtoms[3]);
	energy = CharmmEnergy::instance()->spring(angle, params[0], params[1]);
	return energy;
}
inline double CharmmImproperInteraction::getEnergy(double &_angleDegrees, std::vector<double> *_ad) {
	angle = _angleDegrees * M_PI / 180.0;
	energy = CharmmEnergy::instance()->spring(angle, params[0], params[1],_ad);
	return energy;
}
inline std::string CharmmImproperInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM IMPR %s %s %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), pAtoms[2]->toString().c_str(), pAtoms[3]->toString().c_str(), params[0], params[1], angle, energy); return (std::string)c; };
//inline unsigned int CharmmImproperInteraction::getType() const {return type;}
inline std::string CharmmImproperInteraction::getName() const {return typeName;}
inline bool CharmmImproperInteraction::isSelected(std::string _selection1, std::string _selection2) const {
	if (pAtoms[0]->getSelectionFlag(_selection1) && pAtoms[0]->getSelectionFlag(_selection2)) {
		return true;
	} else {
		return false;
	}
}

inline std::vector<double> CharmmImproperInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],*pAtoms[1],*pAtoms[2], *pAtoms[3],params[0],params[1]);
}

}

#endif

