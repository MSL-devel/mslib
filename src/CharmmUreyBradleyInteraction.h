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

#ifndef CHARMMUREYBRADLEYINTERACTION_H
#define CHARMMUREYBRADLEYINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
class CharmmUreyBradleyInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		CharmmUreyBradleyInteraction();
		CharmmUreyBradleyInteraction(Atom & _a1, Atom & _a2, double _Kub, double _S0);

		// should implement an operator= as well 
		CharmmUreyBradleyInteraction(const CharmmUreyBradleyInteraction & _interaction);
		~CharmmUreyBradleyInteraction();

		/* setting and getting the parameters */
		void setParams(std::vector<double> _params);
		void setParams(double _Kub, double _S0);
		double getMinD() const;
		double getConstant() const;
		
		double getEnergy();
		double getEnergy(double _distance,std::vector<double> *_ad=NULL);
		std::vector<double> getEnergyGrad();

		friend std::ostream & operator<<(std::ostream &_os, CharmmUreyBradleyInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::string getName() const;
		std::pair<double,std::vector<double> > partialDerivative();
		
	private:
		void setup(Atom * _a1, Atom * _a2, double _Kub, double _S0);
		void copy(const CharmmUreyBradleyInteraction & _interaction);

		//static const unsigned int type = 4;
		static const std::string typeName;
		

};

inline void CharmmUreyBradleyInteraction::setParams(std::vector<double> _params) { if (_params.size() != 2) {std::cerr << "ERROR 49123: invalid number of parameters in inline void CharmmUreyBradleyInteraction::setParams(std::vector<double> _params)" << std::endl; exit(49123);} params = _params;}
inline void CharmmUreyBradleyInteraction::setParams(double _Kub, double _S0) {params[0] = _Kub; params[1] = _S0;}
inline double CharmmUreyBradleyInteraction::getMinD() const {return params[1];};
inline double CharmmUreyBradleyInteraction::getConstant() const {return params[0];};
inline double CharmmUreyBradleyInteraction::getEnergy() {
	 return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}
 inline double CharmmUreyBradleyInteraction::getEnergy(double _distance,std::vector<double> *_dd) {
	return CharmmEnergy::instance()->spring(_distance, params[0], params[1],_dd);
}
inline std::string CharmmUreyBradleyInteraction::toString() { 
	char c [1000]; 
	sprintf(c, "CHARMM UREY %s %s %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], pAtoms[0]->distance(*pAtoms[1]), getEnergy()); 
	return (std::string)c;
}
//inline unsigned int CharmmUreyBradleyInteraction::getType() const {return type;}
inline std::string CharmmUreyBradleyInteraction::getName() const {return typeName;}

inline std::vector<double> CharmmUreyBradleyInteraction::getEnergyGrad(){
	std::vector<double> result(6,0.0);
	return result;
}

inline std::pair<double,std::vector<double> > CharmmUreyBradleyInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
	return partials;
}
}

#endif

