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

#ifndef SPRINGCONSTRAINTINTERACTION_H
#define SPRINGCONSTRAINTINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "CharmmEnergy.h"
#include "CartesianPoint.h"


namespace MSL { 

class SpringConstraintInteraction: public OneBodyInteraction {

	/*******************************************************
	 *   Inherits from Interaction (a prototype object
	 *   for interaction )
	 *   Models the spring constraints for constraint minimization
	 *******************************************************/

	public:
		SpringConstraintInteraction();
		SpringConstraintInteraction(Atom & _a1, double _Kb, double _b0);

		// should implement an operator= as well 
		SpringConstraintInteraction(const SpringConstraintInteraction & _interaction);
		~SpringConstraintInteraction();

		/* setting and getting the parameters */
		void setParams(std::vector<double> _params);
		void setParams(double _Kb, double _b0);
		void setSpringConstant(double _Kb);
		double getMinD() const;
		double getConstant() const;
		
		double getEnergy();
		double getEnergy(std::vector<double> *_dd);
		double getEnergy(double _distance, std::vector<double> *_dd=NULL);

		std::vector<double> getEnergyGrad();
		std::vector<double> getEnergyGrad(Atom& a1, double Kb, double b0);

		friend std::ostream & operator<<(std::ostream &_os, SpringConstraintInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		//unsigned int getType() const;
		std::string getName() const;
		CartesianPoint& getReferenceCoor() ;
		bool reset();
		std::pair<double,std::vector<double> > partialDerivative();
		
	private:
		void setup(Atom * _a1, double _Kb, double _b0);
		void copy(const SpringConstraintInteraction & _interaction);

		//static const unsigned int type = 2;
		static const std::string typeName;
		CartesianPoint referenceCoor;
		

};

inline void SpringConstraintInteraction::setParams(std::vector<double> _params) { if (_params.size() != 2) {std::cerr << "ERROR 49123: invalid number of parameters in inline void SpringConstraintInteraction::setParams(std::vector<double> _params)" << std::endl; exit(49123);} params = _params;}
inline void SpringConstraintInteraction::setParams(double _Kb, double _b0) {params[0] = _Kb; params[1] = _b0;}
inline void SpringConstraintInteraction::setSpringConstant(double _Kb) {params[0] = _Kb;}
inline double SpringConstraintInteraction::getMinD() const {return params[1];};
inline double SpringConstraintInteraction::getConstant() const {return params[0];};
inline double SpringConstraintInteraction::getEnergy() { return getEnergy(referenceCoor.distance(pAtoms[0]->getCoor())); }
inline bool SpringConstraintInteraction::reset() { referenceCoor = pAtoms[0]->getCoor(); return true;}  // makes springEnergy 0.0 by moving referenceCoor to current atomCoor 
inline double SpringConstraintInteraction::getEnergy(std::vector<double> *_dd) { 
	if(_dd) {
		double distance = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),referenceCoor,_dd);
		return getEnergy(distance,_dd); 
	}
	return getEnergy();
}
inline double SpringConstraintInteraction::getEnergy(double _distance,std::vector<double> *_dd) { 
	return CharmmEnergy::instance()->spring(_distance, params[0], params[1],_dd); 
}
inline std::string SpringConstraintInteraction::toString() { 
	char c [1000]; 
	sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), referenceCoor.toString().c_str(), params[0], params[1], referenceCoor.distance(pAtoms[0]->getCoor()), getEnergy()); 
	return (std::string)c; 
}
//inline unsigned int SpringConstraintInteraction::getType() const {return type;}
inline std::string SpringConstraintInteraction::getName() const {return typeName;}
inline CartesianPoint& SpringConstraintInteraction::getReferenceCoor() { return referenceCoor;}
inline std::vector<double> SpringConstraintInteraction::getEnergyGrad(){
	return getEnergyGrad(*pAtoms[0],params[0],params[1]);
}

inline std::pair<double,std::vector<double> > SpringConstraintInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),referenceCoor,&(partials.second));
	return partials;
}
}

#endif

