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


#ifndef CHARMMIMM1REFINTERACTION_H
#define CHARMMIMM1REFINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
	class CharmmIMM1RefInteraction: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			CharmmIMM1RefInteraction();
			CharmmIMM1RefInteraction(Atom & _a1, double _refW, double _refZ, double _halfThickness, double exponent);

			// should implement an operator= as well 
			CharmmIMM1RefInteraction(const CharmmIMM1RefInteraction & _interaction);
			~CharmmIMM1RefInteraction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params);
			void setParams(double _refW, double _refC, double _halfThickness, double _exponent);
			//double getMinD() const;
			//double getConstant() const;
			std::vector<double> getParams() const;
			
			double getEnergy();
			double getEnergy(double _param, std::vector<double> *_paramDerivatives=NULL);
			double getEnergy(std::vector<double> *_paramDerivatives);

			std::vector<double> getEnergyGrad(){ std::vector<double> tmp; return tmp;}

			friend std::ostream & operator<<(std::ostream &_os, CharmmIMM1RefInteraction & _term) {_os << _term.toString(); return _os;};
			std::string toString() ;

			//unsigned int getType() const;
			std::string getName() const;
			std::pair<double,std::vector<double> > partialDerivative();
			
		private:
			void setup(Atom * _a1, double _refW, double _refC, double _halfThickness, double _exponent);
			void copy(const CharmmIMM1RefInteraction & _interaction);

			//static const unsigned int type = 2;
			static const std::string typeName;
			double refW;
			double refC;
			double halfThickness;
			double exponent;

			
	};

	inline void CharmmIMM1RefInteraction::setParams(std::vector<double> _params) { if (_params.size() != 4) {std::cerr << "ERROR 31235: invalid number of parameters in inline void CharmmIMM1RefInteraction::setParams(std::vector<double> _params)" << std::endl; exit(31235);} params = _params;}
	inline void CharmmIMM1RefInteraction::setParams(double _refW, double _refC, double _halfThickness, double _exponent) {params[0] = _refW; params[1] = _refC; params[2] = _halfThickness; params[3] = _exponent;}
	inline std::vector<double> CharmmIMM1RefInteraction::getParams() const {return params;};
	inline double CharmmIMM1RefInteraction::getEnergy() {
		if(pAtoms[0]) {
			double fZ = CharmmEnergy::instance()->IMM1ZtransFunction(pAtoms[0]->getZ(),params[2],params[3]);
			return (fZ * params[0] + (1-fZ) * params[1]);
		} else {
			std::cerr << "ERROR 31234: atom pointer NULL in inline void CharmmIMM1RefInteraction::getEnergy()" << std::endl;
			exit(31234);
		}
	}
	inline double CharmmIMM1RefInteraction::getEnergy(double _param, std::vector<double> *_paramDerivatives) { 
		if(_paramDerivatives) {
			_paramDerivatives->resize(pAtoms.size() * 3, 0.0);
		}
		return getEnergy(); 
	}
	inline double CharmmIMM1RefInteraction::getEnergy(std::vector<double> *_paramDerivatives) { 
		if(_paramDerivatives) {
			_paramDerivatives->resize(pAtoms.size() * 3, 0.0);
		}
		return getEnergy(); 
	}
	inline std::string CharmmIMM1RefInteraction::toString() { char c [1000]; sprintf(c, "%s %s %9.4f %9.4f %9.4f %9.4f %9.4f", typeName.c_str(), pAtoms[0]->toString().c_str(), params[0],params[1],params[2],params[3],getEnergy()); return (std::string)c; };
	//inline unsigned int CharmmIMM1RefInteraction::getType() const {return type;}
	inline std::string CharmmIMM1RefInteraction::getName() const {return typeName;}
	inline std::pair<double,std::vector<double> > CharmmIMM1RefInteraction::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = 0;
		partials.second.resize(pAtoms.size() * 3, 0.0);
		return partials;
	}

}

#endif

