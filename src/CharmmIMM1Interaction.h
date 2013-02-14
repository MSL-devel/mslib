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


#ifndef CHARMMIMM1INTERACTION_H
#define CHARMMIMM1INTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "CharmmEEF1Interaction.h"
#include "CharmmEnergy.h"


namespace MSL { 
	class CharmmIMM1Interaction: public TwoBodyInteraction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			CharmmIMM1Interaction();
			CharmmIMM1Interaction(Atom& _pA1, Atom& _pA2, std::vector<double>& _imm1W, std::vector<double>& _imm1C, double _halfThickness, double _exponent);

			// should implement an operator= as well 
			CharmmIMM1Interaction(const CharmmIMM1Interaction & _interaction);
			~CharmmIMM1Interaction();

			/* setting and getting the parameters */
			//double getMinD() const;
			//double getConstant() const;
			
			void setParams(std::vector<double>& _imm1W, std::vector<double>& _imm1C, double _halfThickness, double _exponent);
			double getEnergy();
			double getEnergy(double _dummy, std::vector<double> *_dd=NULL);// gradient not implemented for IMM1 yet
			double getEnergy(std::vector<double> *_dd); // gradient not implemented for IMM1 yet - energy computed by this function doesnot apply the switching function even if cutoffs are in place

			friend std::ostream & operator<<(std::ostream &_os, CharmmIMM1Interaction & _term) {_os << _term.toString(); return _os;};
			std::string toString() ;

			//unsigned int getType() const;
			std::string getName() const;
			void setUseNonBondCutoffs(bool _flag, double _ctonnb=0.0, double _ctofnb=0.0);

			std::vector<double> getEnergyGrad();
			std::pair<double,std::vector<double> > partialDerivative();

		
		private:
			void copy(const CharmmIMM1Interaction & _interaction);
			void setup(Atom* _pA1, Atom* _pA2);
			void deletePointers();

			//static const unsigned int type = 2;
			static const std::string typeName;
			CharmmEEF1Interaction* pEEFW;
			CharmmEEF1Interaction* pEEFC;
			

	};

	inline void CharmmIMM1Interaction::setParams(std::vector<double>& _imm1W, std::vector<double>& _imm1C,double _halfThickness, double _exponent)  {params[0] = _halfThickness; params[1] = _exponent; pEEFW->setParams(_imm1W); pEEFC->setParams(_imm1C);}
	inline double CharmmIMM1Interaction::getEnergy() {
		if(pEEFW) {
			double fz = CharmmEnergy::instance()->IMM1ZtransFunction(pAtoms[0]->getZ(),params[0],params[1]);
			return (fz * pEEFW->getEnergy() + (1-fz) * pEEFC->getEnergy());
		} else {
			std::cerr << "ERROR 12345: CharmmIMM1Interaction::getEnergy() is called without setting up the EEF1 interactions" << std::endl;
			return 0;
		}
	}
	inline double CharmmIMM1Interaction::getEnergy(std::vector<double> *_dd){
		if(_dd != NULL) {
			std::cerr << "WARNING 12234:  CharmmIMM1Interaction::getEnergy(double _distance, std::vector<double> *_dd) is not implemented to get the gradient" << std::endl;
			_dd->resize(pAtoms.size() * 2,0.0);
		}
		return getEnergy();
	}
	inline double CharmmIMM1Interaction::getEnergy(double _distance, std::vector<double> *_dd){
		if(_dd != NULL) {
			std::cerr << "WARNING 12234:  CharmmIMM1Interaction::getEnergy(double _distance, std::vector<double> *_dd) is not implemented to get the gradient" << std::endl;
			_dd->resize(pAtoms.size() * 2,0.0);
		}
		// if cutoffs are in force, this function is not the one to call - use CharmmIMM1Interaction::getEnergy(double _distance, double _groupDistance)
		return getEnergy();
	}
	inline std::string CharmmIMM1Interaction::toString() {
		char c [1000]; 
		sprintf(c, "%s %s %s %9.4f %9.4f %9.4f", typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1],getEnergy()); 
		return (std::string)c; 
	}
	//inline unsigned int CharmmIMM1Interaction::getType() const {return type;}
	inline std::string CharmmIMM1Interaction::getName() const {return typeName;}
	inline void CharmmIMM1Interaction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {
		if(pEEFW && pEEFC) {
			pEEFW->setUseNonBondCutoffs(_flag,_ctonnb,_ctofnb);
			pEEFC->setUseNonBondCutoffs(_flag,_ctonnb,_ctofnb);
		}
	}
	inline std::pair<double,std::vector<double> > CharmmIMM1Interaction::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
		return partials;
	}
	inline std::vector<double> CharmmIMM1Interaction::getEnergyGrad(){
		std::cerr << "WARNING 12234:  CharmmIMM1Interaction::getEnergyGrad() is not implemented to get the gradient" << std::endl;
		std::vector<double> tmp(pAtoms.size() * 3,0.0); 
		return tmp;
	}
}
#endif

