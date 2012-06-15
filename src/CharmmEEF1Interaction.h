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


#ifndef CHARMMEEF1INTERACTION_H
#define CHARMMEEF1INTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
	class CharmmEEF1Interaction: public TwoBodyInteraction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			CharmmEEF1Interaction();
			CharmmEEF1Interaction(Atom & _a1, Atom & _a2, double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j);

			// should implement an operator= as well 
			CharmmEEF1Interaction(const CharmmEEF1Interaction & _interaction);
			~CharmmEEF1Interaction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params);
			void setParams(double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j);
			//double getMinD() const;
			//double getConstant() const;
			std::vector<double> getParams() const;
			
			double getEnergy();
			double getEnergy(double _distance, std::vector<double> *_dd=NULL);// gradient not implemented for EEF1 yet
			double getEnergy(std::vector<double> *_dd); // gradient not implemented for EEF1 yet - energy computed by this function doesnot apply the switching function even if cutoffs are in place
			std::vector<double> getEnergyGrad();
			double getEnergy(double _distance, double _groupDistance);

			friend std::ostream & operator<<(std::ostream &_os, CharmmEEF1Interaction & _term) {_os << _term.toString(); return _os;};
			std::string toString() ;

			//unsigned int getType() const;
			std::string getName() const;

			// use cutoffs for non bonded interactions
			void setUseNonBondCutoffs(bool _flag, double _ctonnb=0.0, double _ctofnb=0.0);
			bool getUseNonBondCutoffs() const;
			double getNonBondCutoffOn() const;
			double getNonBondCutoffOff() const;
			std::pair<double,std::vector<double> > partialDerivative();
		
		private:
			void setup(Atom * _a1, Atom * _a2, double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j);
			void copy(const CharmmEEF1Interaction & _interaction);

			//static const unsigned int type = 2;
			static const std::string typeName;
			
			bool useNonBondCutoffs;
			double nonBondCutoffOn;
			double nonBondCutoffOff;

	};

	inline void CharmmEEF1Interaction::setParams(std::vector<double> _params) { if (_params.size() != 8) {std::cerr << "ERROR 91235: invalid number of parameters in inline void CharmmEEF1Interaction::setParams(std::vector<double> _params)" << std::endl; exit(91235);} params = _params;}
	inline void CharmmEEF1Interaction::setParams(double _V_i, double _Gfree_i, double _Sigw_i, double _rmin_i, double _V_j, double _Gfree_j, double _Sigw_j, double _rmin_j) {params[0] = _V_i; params[1] = _Gfree_i; params[2] = _Sigw_i; params[3] = _rmin_i; params[4] = _V_j; params[5] = _Gfree_j; params[6] = _Sigw_j; params[7] = _rmin_j;}
	inline std::vector<double> CharmmEEF1Interaction::getParams() const {return params;};
	inline double CharmmEEF1Interaction::getEnergy() {
		double distance = 0.0;
		if (useNonBondCutoffs) {
			// with cutoffs
		        return getEnergy(pAtoms[0]->distance(*pAtoms[1]), pAtoms[0]->groupDistance(*pAtoms[1]));
		} else {
			// no cutoffs
			distance = pAtoms[0]->distance(*pAtoms[1]);
		}
		return getEnergy(distance);
	}
	inline double CharmmEEF1Interaction::getEnergy(std::vector<double> *_dd){
		if(_dd != NULL) {
			std::cerr << "WARNING 12234:  CharmmEEF1Interaction::getEnergy(double _distance, std::vector<double> *_dd) is not implemented to get the gradient" << std::endl;
			_dd->resize(pAtoms.size() * 2,0.0);
		}
		return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
	}
	inline double CharmmEEF1Interaction::getEnergy(double _distance, std::vector<double> *_dd){
		if(_dd != NULL) {
			std::cerr << "WARNING 12234:  CharmmEEF1Interaction::getEnergy(double _distance, std::vector<double> *_dd) is not implemented to get the gradient" << std::endl;
			_dd->resize(pAtoms.size() * 2,0.0);
		}
		// if cutoffs are in force, this function is not the one to call - use CharmmEEF1Interaction::getEnergy(double _distance, double _groupDistance)
		return CharmmEnergy::instance()->EEF1Ener(_distance, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7]); 
	}
	inline double CharmmEEF1Interaction::getEnergy(double _distance, double _groupDistance) {
		// called if there are cutoffs
		double factor = 1.0;
		if (_groupDistance  > nonBondCutoffOff) {
			// out of cutofnb, return 0
			return 0.0;
		} else if (_groupDistance > nonBondCutoffOn) {
			// between cutofnb and cutonnb, calculate the switching factor based on the distance
			// between the geometric centers of the atom groups that the two atoms belong to
			factor = CharmmEnergy::instance()->switchingFunction(_groupDistance, nonBondCutoffOn, nonBondCutoffOff);
		}
		return CharmmEnergy::instance()->EEF1Ener(_distance, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7]) * factor;
	}
	inline std::vector<double> CharmmEEF1Interaction::getEnergyGrad(){
		std::cerr << "WARNING 12234:  CharmmEEF1Interaction::getEnergyGrad() is not implemented to get the gradient" << std::endl;
		std::vector<double> tmp(pAtoms.size() * 3,0.0); 
		return tmp;
	}
	inline std::string CharmmEEF1Interaction::toString() {
		char c [1000]; 
		sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], pAtoms[0]->distance(*pAtoms[1]), getEnergy()); 
		return (std::string)c; 
	}
	//inline unsigned int CharmmEEF1Interaction::getType() const {return type;}
	inline std::string CharmmEEF1Interaction::getName() const {return typeName;}
	inline void CharmmEEF1Interaction::setUseNonBondCutoffs(bool _flag, double _ctonnb, double _ctofnb) {useNonBondCutoffs = _flag; nonBondCutoffOn = _ctonnb; nonBondCutoffOff = _ctofnb;}
	inline bool CharmmEEF1Interaction::getUseNonBondCutoffs() const {return useNonBondCutoffs;}
	inline double CharmmEEF1Interaction::getNonBondCutoffOn() const {return nonBondCutoffOn;}
	inline double CharmmEEF1Interaction::getNonBondCutoffOff() const {return nonBondCutoffOff;}

	inline std::pair<double,std::vector<double> > CharmmEEF1Interaction::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = CartesianGeometry::distanceDerivative(pAtoms[0]->getCoor(),pAtoms[1]->getCoor(),&(partials.second));
		return partials;
	}
}

#endif

