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


#ifndef SCWRL4HBONDINTERACTION_H
#define SCWRL4HBONDINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "CharmmEnergy.h"
#include "Atom.h"


namespace MSL { 
	class Scwrl4HBondInteraction: public Interaction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			Scwrl4HBondInteraction();
			// all angles should be in radians
			Scwrl4HBondInteraction(Atom & _d1, Atom & _d2,Atom & _a1, Atom & _a2,Atom& _a3,double _dist, double _ang,double _e1_dihe,double _e2_dihe, double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max, double _scalingFactor = 1.0);

			
			// should implement an operator= as well 
			Scwrl4HBondInteraction(const Scwrl4HBondInteraction & _interaction);
			~Scwrl4HBondInteraction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params); // make sure the angles are in radians
			// all angles should be in radians
			void setParams(double _dist, double _ang, double _e1_dihe, double _e2_dihe,double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max);
			std::vector<double> getParams() const;
			
			double getEnergy();
			double getEnergy(std::vector<double> *paramDerivatives);
			double getEnergy(double _param, std::vector<double> *paramDerivatives=NULL);
			std::vector<double> getEnergyGrad();


			std::string toString() ;

			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, Scwrl4HBondInteraction & _term) {_os << _term.toString(); return _os;};
			bool isSelected (std::string _sele1, std::string _sele2) const;
			bool isActive () const;
			double getW() ;
			void setScalingFactor(double _scalingFactor);
			double getScalingFactor() const;
			void printParameters();
			std::pair<double,std::vector<double> > partialDerivative();

					
		private:
			void setup(Atom * _d1, Atom * _d2,Atom * _a1, Atom * _a2,Atom * _a3,double _d, double _e0,double _e1,double _e2, double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max, double _scalingFactor);
			void copy(const Scwrl4HBondInteraction & _interaction);

			//static const unsigned int type = 2;
			double scalingFactor;
			static const std::string typeName;
	};

	inline void Scwrl4HBondInteraction::setParams(std::vector<double> _params) { 
		if (_params.size() != 9) {
			std::cerr << "ERROR 91235: invalid number of parameters in inline void Scwrl4HBondInteraction::setParams(std::vector<double> _params)" << std::endl;
			 exit(91235);
		} 
		params = _params;
		// precompute cos_alpha_max and cos_beta_max
		params[7] = cos(_params[7]);
		params[8] = cos(_params[8]);
	}
	inline void Scwrl4HBondInteraction::setParams(double _dist, double _ang, double _e1_dihe, double _e2_dihe, double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max) {
		params[0] = _dist; 
		params[1] = _ang; 
		params[2] = _e1_dihe; 
		params[3] = _e2_dihe; 
		params[4] = _d0; 
		params[5] = _sig_d;
		params[6] = _B;
		params[7] = cos(_alpha_max);
		params[8] = cos(_beta_max );
	 }
	inline std::vector<double> Scwrl4HBondInteraction::getParams() const {return params;};
	inline bool Scwrl4HBondInteraction::isSelected(std::string _sele1, std::string _sele2) const {
		if((pAtoms[0]->getSelectionFlag(_sele1) && pAtoms[2]->getSelectionFlag(_sele2)) || (pAtoms[2]->getSelectionFlag(_sele1) && pAtoms[0]->getSelectionFlag(_sele2))) {
			return true;
		} else {
			return false;
		}
	}
	inline bool Scwrl4HBondInteraction::isActive() const {
		return pAtoms[0]->getActive() && pAtoms[1]->getActive() && pAtoms[2]->getActive() && pAtoms[3]->getActive() && pAtoms[4]->getActive();
	}
	inline std::string Scwrl4HBondInteraction::toString() { char c [1000]; sprintf(c, "%s %s %s %s %s %s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(),pAtoms[1]->toString().c_str(),pAtoms[2]->toString().c_str(),pAtoms[3]->toString().c_str(),pAtoms[4]->toString().c_str(), params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], getEnergy()); return (std::string)c; };
	//inline unsigned int Scwrl4HBondInteraction::getType() const {return type;}
	inline std::string Scwrl4HBondInteraction::getName() const {return typeName;}
	inline double Scwrl4HBondInteraction::getScalingFactor() const {return scalingFactor;}
	inline void Scwrl4HBondInteraction::setScalingFactor(double _scalingFactor) {scalingFactor = _scalingFactor;}

	inline std::pair<double,std::vector<double> > Scwrl4HBondInteraction::partialDerivative() {
		std::pair<double, std::vector<double> > partials;
		partials.first = 0.0;
		getEnergy(&(partials.second));
		return partials;
	}
	inline double Scwrl4HBondInteraction::getEnergy(std::vector<double> *_dd) {
		if(_dd) {
			std::cerr << "Scwrl4HBondInteraction::partialDerivative is not implemented" << std::endl;
			_dd->resize(pAtoms.size(),0.0);
		}
		return getEnergy();
	}
}

#endif

