/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
			Scwrl4HBondInteraction(Atom & _d1, Atom & _d2,Atom & _a1, Atom & _a2,Atom& _a3,double _dist, double _ang,double _e1_dihe,double _e2_dihe,double _scalingFactor = 1.0);

			// should implement an operator= as well 
			Scwrl4HBondInteraction(const Scwrl4HBondInteraction & _interaction);
			~Scwrl4HBondInteraction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params);
			void setParams(double _dist, double _ang, double _e1_dihe, double _e2_dihe);
			//double getMinD() const;
			//double getConstant() const;
			std::vector<double> getParams() const;
			
			double getEnergy();
			double getEnergy(double _dummy);

			std::string toString() const;

			//unsigned int getType() const;
			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, Scwrl4HBondInteraction & _term) {_os << _term.toString(); return _os;};
			bool isSelected (std::string _sele1, std::string _sele2) const;
			bool isActive () const;
			double getW() ;
			void setScalingFactor(double _scalingFactor);
			double getScalingFactor() const;
			static void printParameters();
			static void setDebugFlagOn(bool _debugFlagOn = true);

					
		private:
			void setup(Atom * _d1, Atom * _d2,Atom * _a1, Atom * _a2,Atom * _a3,double _d, double _e0,double _e1,double _e2,double _scalingFactor);
			void copy(const Scwrl4HBondInteraction & _interaction);
			double distance;

			//static const unsigned int type = 2;
			static const double d0;
			static const double sig_d;
			static const double B;
			static const double cos_alpha_max;
			static const double cos_beta_max;
			static const double denominator;
			static bool debugFlagOn;
			double scalingFactor;
			static const std::string typeName;
	};

	inline void Scwrl4HBondInteraction::setParams(std::vector<double> _params) { if (_params.size() != 4) {std::cerr << "ERROR 91235: invalid number of parameters in inline void Scwrl4HBondInteraction::setParams(std::vector<double> _params)" << std::endl; exit(91235);} params = _params;}
	inline void Scwrl4HBondInteraction::setParams(double _dist, double _ang, double _e1_dihe, double _e2_dihe) {params[0] = _dist; params[1] = _ang; params[2] = _e1_dihe; params[3] = _e2_dihe; }
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
	inline std::string Scwrl4HBondInteraction::toString() const { char c [1000]; sprintf(c, "SCWRL4 HBOND %s %s %s %s %s %9.4f %9.4f %9.4f %9.4f %20.6f", pAtoms[0]->toString().c_str(),pAtoms[1]->toString().c_str(),pAtoms[2]->toString().c_str(),pAtoms[3]->toString().c_str(),pAtoms[4]->toString().c_str(), params[0], params[1], params[2], params[3],energy); return (std::string)c; };
	//inline unsigned int Scwrl4HBondInteraction::getType() const {return type;}
	inline std::string Scwrl4HBondInteraction::getName() const {return typeName;}
	inline double Scwrl4HBondInteraction::getScalingFactor() const {return scalingFactor;}
	inline void Scwrl4HBondInteraction::setScalingFactor(double _scalingFactor) {scalingFactor = _scalingFactor;}
}

#endif

