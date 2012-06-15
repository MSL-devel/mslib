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


#ifndef EZPOTENTIALINTERACTION_H
#define EZPOTENTIALINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"

/*******************************************************
 *  This is an implementation of the EZ empirical membrane 
 *  insertion potential.  If used please cite:
 *
 *     Senes A, Chadi DC, Law PB, Walters RF, Nanda V, DeGrado WF. 
 *     "E(z), a depth-dependent potential for assessing the energies 
 *     of insertion of amino acid side-chains into membranes, 
 *     derivation and applications to determining the orientation 
 *     of transmembrane and interfacial helices."
 *     J Mol Biol. 2007 366(2), 436-48
 *
 *******************************************************/

namespace MSL { 
	class EZpotentialInteraction: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from OneBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			EZpotentialInteraction();
			EZpotentialInteraction(Atom & _pA1, std::vector<double> _params, bool _sigmoidalFunction);

			// should implement an operator= as well 
			EZpotentialInteraction(const EZpotentialInteraction & _interaction);
			~EZpotentialInteraction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params);
			std::vector<double> getParams() const;
			void setSigmoidalFunction(bool _flag);
			bool getSigmoidalFunction() const;
			
			double getEnergy();
			double getEnergy(double _Zcoor, const std::vector<double> & _param, bool _sigmoidalFunction) const;
			double getEnergy(std::vector<double> *paramDerivatives); // derivatives not implemented
			double getEnergy(double _param, std::vector<double> *paramDerivatives=NULL); // derivatives not implemented
			std::vector<double> getEnergyGrad(); // derivatives not implemented

			friend std::ostream & operator<<(std::ostream &_os, EZpotentialInteraction & _term) {_os << _term.toString(); return _os;};
			std::string toString() ;

			std::string getName() const;
			std::pair<double,std::vector<double> > partialDerivative();
			
		private:
			void setup(Atom * _pA1, std::vector<double> _params, bool _sigmoidalFunction);
			void copy(const EZpotentialInteraction & _interaction);

			//static const unsigned int type = 2;
			static const std::string typeName;
			bool isSigmoidal_flag;
			
	};

	inline void EZpotentialInteraction::setParams(std::vector<double> _params) { if (_params.size() != 3) {std::cerr << "ERROR 12357: invalid number of parameters in inline void EZpotentialInteraction::setParams(std::vector<double> _params, bool _sigmoidalFunction) = " << _params.size() << std::endl; exit(12357);} params = _params;}
	inline std::vector<double> EZpotentialInteraction::getParams() const {return params;};
	inline std::string EZpotentialInteraction::toString() { char c [1000]; sprintf(c, "%s %s %9.4f %9.4f %9.4f %u %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), params[0], params[1], params[2], isSigmoidal_flag, getEnergy()); return (std::string)c; };
	inline void EZpotentialInteraction::setSigmoidalFunction(bool _flag) {
		isSigmoidal_flag = _flag;
	}
	inline bool EZpotentialInteraction::getSigmoidalFunction() const {
		return isSigmoidal_flag;
	}

	inline std::string EZpotentialInteraction::getName() const {return typeName;}

	inline std::pair<double,std::vector<double> > EZpotentialInteraction::partialDerivative() {
		std::cerr << "ERROR 38984: std::pair<double,std::vector<double> > EZpotentialInteraction::partialDerivative() is not implemented" << std::endl;
		exit(38984);
	}
	inline double EZpotentialInteraction::getEnergy(std::vector<double> *_dd) {
		std::cerr << "ERROR 38989: double EZpotentialInteraction::getEnergy(std::vector<double> *_dd) is not implemented" << std::endl;
		exit(38989);
	}
	inline double EZpotentialInteraction::getEnergy(double _param, std::vector<double> *paramDerivatives) {
		std::cerr << "ERROR 38994: double EZpotentialInteraction::getEnergy(double _param, std::vector<double> *paramDerivatives) is not implemented" << std::endl;
		exit(38994);
	}
	inline std::vector<double> EZpotentialInteraction::getEnergyGrad() {
		std::cerr << "ERROR 38999: std::vector<double> EZpotentialInteraction::getEnergyGrad() is not implemented" << std::endl;
		exit(38999);
	}
}

#endif

