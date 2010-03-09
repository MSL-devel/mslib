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


#ifndef CHARMMEEF1REFINTERACTION_H
#define CHARMMEEF1REFINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "CharmmEnergy.h"


namespace MSL { 
	class CharmmEEF1RefInteraction: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from TwoBodyInteraction (a prototype object
		 *   for the interaction of two atoms)
		 *******************************************************/

		public:
			CharmmEEF1RefInteraction();
			CharmmEEF1RefInteraction(Atom & _a1, double _ref);

			// should implement an operator= as well 
			CharmmEEF1RefInteraction(const CharmmEEF1RefInteraction & _interaction);
			~CharmmEEF1RefInteraction();

			/* setting and getting the parameters */
			void setParams(std::vector<double> _params);
			void setParams(double _ref);
			//double getMinD() const;
			//double getConstant() const;
			std::vector<double> getParams() const;
			
			double getEnergy();
			double getEnergy(double _phony);

			friend std::ostream & operator<<(std::ostream &_os, CharmmEEF1RefInteraction & _term) {_os << _term.toString(); return _os;};
			std::string toString() const;

			//unsigned int getType() const;
			std::string getName() const;
			
		private:
			void setup(Atom * _a1, double _ref);
			void copy(const CharmmEEF1RefInteraction & _interaction);

			//static const unsigned int type = 2;
			static const std::string typeName;
			
	};

	inline void CharmmEEF1RefInteraction::setParams(std::vector<double> _params) { if (_params.size() != 1) {std::cerr << "ERROR 31235: invalid number of parameters in inline void CharmmEEF1RefInteraction::setParams(std::vector<double> _params)" << std::endl; exit(31235);} params = _params;}
	inline void CharmmEEF1RefInteraction::setParams(double _ref) {params[0] = _ref;}
	inline std::vector<double> CharmmEEF1RefInteraction::getParams() const {return params;};
	inline double CharmmEEF1RefInteraction::getEnergy() {return params[0];}
	inline double CharmmEEF1RefInteraction::getEnergy(double _phony) {return getEnergy();}
	inline std::string CharmmEEF1RefInteraction::toString() const { char c [1000]; sprintf(c, "CHARMM EEF1REF %s %9.4f", pAtoms[0]->toString().c_str(), params[0]); return (std::string)c; };
	//inline unsigned int CharmmEEF1RefInteraction::getType() const {return type;}
	inline std::string CharmmEEF1RefInteraction::getName() const {return typeName;}

}

#endif

