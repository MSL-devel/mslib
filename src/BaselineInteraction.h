/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2011 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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


#ifndef BASELINEINTERACTION_H
#define BASELINEINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "OneBodyInteraction.h"
#include "Atom.h"


namespace MSL { 
	class BaselineInteraction: public OneBodyInteraction {

		/*******************************************************
		 *   Inherits from OneBodyInteraction (a prototype object
		 *   for the interaction of one atom)
		 *******************************************************/

		public:
			BaselineInteraction();
			BaselineInteraction(Atom & _d1, double _energy);

			// should implement an operator= as well 
			BaselineInteraction(const BaselineInteraction & _interaction);
			~BaselineInteraction();

			double getEnergy(double _param, std::vector<double> *paramDerivatives=NULL);
			std::vector<double> getEnergyGrad();

			double getEnergy();
			std::string toString() const;
			void printParameters();

			//unsigned int getType() const;
			std::string getName() const;
			friend std::ostream & operator<<(std::ostream &_os, BaselineInteraction & _term) {_os << _term.toString(); return _os;};

					
		private:
			void setup(Atom * _d1, double _energy);
			void copy(const BaselineInteraction & _interaction);
			static const std::string typeName;
	};

	inline std::string BaselineInteraction::toString() const { 
		if(pAtoms.size() && pAtoms[0]) {
			char c [1000]; 
			sprintf(c, "BASELINE %s %s %9.4f", pAtoms[0]->toString().c_str(),pAtoms[0]->getResidueName().c_str(),energy); 
			return (std::string)c; 
		}
		return "";
	};
	inline std::string BaselineInteraction::getName() const {return typeName;}
}

#endif

