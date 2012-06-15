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

#ifndef USERDEFINEDINTERACTION_H
#define USERDEFINEDINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "TwoBodyInteraction.h"
#include "UserDefinedEnergy.h"


namespace MSL { 
class UserDefinedInteraction: public TwoBodyInteraction {

	/*******************************************************
	 *   Inherits from TwoBodyInteraction (a prototype object
	 *   for the interaction of two atoms)
	 *******************************************************/

	public:
		UserDefinedInteraction();
		UserDefinedInteraction(Atom & _a1, Atom & _a2, std::string _type);

		// add an operator= 
		UserDefinedInteraction(const UserDefinedInteraction & _interaction);
		~UserDefinedInteraction();




		double getEnergy();
		double getEnergy(std::vector<double> *_dd);
		double getEnergy(double _distance, std::vector<double> *_dd=NULL);
		std::vector<double> getEnergyGrad();

		friend std::ostream & operator<<(std::ostream &_os, UserDefinedInteraction & _term) {_os << _term.toString(); return _os;};
		std::string toString() ;

		std::string getName() const;
		void setName(std::string _type);
		std::pair<double,std::vector<double> > partialDerivative();
		
	private:
		void setup(Atom * _a1, Atom * _a2, std::string _type);
		void copy(const UserDefinedInteraction & _interaction);

		std::string typeName;
		

};

inline double UserDefinedInteraction::getEnergy() {
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}

inline double UserDefinedInteraction::getEnergy(std::vector<double>* _dd) {
	return getEnergy(pAtoms[0]->distance(*pAtoms[1]));
}


inline std::string UserDefinedInteraction::toString() { char c [1000]; sprintf(c, "%s %s %s %9.4f %9.4f %9.4f %20.6f", typeName.c_str(), pAtoms[0]->toString().c_str(), pAtoms[1]->toString().c_str(), params[0], params[1], pAtoms[0]->distance(*pAtoms[1]), getEnergy()); return (std::string)c; };
inline std::string UserDefinedInteraction::getName() const {return typeName;}
inline void UserDefinedInteraction::setName(std::string _type) { typeName = _type; }

inline std::vector<double> UserDefinedInteraction::getEnergyGrad(){
   std::vector<double> tmp;
	return tmp;
}

inline std::pair<double,std::vector<double> > UserDefinedInteraction::partialDerivative() {
	std::pair<double, std::vector<double> > partials;
	std::cerr << "UserDefinedInteraction::partialDerivative is not implemented" << std::endl;
	partials.first = 0;
	partials.second.resize(pAtoms.size() * 3, 0.0);
	return partials;
}
}

#endif

