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

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <vector>
#include <iostream>
#include <stdio.h>
#include <cstdlib>



/*! \brief This class is used to enumerate combinations of states
 */


namespace MSL { 
class Enumerator {
	public:
		Enumerator();
		Enumerator(std::vector<unsigned int> _states);
		Enumerator(std::vector<std::vector<unsigned int> > & _values);
		Enumerator(Enumerator & _enum);
		~Enumerator();

		std::vector<unsigned int> & operator[](size_t n); // returns value of state
		std::vector<unsigned int> & operator()(size_t n); // returns state


		std::vector<unsigned int> & getValue(size_t n); // returns value of state
		std::vector<unsigned int> & getState(size_t n); // returns state

		void setStates(std::vector<unsigned int> states);
		
		unsigned int size() const;
		unsigned int numberOfVariables() const;
		unsigned int getStateIndex(std::vector<unsigned int> _states) const; /*! \brief Returns the index of a state (i.e. 3, 5, 8, 7 => 385992) */

		void setMaxLimit(unsigned int _limit);
		unsigned int getMaxLimit() const;

		void setValues(std::vector<std::vector<unsigned int> > & _values);

	private:
		void calcEnumeration();
		void calcEnumerationValues();

		std::vector<unsigned int> statesPerElement;
		std::vector<std::vector<unsigned int> > enumerations;
		std::vector<std::vector<unsigned int > > values;
		std::vector<std::vector<unsigned int > > enumeratedValues;
		unsigned int combinatorialSize;
		unsigned int maxLimit;
		bool valueSet_flag;
};

inline std::vector<unsigned int> & Enumerator::operator[](size_t _n) {
	return getValue(_n);
}

inline std::vector<unsigned int> & Enumerator::operator()(size_t _n) {
	return getState(_n);
}

inline std::vector<unsigned int> & Enumerator::getValue(size_t _n) {
	if(!valueSet_flag) {
		return enumerations[_n];
	} 
	return enumeratedValues[_n];
}

inline std::vector<unsigned int> & Enumerator::getState(size_t _n) {
	return enumerations[_n];
}

inline void Enumerator::setMaxLimit(unsigned int _limit) {
	maxLimit = _limit;
}

inline unsigned int Enumerator::getMaxLimit() const {
	return maxLimit;
}

inline unsigned int Enumerator::size() const {
	return combinatorialSize;
}

inline unsigned int Enumerator::numberOfVariables() const {
	return statesPerElement.size();
}



}

#endif

