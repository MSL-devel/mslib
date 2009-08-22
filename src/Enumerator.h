/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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


using namespace std;

/*! \brief This class is used to enumerate combinations of states
 */


class Enumerator {
	public:
		Enumerator();
		Enumerator(vector<unsigned int> _states);
		Enumerator(Enumerator & _enum);
		~Enumerator();

		vector<unsigned int> & operator[](size_t n);
		void setStates(vector<unsigned int> states);
		
		unsigned int size() const;
		unsigned int numberOfVariables() const;
		unsigned int getStateIndex(vector<unsigned int> _states) const; /*! \brief Returns the index of a state (i.e. 3, 5, 8, 7 => 385992) */

		void setMaxLimit(unsigned int _limit);
		unsigned int getMaxLimit() const;

	private:
		void calcEnumeration();

		vector<unsigned int> statesPerElement;
		vector<vector<unsigned int> > enumerations;
		unsigned int combinatorialSize;
		unsigned int maxLimit;
};

inline vector<unsigned int> & Enumerator::operator[](size_t _n) {
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



#endif

