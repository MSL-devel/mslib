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

#include "Enumerator.h"

using namespace MSL;
using namespace std;


Enumerator::Enumerator() {
	maxLimit = 10000000;
	valueSet_flag = false;
}

Enumerator::Enumerator(vector<unsigned int> _states) {
	maxLimit = 10000000;
	valueSet_flag = false;
	setStates(_states);
}

Enumerator::Enumerator(vector<vector<unsigned int> > & _values) {
	maxLimit = 10000000;
	valueSet_flag = false;
	setValues(_values);
}

Enumerator::Enumerator(Enumerator & _enum) {
	maxLimit = _enum.maxLimit;
	enumerations = _enum.enumerations;
	statesPerElement = _enum.statesPerElement;
	combinatorialSize = _enum.combinatorialSize;
	values = _enum.values;
	valueSet_flag = _enum.valueSet_flag;
}

Enumerator::~Enumerator() {
}

void Enumerator::setStates(vector<unsigned int> _states) {
	statesPerElement = _states;
	enumerations.clear();
	combinatorialSize = 1;
	if (statesPerElement.size() == 0) {
		combinatorialSize = 0;
		return;
	}
	for (unsigned int i=0; i<statesPerElement.size(); i++) {
		combinatorialSize *= statesPerElement[i];
	}
	calcEnumeration();
}

void Enumerator::calcEnumeration() {

	if (combinatorialSize > maxLimit) {
		cerr << "ERROR 3175: number of combinations " << combinatorialSize << " exceeds the max limit " << maxLimit << " in void Enumerator::calcEnumeration()" << endl;
		exit(3175);
	} else if (combinatorialSize == 0) {
		return;
	}

	vector<unsigned int> tmp(statesPerElement.size(), 0);
	for (unsigned int i=0; i<combinatorialSize; i++) {
		enumerations.push_back(tmp);
	}

	unsigned int prev = 1;
	for (unsigned int i=0; i<statesPerElement.size(); i++) {
		unsigned int steps = combinatorialSize / (statesPerElement[i] * prev);
		for (unsigned int j=0; j<steps; j++) {
			for (unsigned int k=0; k<statesPerElement[i]; k++) {
				for (unsigned int l=0; l<prev; l++) {
					enumerations[prev*(statesPerElement[i]*j+k)+l][i] = k;
				}
			}
		}
		prev *= statesPerElement[i];
	}
}
		
unsigned int Enumerator::getStateIndex(vector<unsigned int> _states) const {
	/************************************
	 *  Give a state it returns its index
	 ************************************/
	int out = 0;;
	int factor = 1;
	if (_states.size() != statesPerElement.size()) {
		cerr << "ERROR 3179: the number of input positions (" << _states.size() << ") is different that the number of positions of the Enumerator (" << statesPerElement.size() << ") in unsigned int Enumerator:getStateIndex(vector<int> _states) const" << endl;
		exit(3179);
	}
	for (int i=0; i<_states.size(); i++) {
		if (_states[i] >= statesPerElement[i]) {
			cerr << "ERROR 3184: the input state at position " << i << " (" << _states[i] << ") is  larger than the number of possible states of the Enumerator at that position (" << statesPerElement[i] << ") in unsigned int Enumerator:getStateIndex(vector<int> _states) const" << endl;
			exit(3184);
		}
		out += _states[i] * factor;
		factor *= statesPerElement[i];
	}
	return out;
}

void Enumerator::calcEnumerationValues() {

	enumeratedValues.clear();
	for (unsigned int i=0; i<enumerations.size(); i++) {
		enumeratedValues.push_back(vector<unsigned int>());
		for (unsigned int j=0; j<enumerations[i].size(); j++) {
			enumeratedValues[i].push_back(values[j][enumerations[i][j]]);
		}
	}
}

void Enumerator::setValues(vector<vector<unsigned int> > & _values) {
	values = _values;
	vector<unsigned int> states;
	for (int i=0; i<values.size(); i++) {
		states.push_back(values[i].size());
	}
	setStates(states);		
	valueSet_flag = true;
	calcEnumerationValues();
}

