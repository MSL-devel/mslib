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


#include "BaselineInteraction.h"

using namespace MSL;
using namespace std;

const string BaselineInteraction::typeName = "BASELINE";

BaselineInteraction::BaselineInteraction() {
	setup(NULL,0.0);
}

BaselineInteraction::BaselineInteraction(Atom & _d1, double _energy) {
	setup (&_d1,_energy);
}

BaselineInteraction::BaselineInteraction(const BaselineInteraction & _interaction) {
	setup(NULL, 0.0);
	copy(_interaction);
}

BaselineInteraction::~BaselineInteraction() {
}


void BaselineInteraction::setup(Atom * _pA1,double _energy) {
	pAtoms.push_back(_pA1);
	params.push_back(_energy);
}

void BaselineInteraction::copy(const BaselineInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;	
}

void BaselineInteraction::printParameters() {
	if(pAtoms.size() > 0 && pAtoms[0] && params.size() > 0) {
		cout << " ResName " << pAtoms[0]->getResidueName() << endl;
		cout << " atomName " << pAtoms[0]->getName() << endl;
		cout << " energy " << params[0] << endl;
	} else {
		cout << "No atoms defined in interaction" << endl;
	}
}

double BaselineInteraction::getEnergy(double _dummy, std::vector<double> *paramDerivatives) {
	cerr << "WARNING 12334: BaselineInteraction::getEnergy(double _dummy) is not implemented" << endl;
	return 100.0;

}
std::vector<double> BaselineInteraction::getEnergyGrad(){
	cerr << "WARNING 12334: BaselineInteraction::getEnergy(double _dummy) is not implemented" << endl;
	std::vector<double> foo;
	return foo;
}

double BaselineInteraction::getEnergy() {
	return(params[0]);
}

