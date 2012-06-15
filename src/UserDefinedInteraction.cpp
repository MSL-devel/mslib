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

#include "UserDefinedInteraction.h"

using namespace MSL;
using namespace std;


//const string UserDefinedInteraction::typeName = "USER_DEFINED";

UserDefinedInteraction::UserDefinedInteraction() {
	setup(NULL, NULL, "USER_DEFINED");
}

UserDefinedInteraction::UserDefinedInteraction(Atom & _a1, Atom & _a2, string _type) {
	setup (&_a1, &_a2, _type);
}

UserDefinedInteraction::UserDefinedInteraction(const UserDefinedInteraction & _interaction) {
	setup(NULL, NULL, "USER_DEFINED");
	copy(_interaction);
}

UserDefinedInteraction::~UserDefinedInteraction() {
}

void UserDefinedInteraction::setup(Atom * _pA1, Atom * _pA2, string _type) {

	typeName = _type;
	pAtoms = vector<Atom*> (2);
	pAtoms[0] = NULL; pAtoms[1] = NULL;
	setAtoms(*_pA1, *_pA2);	
}

void UserDefinedInteraction::copy(const UserDefinedInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	typeName = _interaction.typeName;
}


double UserDefinedInteraction::getEnergy(double _distance,std::vector<double> *_dd) {
	// GRADIENT is not implemented
	string name1 = pAtoms[0]->getResidueName()+":"+pAtoms[0]->getName();
	string name2 = pAtoms[1]->getResidueName()+":"+pAtoms[1]->getName();
	return UserDefinedEnergy::instance()->getTwoBodyPotentialValue(typeName,name1,name2,_distance);
}
