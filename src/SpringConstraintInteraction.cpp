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

#include "SpringConstraintInteraction.h"

using namespace MSL;
using namespace std;

const string SpringConstraintInteraction::typeName = "SPRING_CONSTRAINT";

SpringConstraintInteraction::SpringConstraintInteraction() {
	setup(NULL,0.0, 0.0);
}

SpringConstraintInteraction::SpringConstraintInteraction(Atom & _a1, double _Kb, double _b0) {
	setup (&_a1, _Kb, _b0);
}

SpringConstraintInteraction::SpringConstraintInteraction(const SpringConstraintInteraction & _interaction) {
	setup(NULL,0.0, 0.0);
	copy(_interaction);
}

SpringConstraintInteraction::~SpringConstraintInteraction() {
}

void SpringConstraintInteraction::setup(Atom * _pA1, double _Kb, double _b0) {
	pAtoms = vector<Atom*> (1, (Atom*)NULL);
	pAtoms[0] = _pA1;
	if(_pA1) {
		referenceCoor = _pA1->getCoor();
	} else {
		referenceCoor.setCoor(0.0,0.0,0.0);
	}
	params = vector<double>(2, 0.0);
	setParams(_Kb, _b0);
}

void SpringConstraintInteraction::copy(const SpringConstraintInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	referenceCoor = _interaction.referenceCoor;
}

std::vector<double> SpringConstraintInteraction::getEnergyGrad(Atom& a1, double Kb, double b0) {
	std::vector<double> dd;
	double distance = CartesianGeometry::distanceDerivative(a1.getCoor(), referenceCoor,&dd);
	CharmmEnergy::instance()->springGrad(dd, distance, Kb, b0);
	return dd;
}
