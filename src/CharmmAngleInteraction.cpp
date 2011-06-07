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

#include "CharmmAngleInteraction.h"

using namespace MSL;
using namespace std;


const string CharmmAngleInteraction::typeName = "CHARMM_ANGL";

CharmmAngleInteraction::CharmmAngleInteraction() {
	setup(NULL, NULL, NULL, 0.0, 0.0);
}

CharmmAngleInteraction::CharmmAngleInteraction(Atom & _a1, Atom & _a2, Atom & _a3, double _Ktheta, double _Theta0Radians) {
	setup (&_a1, &_a2, &_a3, _Ktheta, _Theta0Radians);
}

CharmmAngleInteraction::CharmmAngleInteraction(const CharmmAngleInteraction & _interaction) {
	setup(NULL, NULL, NULL, 0.0, 0.0);
	copy(_interaction);
}

CharmmAngleInteraction::~CharmmAngleInteraction() {
}

void CharmmAngleInteraction::setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, double _Ktheta, double _Theta0Radians) {
	pAtoms = vector<Atom*>(3, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2, *_pA3);	
	params = vector<double>(2, 0.0);
	setParams(_Ktheta, _Theta0Radians);
}

void CharmmAngleInteraction::copy(const CharmmAngleInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
}

std::vector<double> CharmmAngleInteraction::getEnergyGrad(Atom& a1, Atom& a2, Atom& a3, double Ktheta, double Theta0Radians) {
	std::vector<double> ad;
	double angleRadians = CartesianGeometry::angleDerivative(a1.getCoor(), a2.getCoor(), a3.getCoor(),&ad);
	CharmmEnergy::instance()->springGrad(ad, angleRadians, Ktheta, Theta0Radians);
	return ad;
}
