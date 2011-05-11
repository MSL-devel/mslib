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

#include "CharmmVdwInteraction.h"

using namespace MSL;
using namespace std;



CharmmVdwInteraction::CharmmVdwInteraction() {
	setup(NULL, NULL, 0.0, 0.0);
}

CharmmVdwInteraction::CharmmVdwInteraction(Atom & _a1, Atom & _a2, double _rmin, double _Emin) {
	setup (&_a1, &_a2, _rmin, _Emin);
}

CharmmVdwInteraction::CharmmVdwInteraction(const CharmmVdwInteraction & _interaction) {
	setup(NULL, NULL, 0.0, 0.0);
	copy(_interaction);
}

CharmmVdwInteraction::~CharmmVdwInteraction() {
}

void CharmmVdwInteraction::setup(Atom * _pA1, Atom * _pA2, double _rmin, double _Emin) {
	pAtoms = vector<Atom*>(2, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2);	
	params = vector<double>(2, 0.0);
	setParams(_rmin, _Emin);
	useNonBondCutoffs = false;
	nonBondCutoffOn = 997;
	nonBondCutoffOff = 998;
	typeName = "CHARMM_VDW";
}

void CharmmVdwInteraction::copy(const CharmmVdwInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	typeName = _interaction.typeName;
	useNonBondCutoffs = _interaction.useNonBondCutoffs;
	nonBondCutoffOn = _interaction.nonBondCutoffOn;
	nonBondCutoffOff = _interaction.nonBondCutoffOff;
}

std::vector<double> CharmmVdwInteraction::getEnergyGrad(Atom& a1, Atom& a2, double rmin, double Emin) {
	std::vector<double> dd;
	CartesianGeometry::distanceDerivative(a1.getCoor(), a2.getCoor(),&dd);
	CharmmEnergy::instance()->LJGrad(dd, a1.distance(a2), rmin, Emin);
	return dd;
}
