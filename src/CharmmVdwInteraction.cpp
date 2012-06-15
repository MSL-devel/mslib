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

#include "CharmmVdwInteraction.h"

using namespace MSL;
using namespace std;


const string CharmmVdwInteraction::typeName = "CHARMM_VDW";

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
}

void CharmmVdwInteraction::copy(const CharmmVdwInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	useNonBondCutoffs = _interaction.useNonBondCutoffs;
	nonBondCutoffOn = _interaction.nonBondCutoffOn;
	nonBondCutoffOff = _interaction.nonBondCutoffOff;
}

vector<double> CharmmVdwInteraction::getEnergyGrad(Atom& _a1, Atom& _a2, double _rmin, double _emin) {
	vector<double> dd;
	double distance = CartesianGeometry::distanceDerivative(_a1.getCoor(), _a2.getCoor(),&dd);
	CharmmEnergy::instance()->LJGrad(dd, distance, _rmin, _emin);
	return dd;
}
