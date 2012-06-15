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

#include "CharmmImproperInteraction.h"

using namespace MSL;
using namespace std;


const string CharmmImproperInteraction::typeName = "CHARMM_IMPR";

CharmmImproperInteraction::CharmmImproperInteraction() {
	setup(NULL, NULL, NULL, NULL, 0.0, 0.0);
}

CharmmImproperInteraction::CharmmImproperInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, double _Kpsi, double _Psi0Radians) {
	setup (&_a1, &_a2, &_a3, &_a4, _Kpsi, _Psi0Radians);
}

CharmmImproperInteraction::CharmmImproperInteraction(const CharmmImproperInteraction & _interaction) {
	setup(NULL, NULL, NULL, NULL, 0.0, 0.0);
	copy(_interaction);
}

CharmmImproperInteraction::~CharmmImproperInteraction() {
}

void CharmmImproperInteraction::setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, double _Kpsi, double _Psi0Radians) {
	pAtoms = vector<Atom*>(4, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2, *_pA3, *_pA4);	
	params = vector<double>(2, 0.0);
	setParams(_Kpsi, _Psi0Radians);
}

void CharmmImproperInteraction::copy(const CharmmImproperInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
}

std::vector<double> CharmmImproperInteraction::getEnergyGrad(Atom& a1, Atom& a2, Atom& a3, Atom& a4, double Kpsi, double Psi0Radians) {
	std::vector<double> id;
	double angleRadians = CartesianGeometry::dihedralDerivative(a1.getCoor(), a2.getCoor(), a3.getCoor(), a4.getCoor(), &id);
	CharmmEnergy::instance()->springGrad(id, angleRadians, Kpsi, Psi0Radians);
	return id;
}
