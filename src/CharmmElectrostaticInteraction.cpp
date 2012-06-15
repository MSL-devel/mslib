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


#include "CharmmElectrostaticInteraction.h"

using namespace MSL;
using namespace std;


const string CharmmElectrostaticInteraction::typeName = "CHARMM_ELEC";

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction() {
	setup(NULL, NULL, 1.0, 1.0, false);
}

//CharmmElectrostaticInteraction::CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, bool _is14) {
//	setup (&_a1, &_a2, _is14);
//}

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, double _dielectricConstant, double _14rescaling, bool _useRdielectric) {
	setup (&_a1, &_a2, _dielectricConstant, _14rescaling, _useRdielectric);
}

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction(const CharmmElectrostaticInteraction & _interaction) {
	setup(NULL, NULL, 1.0, 1.0, false);
	copy(_interaction);
}

CharmmElectrostaticInteraction::~CharmmElectrostaticInteraction() {
}

void CharmmElectrostaticInteraction::setup(Atom * _pA1, Atom * _pA2, double _dielectricConstant, double _14rescaling, bool _useRdielectric) {
	//is14 = _is14;
	pAtoms = vector<Atom*>(2, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2);	
	params = vector<double>(2, 1.0);
	params[0] = _dielectricConstant;
	params[1] = _14rescaling;
	useRiel = _useRdielectric;
	useNonBondCutoffs = false;
	update();
	nonBondCutoffOn = 0.0;
	nonBondCutoffOff = 0.0;
}
/*
void CharmmElectrostaticInteraction::setup(Atom * _pA1, Atom * _pA2, bool _is14) {
	is14 = _is14;
	pAtoms = vector<Atom*>(2, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2);	
	update();
}
*/

void CharmmElectrostaticInteraction::copy(const CharmmElectrostaticInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	useRiel = _interaction.useRiel;
	Kq_q1_q1_rescal_over_diel = _interaction.Kq_q1_q1_rescal_over_diel;
}

void CharmmElectrostaticInteraction::update() {
	//params = vector<double>(2, 1.0);
	//params[0] = CharmmEnergy::instance()->getDielectricConstant(); //dielectric
	//useRiel = CharmmEnergy::instance()->getUseRdielectric();
	if (pAtoms[0] != NULL && pAtoms[1] != NULL) {
		//if (is14) {
		//	// rescaling factor for 1-4 interactions, it was preset to 1
		//	params[1] = CharmmEnergy::instance()->getElec14factor();
		//}
		Kq_q1_q1_rescal_over_diel = CharmmEnergy::Kq * pAtoms[0]->getCharge() * pAtoms[1]->getCharge() * params[1] / params[0];
	} else {
		Kq_q1_q1_rescal_over_diel = 0.0;
	}
}
std::vector<double> CharmmElectrostaticInteraction::getEnergyGrad(Atom& _a1, Atom& _a2, bool _is14) {
	std::vector<double> dd;
	double distance = CartesianGeometry::distanceDerivative(_a1.getCoor(), _a2.getCoor(),&dd);
	//CharmmEnergy::instance()->coulombEnerGrad(dd, a1.distance(a2), a1.getCharge(), a2.getCharge(), _is14 ? CharmmEnergy::instance()->getElec14factor() : 1);
	CharmmEnergy::instance()->coulombEnerGrad(dd, distance, Kq_q1_q1_rescal_over_diel, useRiel);
	return dd;
}


