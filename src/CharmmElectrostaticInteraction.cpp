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

#include "CharmmElectrostaticInteraction.h"

const string CharmmElectrostaticInteraction::typeName = "CHARMM_ELEC";

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction() {
	setup(NULL, NULL, false);
}

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction(Atom & _a1, Atom & _a2, bool _is14) {
	setup (&_a1, &_a2, _is14);
}

CharmmElectrostaticInteraction::CharmmElectrostaticInteraction(const CharmmElectrostaticInteraction & _interaction) {
	setup(NULL, NULL, false);
	copy(_interaction);
}

CharmmElectrostaticInteraction::~CharmmElectrostaticInteraction() {
}

void CharmmElectrostaticInteraction::setup(Atom * _pA1, Atom * _pA2, bool _is14) {
	is14 = _is14;

	pAtoms = vector<Atom*> (2);
	pAtoms[0] = NULL; pAtoms[1] = NULL;
	setAtoms(*_pA1, *_pA2);	
	distance = 0.0;
	update();
}

void CharmmElectrostaticInteraction::copy(const CharmmElectrostaticInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	update();
}

void CharmmElectrostaticInteraction::update() {
	params = vector<double>(2, 1.0);
	params[0] = CharmmEnergy::instance()->getDielectricConstant(); //dielectric
	useRiel = CharmmEnergy::instance()->getUseRdielectric();
	if (pAtoms[0] != NULL && pAtoms[1] != NULL) {
		if (is14) {
			// rescaling factor for 1-4 interactions, it was preset to 1
			params[1] = CharmmEnergy::instance()->getElec14factor();
		}
		Kq_q1_q1_rescal_over_diel = CharmmEnergy::Kq * pAtoms[0]->getCharge() * pAtoms[1]->getCharge() * params[1] / params[0];
	} else {
		Kq_q1_q1_rescal_over_diel = 0.0;
	}
}
