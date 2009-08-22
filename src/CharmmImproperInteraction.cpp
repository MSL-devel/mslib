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

#include "CharmmImproperInteraction.h"

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
	pAtoms = vector<Atom*> (4);
	pAtoms[0] = NULL; pAtoms[1] = NULL;pAtoms[2] = NULL;pAtoms[3] = NULL;
	setAtoms(*_pA1, *_pA2, *_pA3, *_pA4);	
	params = vector<double>(2);
	params[0] = 0.0;	params[1] = 0.0;	
	setParams(_Kpsi, _Psi0Radians);
	angle = 0.0;
}

void CharmmImproperInteraction::copy(const CharmmImproperInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
}

