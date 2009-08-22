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

#include "CharmmUreyBradleyInteraction.h"

const string CharmmUreyBradleyInteraction::typeName = "CHARMM_U-BR";

CharmmUreyBradleyInteraction::CharmmUreyBradleyInteraction() {
	setup(NULL, NULL, 0.0, 0.0);
}

CharmmUreyBradleyInteraction::CharmmUreyBradleyInteraction(Atom & _a1, Atom & _a2, double _Kub, double _S0) {
	setup (&_a1, &_a2, _Kub, _S0);
}

CharmmUreyBradleyInteraction::CharmmUreyBradleyInteraction(const CharmmUreyBradleyInteraction & _interaction) {
	setup(NULL, NULL, 0.0, 0.0);
	copy(_interaction);
}

CharmmUreyBradleyInteraction::~CharmmUreyBradleyInteraction() {
}

void CharmmUreyBradleyInteraction::setup(Atom * _pA1, Atom * _pA2, double _Kub, double _S0) {
	pAtoms = vector<Atom*> (2);
	pAtoms[0] = NULL; pAtoms[1] = NULL;
	setAtoms(*_pA1, *_pA2);	
	params = vector<double>(2);
	params[0] = 0.0;	params[1] = 0.0;
	setParams(_Kub, _S0);
	distance = 0.0;
}

void CharmmUreyBradleyInteraction::copy(const CharmmUreyBradleyInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
}

