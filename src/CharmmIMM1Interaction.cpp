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


#include "CharmmIMM1Interaction.h"

using namespace MSL;
using namespace std;


const string CharmmIMM1Interaction::typeName = "CHARMM_IMM1";

CharmmIMM1Interaction::CharmmIMM1Interaction() {
	setup(NULL, NULL);
}

CharmmIMM1Interaction::CharmmIMM1Interaction(Atom& _pA1, Atom& _pA2, vector<double>& _imm1W, vector<double>& _imm1C, double _halfThickness, double _exponent) {
	setup (&_pA1,&_pA2);
	setParams (_imm1W,_imm1C,_halfThickness,_exponent);
}

CharmmIMM1Interaction::CharmmIMM1Interaction(const CharmmIMM1Interaction & _interaction) {
	deletePointers();
	copy(_interaction);
}

void CharmmIMM1Interaction::deletePointers() {
	if(pEEFW) {
		delete pEEFW;
		pEEFW = NULL;
	}
	if(pEEFC) {
		delete pEEFC;
		pEEFC = NULL;
	}
}

CharmmIMM1Interaction::~CharmmIMM1Interaction() {
	deletePointers();
}




void CharmmIMM1Interaction::setup( Atom* _pA1, Atom* _pA2) {
	pAtoms = vector<Atom*> (2, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2);	
	params = vector<double>(2, 0.0);
	pEEFW = new CharmmEEF1Interaction;
	pEEFW->setAtoms(*_pA1, *_pA2);	
	pEEFC = new CharmmEEF1Interaction;
	pEEFC->setAtoms(*_pA1, *_pA2);	
	
}

void CharmmIMM1Interaction::copy(const CharmmIMM1Interaction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	pEEFW = new CharmmEEF1Interaction(*(_interaction.pEEFW));
	pEEFC = new CharmmEEF1Interaction(*(_interaction.pEEFC));
}

