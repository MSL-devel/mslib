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


#include "EZpotentialInteraction.h"

using namespace MSL;
using namespace std;


const string EZpotentialInteraction::typeName = "EZ_POTENTIAL";

EZpotentialInteraction::EZpotentialInteraction() {
	setup(NULL, vector<double>(3, 0.0), true);
}

EZpotentialInteraction::EZpotentialInteraction(Atom & _pA1, std::vector<double> _params, bool _sigmoidalFunction) {
	setup (&_pA1, _params, _sigmoidalFunction);
}

EZpotentialInteraction::EZpotentialInteraction(const EZpotentialInteraction & _interaction) {
	setup(NULL, vector<double>(3, 0.0), true);
	copy(_interaction);
}

EZpotentialInteraction::~EZpotentialInteraction() {
}

void EZpotentialInteraction::setup(Atom * _pA1, vector<double> _params, bool _sigmoidalFunction) {
	setParams(_params);
	pAtoms = vector<Atom*> (1, (Atom*)NULL);
	setAtoms(*_pA1);	
	isSigmoidal_flag = _sigmoidalFunction;
}

void EZpotentialInteraction::copy(const EZpotentialInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	isSigmoidal_flag = _interaction.isSigmoidal_flag;
}

double EZpotentialInteraction::getEnergy() {
	return getEnergy(pAtoms[0]->getZ(), params, isSigmoidal_flag);
}

double EZpotentialInteraction::getEnergy(double _Zcoor, const vector<double> & _param, bool _sigmoidalFunction) const { 
	if(_Zcoor < 0) 
		_Zcoor = -1*_Zcoor;
	if (_Zcoor > 25.0) {
		_Zcoor = 25.0;
	}
	if (_sigmoidalFunction) {
		return (_param[0] / ( 1 + (pow((_Zcoor/_param[1]), _param[2]))));
	} else {
		return _param[2]*exp(-((_Zcoor-_param[1])*(_Zcoor-_param[1])) / (2*_param[0]*_param[0]));
	}
}


