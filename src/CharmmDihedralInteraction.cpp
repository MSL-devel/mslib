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

#include "CharmmDihedralInteraction.h"

using namespace MSL;
using namespace std;


const string CharmmDihedralInteraction::typeName = "CHARMM_DIHE";

CharmmDihedralInteraction::CharmmDihedralInteraction() {
	vector<vector<double> > multi(1, vector<double>(3, 0.0));
	multi[0][1] = 1;
	setup(NULL, NULL, NULL, NULL, multi);
}

CharmmDihedralInteraction::CharmmDihedralInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, double _Kchi, double _N, double _deltaRadians) {
	vector<vector<double> > multi(1, vector<double>(3, 0.0));
	multi[0][0] = _Kchi;
	multi[0][1] = _N;
	multi[0][2] = _deltaRadians;
	setup (&_a1, &_a2, &_a3, &_a4, multi);
}

CharmmDihedralInteraction::CharmmDihedralInteraction(Atom & _a1, Atom & _a2, Atom & _a3, Atom & _a4, vector<vector <double > > _multipleParams) {
	setup (&_a1, &_a2, &_a3, &_a4, _multipleParams);
}

CharmmDihedralInteraction::CharmmDihedralInteraction(const CharmmDihedralInteraction & _interaction) {
	setup(NULL, NULL, NULL, NULL, vector<vector<double> >(1, vector<double>(3, 0.0)));
	copy(_interaction);
}

CharmmDihedralInteraction::~CharmmDihedralInteraction() {
}

void CharmmDihedralInteraction::setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, vector<vector<double> >  _params) {
	pAtoms = vector<Atom*>(4, (Atom*)NULL);
	setAtoms(*_pA1, *_pA2, *_pA3, *_pA4);	
	setParams(_params);
}

void CharmmDihedralInteraction::copy(const CharmmDihedralInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
}

std::vector<double> CharmmDihedralInteraction::getEnergyGrad(Atom& a1, Atom& a2, Atom& a3, Atom& a4, double Kchi, double N, double deltaRadians) {
	std::vector<double> dd;
	double angleRadians = CartesianGeometry::dihedralDerivative(a1.getCoor(), a2.getCoor(), a3.getCoor(), a4.getCoor(),&dd);
	CharmmEnergy::instance()->dihedralEnerGrad(dd, angleRadians, Kchi, N, deltaRadians);
	return dd;
}
