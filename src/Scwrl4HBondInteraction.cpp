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


#include "Scwrl4HBondInteraction.h"
#include <math.h>

using namespace MSL;
using namespace std;


const string Scwrl4HBondInteraction::typeName = "SCWRL4_HBOND";
// parameters from "G.G.Krivov et al,Improved prediction of protein side-chain conformations with SCWRL4"
bool Scwrl4HBondInteraction::debugFlagOn = false;


Scwrl4HBondInteraction::Scwrl4HBondInteraction() {
	setup(NULL, NULL, NULL, NULL, NULL,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
}

Scwrl4HBondInteraction::Scwrl4HBondInteraction(Atom & _d1, Atom & _d2,Atom & _a1, Atom & _a2,Atom& _a3,double _dist, double _ang,double _e1_dihe,double _e2_dihe,
	double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max,double _scalingFactor) {
	setup (&_d1, &_d2,&_a1, &_a2,&_a3,_dist, _ang,_e1_dihe,_e2_dihe,_d0,_sig_d,_B,_alpha_max,_beta_max,_scalingFactor);
}

Scwrl4HBondInteraction::Scwrl4HBondInteraction(const Scwrl4HBondInteraction & _interaction) {
	setup(NULL, NULL, NULL, NULL, NULL,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	copy(_interaction);
}

Scwrl4HBondInteraction::~Scwrl4HBondInteraction() {
}


void Scwrl4HBondInteraction::setup(Atom * _pA1, Atom * _pA2, Atom * _pA3, Atom * _pA4, Atom * _pA5, double _dist, double _ang, double _e1_dihe, double _e2_dihe,double _d0, double _sig_d, double _B, double _alpha_max, double _beta_max,double _scalingFactor) {
	pAtoms.push_back(_pA1); // actual donor 
	pAtoms.push_back(_pA2);  // bonded to donor
	pAtoms.push_back(_pA3);  // actual acceptor
	pAtoms.push_back(_pA4);  // bonded to acceptor
	pAtoms.push_back(_pA5);  // bonded to acceptor_2
	params.push_back(_dist);
	params.push_back(_ang);
	params.push_back(_e1_dihe);
	params.push_back(_e2_dihe);
	params.push_back(_d0);
	params.push_back(_sig_d);
	params.push_back(_B);
	params.push_back(cos(_alpha_max)); // store cos_alpha_max
	params.push_back(cos(_beta_max )); // store cos_beta_max 
	scalingFactor = _scalingFactor;
}

void Scwrl4HBondInteraction::copy(const Scwrl4HBondInteraction & _interaction) {
	pAtoms = _interaction.pAtoms;
	params = _interaction.params;
	scalingFactor = _interaction.scalingFactor;
}
double Scwrl4HBondInteraction::getEnergy(double _d, std::vector<double> *_paramDerivatives) {
		if(_paramDerivatives) {
			return getEnergy(_paramDerivatives);
		}
		return getEnergy();
}
std::vector<double> Scwrl4HBondInteraction::getEnergyGrad(){
		std::vector<double> foo;
		getEnergy(&foo);
		return foo;
}

void Scwrl4HBondInteraction::printParameters() {
	if(params.size() == 9) {
		cout << " d: " << params[0] << endl;
		cout << " ang: " << params[1] << endl;
		cout << " dihe1: " << params[2] << endl;
		cout << " dihe2: " << params[3] << endl;
		cout << " d0: " << params[4] << endl;
		cout << " sig_d: " << params[5] << endl;
		cout << " B: " << params[6] << endl;
		cout << " cos_alpha_max " << params[7] << endl;
		cout << " cos_beta_max " << params[8] << endl;
	} else {
		cerr << "ERROR 12247: Parameters not set" << endl;
	}
}

void Scwrl4HBondInteraction::setDebugFlagOn(bool _debugFlagOn) {
	debugFlagOn = _debugFlagOn;
}

double Scwrl4HBondInteraction::getW() {
	CartesianPoint e0(0,0,0); 
	CartesianPoint e1(0,0,0); 
	CartesianPoint e2(0,0,0); 
	CartesianPoint n(0,0,0); 
	
	// atoms[0] actual donor 
	// atoms[1] bonded to donor
	// atoms[2] actual acceptor
	// atoms[3] bonded to acceptor
	// atoms[4] bonded to acceptor_2 

	e0 = (pAtoms[0]->getCoor() - pAtoms[1]->getCoor());
	e1 = (CartesianGeometry::buildRadians(pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),pAtoms[4]->getCoor(),params[0],params[1],params[2])) - pAtoms[2]->getCoor();
	n = (pAtoms[0]->getCoor()-pAtoms[2]->getCoor());


	double d = pAtoms[0]->distance(*pAtoms[2]);

	double t1 = (params[5] *params[5]) - (d - params[4]) * (d -params[4]);

	double cos_alpha = cos(CartesianGeometry::angleRadians((n * -1.0),e0));
	double t2 = cos_alpha - params[7];
	double cos_beta = cos(CartesianGeometry::angleRadians(n,e1));

	double denominator = params[5] * sqrt((1-params[7]) * (1-params[8])); 
	double w = 0;
	if(t1 > 0 && t2 > 0) {
		double t3 = cos_beta - params[8];
		//cout << "UUU t3 " << t3 << endl;
		//cout << "UUU cos_beta " << cos_beta << endl;
		if(t3 > 0) {
			w = sqrt(t1 * t2 * t3)/denominator; 
		} else {
			e2 = (CartesianGeometry::buildRadians(pAtoms[2]->getCoor(),pAtoms[3]->getCoor(),pAtoms[4]->getCoor(),params[0],params[1],params[3])) - pAtoms[2]->getCoor();
			cos_beta = cos(CartesianGeometry::angleRadians(n,e2));
			t3 = cos_beta - params[8];
			//cout <<	"e2:" << e2 << endl; 
			//cout << "UUU cos_beta " << cos_beta << endl;
			//cout << "UUU t3 " << t3 << endl;
			if(t3 > 0) {
				w = sqrt(t1 * t2 * t3)/denominator; 
			}
		}
	}        
	
	if(w > 0 && debugFlagOn) {
		cout << " UUU w: " << w << endl; 
		cout << " UUU hydrogen " << *pAtoms[0] << endl;
		cout << " UUU Donor " << *pAtoms[1] << endl;
		cout << " UUU Acceptor " << *pAtoms[2] << endl;
		cout << " UUU Bonded to Acceptor " <<*pAtoms[3] << endl;
		cout << " UUU Angled to Acceptor " <<  *pAtoms[4] << endl;
		cout << " UUU e0: " << e0 + pAtoms[0]->getCoor() << endl; 
		cout << " UUU e1: " << e1 + pAtoms[2]->getCoor() << endl; 
		cout << " UUU e2: " << e2 + pAtoms[2]->getCoor() << endl; 
		cout << " UUU n: " << n + pAtoms[0]->getCoor() << endl; 
		cout << " UUU Energy: " << (scalingFactor * w * pAtoms[0]->getCharge() * pAtoms[2]->getCharge() * params[6]) << endl; 
		cout << "charge 1: " << pAtoms[0]->getCharge() << endl;
		cout << "charge 2: " << pAtoms[2]->getCharge() << endl;
		cout << "B: " << params[6] << endl;
		cout << " UUU d: " << d << endl;
		cout << " UUU t1 " << t1 << endl;
		cout << " UUU cos_alpha " << cos_alpha << endl;
		cout << " UUU t2 " << t2 << endl;
		
	}
	return w;
}	
 
double Scwrl4HBondInteraction::getEnergy() {
	return(scalingFactor * getW() * pAtoms[0]->getCharge() * pAtoms[2]->getCharge() * params[6]);
}

