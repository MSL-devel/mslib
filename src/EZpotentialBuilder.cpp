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


#include "EZpotentialBuilder.h"

using namespace MSL;
using namespace std;


EZpotentialBuilder::EZpotentialBuilder() {
	setup();
}

EZpotentialBuilder::EZpotentialBuilder(System & _system) {
	setup();
	pSystem = &_system;
}

EZpotentialBuilder::EZpotentialBuilder( EZpotentialBuilder & _sysBuild) {
	setup();
	copy(_sysBuild);
}

EZpotentialBuilder::~EZpotentialBuilder() {
	deletePointers();
}

void EZpotentialBuilder::operator=( EZpotentialBuilder & _sysBuild) {
	copy(_sysBuild);
}


void EZpotentialBuilder::setup() {
	pSystem = NULL;
	setParams();
	useCB_flag = false;
	addTermini_flag = false;
}

void EZpotentialBuilder::copy( EZpotentialBuilder & _sysBuild) {
	pSystem = _sysBuild.pSystem;
}

void EZpotentialBuilder::deletePointers() {
}

bool EZpotentialBuilder::buildInteractions() {

	EnergySet* ESet = pSystem->getEnergySet();
	// delete all existing scwrl4HBondinteractions
	ESet->eraseTerm("EZ_POTENTIAL");

	AtomPointerVector atoms = pSystem->getAllAtomPointers();

	for (unsigned int i=0; i<atoms.size(); i++) {
		string resName = atoms[i]->getResidueName();
		if (parameters.find(resName) != parameters.end()) {
			string atomName = atoms[i]->getName();
			if (parameters[resName].find(atomName) != parameters[resName].end()) {
				EZpotentialInteraction * interaction = new EZpotentialInteraction(*atoms[i], parameters[resName][atomName].params, parameters[resName][atomName].sigmoidalFunc);
				ESet->addInteraction(interaction);
			}
			if (addTermini_flag && parameters["_NTER_"].find(atomName) != parameters["_NTER_"].end() && atoms[i]->isPositionNterminal()) {
				EZpotentialInteraction * interaction = new EZpotentialInteraction(*atoms[i], parameters["_NTER_"][atomName].params, parameters["_NTER_"][atomName].sigmoidalFunc);
				ESet->addInteraction(interaction);
			}
			if (addTermini_flag && parameters["_CTER_"].find(atomName) != parameters["_CTER_"].end() && atoms[i]->isPositionCterminal()) {
				EZpotentialInteraction * interaction = new EZpotentialInteraction(*atoms[i], parameters["_CTER_"][atomName].params, parameters["_CTER_"][atomName].sigmoidalFunc);
				ESet->addInteraction(interaction);
			}
		}
	}

	return true;
}


void EZpotentialBuilder::setParams() {

	string atom = "";

	// ALA
	atom = "CB";
	parameters["ALA"][atom].sigmoidalFunc = true;
	parameters["ALA"][atom].params = vector<double>(3, 0.0);
	parameters["ALA"][atom].params[0] = -0.28750;
	parameters["ALA"][atom].params[1] = 10.22;
	parameters["ALA"][atom].params[2] = 4.67200;
	parameters["ALA"][atom].Nterminal = false;
	parameters["ALA"][atom].Cterminal = false;

	// ASP
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CG";
	}
	parameters["ASP"][atom].sigmoidalFunc = true;
	parameters["ASP"][atom].params = vector<double>(3, 0.0);
	parameters["ASP"][atom].params[0] = 1.19280;
	parameters["ASP"][atom].params[1] = 14.25;
	parameters["ASP"][atom].params[2] = 8.98150;
	parameters["ASP"][atom].Nterminal = false;
	parameters["ASP"][atom].Cterminal = false;

	// GLU
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CD";
	}
	parameters["GLU"][atom].sigmoidalFunc = true;
	parameters["GLU"][atom].params = vector<double>(3, 0.0);
	parameters["GLU"][atom].params[0] = 1.29640;
	parameters["GLU"][atom].params[1] = 14.66;
	parameters["GLU"][atom].params[2] = 4.16260;
	parameters["GLU"][atom].Nterminal = false;
	parameters["GLU"][atom].Cterminal = false;

	// PHE
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CZ";
	}
	parameters["PHE"][atom].sigmoidalFunc = true;
	parameters["PHE"][atom].params = vector<double>(3, 0.0);
	parameters["PHE"][atom].params[0] = -0.80007;
	parameters["PHE"][atom].params[1] = 19.66600;
	parameters["PHE"][atom].params[2] = 7.12320;
	parameters["PHE"][atom].Nterminal = false;
	parameters["PHE"][atom].Cterminal = false;


	// GLY
	atom = "CA";
	parameters["GLY"][atom].sigmoidalFunc = true;
	parameters["GLY"][atom].params = vector<double>(3, 0.0);
	parameters["GLY"][atom].params[0] = -0.01;
	parameters["GLY"][atom].params[1] = 13.86;
	parameters["GLY"][atom].params[2] = 6;
	parameters["GLY"][atom].Nterminal = false;
	parameters["GLY"][atom].Cterminal = false;

	// HIS
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CE1";
	}
	parameters["HIS"][atom].sigmoidalFunc = true;
	parameters["HIS"][atom].params = vector<double>(3, 0.0);
	parameters["HIS"][atom].params[0] = 0.75243;
	parameters["HIS"][atom].params[1] = 12.26;
	parameters["HIS"][atom].params[2] = 2.77;
	parameters["HIS"][atom].Nterminal = false;
	parameters["HIS"][atom].Cterminal = false;

	// HSD (CHRMM 22)
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CE1";
	}
	parameters["HSD"][atom].sigmoidalFunc = true;
	parameters["HSD"][atom].params = vector<double>(3, 0.0);
	parameters["HSD"][atom].params[0] = 0.75243;
	parameters["HSD"][atom].params[1] = 12.26;
	parameters["HSD"][atom].params[2] = 2.77;
	parameters["HSD"][atom].Nterminal = false;
	parameters["HSD"][atom].Cterminal = false;

	// HSE (CHRMM 22)
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CE1";
	}
	parameters["HSE"][atom].sigmoidalFunc = true;
	parameters["HSE"][atom].params = vector<double>(3, 0.0);
	parameters["HSE"][atom].params[0] = 0.75243;
	parameters["HSE"][atom].params[1] = 12.26;
	parameters["HSE"][atom].params[2] = 2.77;
	parameters["HSE"][atom].Nterminal = false;
	parameters["HSE"][atom].Cterminal = false;

	// HSP (CHRMM 22)
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CE1";
	}
	parameters["HSP"][atom].sigmoidalFunc = true;
	parameters["HSP"][atom].params = vector<double>(3, 0.0);
	parameters["HSP"][atom].params[0] = 0.75243;
	parameters["HSP"][atom].params[1] = 12.26;
	parameters["HSP"][atom].params[2] = 2.77;
	parameters["HSP"][atom].Nterminal = false;
	parameters["HSP"][atom].Cterminal = false;

	// ILE
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CG1";
	}
	parameters["ILE"][atom].sigmoidalFunc = true;
	parameters["ILE"][atom].params = vector<double>(3, 0.0);
	parameters["ILE"][atom].params[0] = -0.55854;
	parameters["ILE"][atom].params[1] = 14.34;
	parameters["ILE"][atom].params[2] = 10.6860;
	parameters["ILE"][atom].Nterminal = false;
	parameters["ILE"][atom].Cterminal = false;

	// LYS
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CE";
	}
	parameters["LYS"][atom].sigmoidalFunc = true;
	parameters["LYS"][atom].params = vector<double>(3, 0.0);
	parameters["LYS"][atom].params[0]   = 1.66430;
	parameters["LYS"][atom].params[1] = 11.11;
	parameters["LYS"][atom].params[2] = 2.08800;
	parameters["LYS"][atom].Nterminal = false;
	parameters["LYS"][atom].Cterminal = false;

	// LEU
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CG";
	}
	parameters["LEU"][atom].sigmoidalFunc = true;
	parameters["LEU"][atom].params = vector<double>(3, 0.0);
	parameters["LEU"][atom].params[0] = -0.64267;
	parameters["LEU"][atom].params[1] = 17.34;
	parameters["LEU"][atom].params[2] = 8.61;
	parameters["LEU"][atom].Nterminal = false;
	parameters["LEU"][atom].Cterminal = false;

	// MET
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "SD";
	}
	parameters["MET"][atom].sigmoidalFunc = true;
	parameters["MET"][atom].params = vector<double>(3, 0.0);
	parameters["MET"][atom].params[0] = -0.28;
	parameters["MET"][atom].params[1] = 18.04;
	parameters["MET"][atom].params[2] = 7.13;
	parameters["MET"][atom].Nterminal = false;
	parameters["MET"][atom].Cterminal = false;

	// ASN
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CG";
	}
	parameters["ASN"][atom].sigmoidalFunc = true;
	parameters["ASN"][atom].params = vector<double>(3, 0.0);
	parameters["ASN"][atom].params[0] = 0.89;
	parameters["ASN"][atom].params[1] = 12.78;
	parameters["ASN"][atom].params[2] = 6.28;
	parameters["ASN"][atom].Nterminal = false;
	parameters["ASN"][atom].Cterminal = false;

	// PRO
	atom = "CB";
	parameters["PRO"][atom].sigmoidalFunc = true;
	parameters["PRO"][atom].params = vector<double>(3, 0.0);
	parameters["PRO"][atom].params[0] = 0.82570;
	parameters["PRO"][atom].params[1] = 18.09;
	parameters["PRO"][atom].params[2] = 3.53140;
	parameters["PRO"][atom].Nterminal = false;
	parameters["PRO"][atom].Cterminal = false;

	// GLN
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CD";
	}
	parameters["GLN"][atom].sigmoidalFunc = true;
	parameters["GLN"][atom].params = vector<double>(3, 0.0);
	parameters["GLN"][atom].params[0]   = 1.20550;
	parameters["GLN"][atom].params[1] = 10.46;
	parameters["GLN"][atom].params[2] = 2.59380;
	parameters["GLN"][atom].Nterminal = false;
	parameters["GLN"][atom].Cterminal = false;

	// ARG
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CZ";
	}
	parameters["ARG"][atom].sigmoidalFunc = true;
	parameters["ARG"][atom].params = vector<double>(3, 0.0);
	parameters["ARG"][atom].params[0] = 1.54550;
	parameters["ARG"][atom].params[1] = 9.34;
	parameters["ARG"][atom].params[2] = 4.68;
	parameters["ARG"][atom].Nterminal = false;
	parameters["ARG"][atom].Cterminal = false;

	// SER
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "OG";
	}
	parameters["SER"][atom].sigmoidalFunc = true;
	parameters["SER"][atom].params = vector<double>(3, 0.0);
	parameters["SER"][atom].params[0] = 0.1;
	parameters["SER"][atom].params[1] = 13.86;
	parameters["SER"][atom].params[2] = 6.0;
	parameters["SER"][atom].Nterminal = false;
	parameters["SER"][atom].Cterminal = false;

	// THR
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "OG1";
	}
	parameters["THR"][atom].sigmoidalFunc = true;
	parameters["THR"][atom].params = vector<double>(3, 0.0);
	parameters["THR"][atom].params[0] = 0.01;
	parameters["THR"][atom].params[1] = 13.86;
	parameters["THR"][atom].params[2] = 6.0;
	parameters["THR"][atom].Nterminal = false;
	parameters["THR"][atom].Cterminal = false;

	// VAL";
	atom = "CB";
	parameters["VAL"][atom].sigmoidalFunc = true;
	parameters["VAL"][atom].params = vector<double>(3, 0.0);
	parameters["VAL"][atom].params[0] = -0.47227;
	parameters["VAL"][atom].params[1] = 11.35;
	parameters["VAL"][atom].params[2] = 4.97;
	parameters["VAL"][atom].Nterminal = false;
	parameters["VAL"][atom].Cterminal = false;

	// TRP
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "CZ3";
	}
	parameters["TRP"][atom].sigmoidalFunc = false;
	parameters["TRP"][atom].params = vector<double>(3, 0.0);
	parameters["TRP"][atom].params[0] = 7.20;
	parameters["TRP"][atom].params[1] = 11.65;
	parameters["TRP"][atom].params[2] = -0.85;
	parameters["TRP"][atom].Nterminal = false;
	parameters["TRP"][atom].Cterminal = false;

	// TYR
	if (useCB_flag) {
		atom = "CB";
	} else {
		atom = "OH";
	}
	parameters["TYR"][atom].sigmoidalFunc = false;
	parameters["TYR"][atom].params = vector<double>(3, 0.0);
	parameters["TYR"][atom].params[0] = 6.20;
	parameters["TYR"][atom].params[1] = 13.04;
	parameters["TYR"][atom].params[2] = -0.42;
	parameters["TYR"][atom].Nterminal = false;
	parameters["TYR"][atom].Cterminal = false;

	// _ANY_
	atom = "N";
	parameters["_NTER_"][atom].sigmoidalFunc = true;
	parameters["_NTER_"][atom].params = vector<double>(3, 0.0);
	parameters["_NTER_"][atom].params[0] = 1.5742;
	parameters["_NTER_"][atom].params[1] = 11.4590;
	parameters["_NTER_"][atom].params[2] = 7.1879;
	parameters["_NTER_"][atom].Nterminal = true;
	parameters["_NTER_"][atom].Cterminal = false;

	// _ANY_
	atom = "C";
	parameters["_CTER_"][atom].sigmoidalFunc = true;
	parameters["_CTER_"][atom].params = vector<double>(3, 0.0);
	parameters["_CTER_"][atom].params[0] = 0.6095;
	parameters["_CTER_"][atom].params[1] = 13.1330;
	parameters["_CTER_"][atom].params[2] = 6.8147;
	parameters["_CTER_"][atom].Nterminal = false;
	parameters["_CTER_"][atom].Cterminal = true;

}



