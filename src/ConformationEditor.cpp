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


#include "ConformationEditor.h"

using namespace MSL;
using namespace std;

ConformationEditor::ConformationEditor() {
	pSys = NULL;
	pIcTable = NULL;
	defineDegreesOfFreedom();
}

ConformationEditor::ConformationEditor(System & _sys) {
	pSys = &_sys;
	pIcTable = &(pSys->getIcTable());
	defineDegreesOfFreedom();
}

ConformationEditor::~ConformationEditor() {
}


bool ConformationEditor::editIC(string _positionId, string _deegreOfFreedom, double _value) {

	// this function checks if 2-4 atom names were given (like "N,CA,CB,CG"), in which
	// case it adds them to the residueId ("A,27,N A,27,CA A,27,CB A,27,CG"). 
	// If only one name was given, it is intepreted as a label (like "chi1"), and processed
	// to convert it to atom names
	if (pSys == NULL || pIcTable == NULL) {
		return false;
	}

	vector<string> tokens = MslTools::tokenizeAndTrim(_deegreOfFreedom, ",");

	if (tokens.size() < 1 || tokens.size() > 4) {
		// syntax error
		return false;
	}

	// parse the _positionId
	string chain;
	int resnum;
	string iCode;
	if (!MslTools::parsePositionId(_positionId, chain, resnum, iCode)) {
		cerr << "WARNING 34828: unknown position " << _positionId << endl;
		return false;
	}


	if (tokens.size() == 1) {
		// it is a label, let's translate it into atom names
		Position * pPos = NULL;
		if (pSys->positionExists(_positionId)) {
			pPos = &(pSys->getLastFoundPosition());
		} else {
			cerr << "WARNING 34833: unknown position " << _positionId << endl;
			return false;
		}

		string resname = pPos->getResidueName();	

		if (degOfFreedomlabels.find(resname) == degOfFreedomlabels.end()) {
			// residue not found in labels
			return false;
		}
		
		map<string, vector<string> >::iterator found = degOfFreedomlabels[resname].find(_deegreOfFreedom);
		if (found == degOfFreedomlabels[resname].end()) {
			// label not found  for residue
			return false;
		}
		tokens = degOfFreedomlabels[resname][_deegreOfFreedom];
		if (tokens.size() < 1 || tokens.size() > 4) {
			// syntax error
			return false;
		}
	}

	// atom names might have + or - signs to indicate that they belong to the 
	// previous or next residue ("-N" or "+C"), let's deal with it
	for (unsigned int i=0; i<tokens.size(); i++) {
		if (tokens[i].size() > 1 && tokens[i].substr(0, 1) == "-") {
			tokens[i] = MslTools::getAtomId(chain, resnum-1, iCode, tokens[i].substr(1, tokens[i].size()-1));
		} else if (tokens[i].substr(0, 1) == "+") {
			tokens[i] = MslTools::getAtomId(chain, resnum+1, iCode, tokens[i].substr(1, tokens[i].size()-1));
		} else {
			tokens[i] = MslTools::getAtomId(chain, resnum, iCode, tokens[i]);
		}
	}


	// find the atoms
	vector<Atom*> pAtoms;
	for (unsigned int i=0; i<tokens.size(); i++) {
		if (pSys->atomExists(tokens[i])) {
			pAtoms.push_back(&pSys->getLastFoundAtom());
		} else {
			// atom not found
			cerr << "WARNING 34837: atom " << tokens[i] << " NOT found" << endl;
			return false;
		}
	}

	return editIC(pAtoms, _value);

}

bool ConformationEditor::editIC(string _deegreOfFreedom, double _value) {

	if (pSys == NULL || pIcTable == NULL) {
		return false;
	}

	// find the atoms
	vector<Atom*> pAtoms;
	vector<string> tokens = MslTools::tokenizeAndTrim(_deegreOfFreedom);
	for (unsigned int i=0; i<tokens.size(); i++) {
		if (pSys->atomExists(tokens[i])) {
			pAtoms.push_back(&pSys->getLastFoundAtom());
		} else {
			// atom not found
			cerr << "WARNING 34837: atom " << tokens[i] << " NOT found" << endl;
			return false;
		}
	}

	return editIC(pAtoms, _value);
	
}

bool ConformationEditor::editIC(vector<Atom*> _pAtoms, double _value) {
	bool out = false;

	set<Atom*> excluded;
	if (_pAtoms.size() < 2 || _pAtoms.size() > 4) {	
		// incorrect number of atoms (ADD A WARNING?)
		return false;
	}
	for (unsigned int i=0; i<_pAtoms.size(); i++) {
		if (i<_pAtoms.size()-1) {
			excluded.insert(_pAtoms[i]);
		}
		if (_pAtoms[i] == NULL) {
			return false;
		}
	}

	if (_pAtoms.size() == 2) {
		out = pIcTable->editBond(_pAtoms[0], _pAtoms[1], _value);
	} else if (_pAtoms.size() == 3) {
		out = pIcTable->editAngle(_pAtoms[0], _pAtoms[1], _pAtoms[2], _value);
	} else if (_pAtoms.size() == 4) {
		out = pIcTable->editDihedral(_pAtoms[0], _pAtoms[1], _pAtoms[2], _pAtoms[3], _value);
	}

	// wipe the coordinates of the last atom in the IC, it will
	// be rebuilt
	_pAtoms.back()->wipeCoordinates();
	buildAtoms.insert(_pAtoms.back());

	// find the atoms bound to the last atom (excluding those of the IC)
	// and wipe their coordinates
	set<Atom*> boundAtoms = _pAtoms.back()->findLinkedAtoms(excluded);
	for (set<Atom*>::iterator k=boundAtoms.begin(); k!=boundAtoms.end(); k++) {
		(*k)->wipeCoordinates();
		buildAtoms.insert(*k);
		//cout << "Wiping coordinates of " << **k << endl;
	}

	return out;

}

void ConformationEditor::defineDegreesOfFreedom() {

	// by default we encode the names for PDB v.2.3

	degOfFreedomlabels.clear();

	// ---------------  ALA  ---------------------- 
	// phi psi
	degOfFreedomlabels["ALA"]["phi"].push_back("-C");
	degOfFreedomlabels["ALA"]["phi"].push_back("N");
	degOfFreedomlabels["ALA"]["phi"].push_back("CA");
	degOfFreedomlabels["ALA"]["phi"].push_back("C");

	degOfFreedomlabels["ALA"]["psi"].push_back("N");
	degOfFreedomlabels["ALA"]["psi"].push_back("CA");
	degOfFreedomlabels["ALA"]["psi"].push_back("C");
	degOfFreedomlabels["ALA"]["psi"].push_back("+N");

	// ---------------  ARG  ---------------------- 
	// phi psi
	degOfFreedomlabels["ARG"]["phi"].push_back("-C");
	degOfFreedomlabels["ARG"]["phi"].push_back("N");
	degOfFreedomlabels["ARG"]["phi"].push_back("CA");
	degOfFreedomlabels["ARG"]["phi"].push_back("C");

	degOfFreedomlabels["ARG"]["psi"].push_back("N");
	degOfFreedomlabels["ARG"]["psi"].push_back("CA");
	degOfFreedomlabels["ARG"]["psi"].push_back("C");
	degOfFreedomlabels["ARG"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["ARG"]["chi1"].push_back("N");
	degOfFreedomlabels["ARG"]["chi1"].push_back("CA");
	degOfFreedomlabels["ARG"]["chi1"].push_back("CB");
	degOfFreedomlabels["ARG"]["chi1"].push_back("CG");

	degOfFreedomlabels["ARG"]["chi2"].push_back("CA");
	degOfFreedomlabels["ARG"]["chi2"].push_back("CB");
	degOfFreedomlabels["ARG"]["chi2"].push_back("CG");
	degOfFreedomlabels["ARG"]["chi2"].push_back("CD");

	degOfFreedomlabels["ARG"]["chi3"].push_back("CB");
	degOfFreedomlabels["ARG"]["chi3"].push_back("CG");
	degOfFreedomlabels["ARG"]["chi3"].push_back("CD");
	degOfFreedomlabels["ARG"]["chi3"].push_back("NE");

	degOfFreedomlabels["ARG"]["chi4"].push_back("CG");
	degOfFreedomlabels["ARG"]["chi4"].push_back("CD");
	degOfFreedomlabels["ARG"]["chi4"].push_back("NE");
	degOfFreedomlabels["ARG"]["chi4"].push_back("CZ");

	// ---------------  ASN  ---------------------- 
	// phi psi
	degOfFreedomlabels["ASN"]["phi"].push_back("-C");
	degOfFreedomlabels["ASN"]["phi"].push_back("N");
	degOfFreedomlabels["ASN"]["phi"].push_back("CA");
	degOfFreedomlabels["ASN"]["phi"].push_back("C");

	degOfFreedomlabels["ASN"]["psi"].push_back("N");
	degOfFreedomlabels["ASN"]["psi"].push_back("CA");
	degOfFreedomlabels["ASN"]["psi"].push_back("C");
	degOfFreedomlabels["ASN"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["ASN"]["chi1"].push_back("N");
	degOfFreedomlabels["ASN"]["chi1"].push_back("CA");
	degOfFreedomlabels["ASN"]["chi1"].push_back("CB");
	degOfFreedomlabels["ASN"]["chi1"].push_back("CG");

	degOfFreedomlabels["ASN"]["chi2"].push_back("CA");
	degOfFreedomlabels["ASN"]["chi2"].push_back("CB");
	degOfFreedomlabels["ASN"]["chi2"].push_back("CG");
	degOfFreedomlabels["ASN"]["chi2"].push_back("OD1");

	// ---------------  ASP  ---------------------- 
	// phi psi
	degOfFreedomlabels["ASP"]["phi"].push_back("-C");
	degOfFreedomlabels["ASP"]["phi"].push_back("N");
	degOfFreedomlabels["ASP"]["phi"].push_back("CA");
	degOfFreedomlabels["ASP"]["phi"].push_back("C");

	degOfFreedomlabels["ASP"]["psi"].push_back("N");
	degOfFreedomlabels["ASP"]["psi"].push_back("CA");
	degOfFreedomlabels["ASP"]["psi"].push_back("C");
	degOfFreedomlabels["ASP"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["ASP"]["chi1"].push_back("N");
	degOfFreedomlabels["ASP"]["chi1"].push_back("CA");
	degOfFreedomlabels["ASP"]["chi1"].push_back("CB");
	degOfFreedomlabels["ASP"]["chi1"].push_back("CG");

	degOfFreedomlabels["ASP"]["chi2"].push_back("CA");
	degOfFreedomlabels["ASP"]["chi2"].push_back("CB");
	degOfFreedomlabels["ASP"]["chi2"].push_back("CG");
	degOfFreedomlabels["ASP"]["chi2"].push_back("OD1");

	// ---------------  CYS  ---------------------- 
	// phi psi
	degOfFreedomlabels["CYS"]["phi"].push_back("-C");
	degOfFreedomlabels["CYS"]["phi"].push_back("N");
	degOfFreedomlabels["CYS"]["phi"].push_back("CA");
	degOfFreedomlabels["CYS"]["phi"].push_back("C");

	degOfFreedomlabels["CYS"]["psi"].push_back("N");
	degOfFreedomlabels["CYS"]["psi"].push_back("CA");
	degOfFreedomlabels["CYS"]["psi"].push_back("C");
	degOfFreedomlabels["CYS"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["CYS"]["chi1"].push_back("N");
	degOfFreedomlabels["CYS"]["chi1"].push_back("CA");
	degOfFreedomlabels["CYS"]["chi1"].push_back("CB");
	degOfFreedomlabels["CYS"]["chi1"].push_back("SG");

	// chis
	degOfFreedomlabels["CYS"]["chi1"].push_back("CA");
	degOfFreedomlabels["CYS"]["chi1"].push_back("CB");
	degOfFreedomlabels["CYS"]["chi1"].push_back("SG");
	degOfFreedomlabels["CYS"]["chi2"].push_back("HG");

	// ---------------  GLN  ---------------------- 
	// phi psi
	degOfFreedomlabels["GLN"]["phi"].push_back("-C");
	degOfFreedomlabels["GLN"]["phi"].push_back("N");
	degOfFreedomlabels["GLN"]["phi"].push_back("CA");
	degOfFreedomlabels["GLN"]["phi"].push_back("C");

	degOfFreedomlabels["GLN"]["psi"].push_back("N");
	degOfFreedomlabels["GLN"]["psi"].push_back("CA");
	degOfFreedomlabels["GLN"]["psi"].push_back("C");
	degOfFreedomlabels["GLN"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["GLN"]["chi1"].push_back("N");
	degOfFreedomlabels["GLN"]["chi1"].push_back("CA");
	degOfFreedomlabels["GLN"]["chi1"].push_back("CB");
	degOfFreedomlabels["GLN"]["chi1"].push_back("CG");

	degOfFreedomlabels["GLN"]["chi2"].push_back("CA");
	degOfFreedomlabels["GLN"]["chi2"].push_back("CB");
	degOfFreedomlabels["GLN"]["chi2"].push_back("CG");
	degOfFreedomlabels["GLN"]["chi2"].push_back("CD");

	degOfFreedomlabels["GLN"]["chi3"].push_back("CB");
	degOfFreedomlabels["GLN"]["chi3"].push_back("CG");
	degOfFreedomlabels["GLN"]["chi3"].push_back("CD");
	degOfFreedomlabels["GLN"]["chi3"].push_back("OE1");

	// ---------------  GLU  ---------------------- 
	// phi psi
	degOfFreedomlabels["GLU"]["phi"].push_back("-C");
	degOfFreedomlabels["GLU"]["phi"].push_back("N");
	degOfFreedomlabels["GLU"]["phi"].push_back("CA");
	degOfFreedomlabels["GLU"]["phi"].push_back("C");

	degOfFreedomlabels["GLU"]["psi"].push_back("N");
	degOfFreedomlabels["GLU"]["psi"].push_back("CA");
	degOfFreedomlabels["GLU"]["psi"].push_back("C");
	degOfFreedomlabels["GLU"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["GLU"]["chi1"].push_back("N");
	degOfFreedomlabels["GLU"]["chi1"].push_back("CA");
	degOfFreedomlabels["GLU"]["chi1"].push_back("CB");
	degOfFreedomlabels["GLU"]["chi1"].push_back("CG");

	degOfFreedomlabels["GLU"]["chi2"].push_back("CA");
	degOfFreedomlabels["GLU"]["chi2"].push_back("CB");
	degOfFreedomlabels["GLU"]["chi2"].push_back("CG");
	degOfFreedomlabels["GLU"]["chi2"].push_back("CD");

	degOfFreedomlabels["GLU"]["chi3"].push_back("CB");
	degOfFreedomlabels["GLU"]["chi3"].push_back("CG");
	degOfFreedomlabels["GLU"]["chi3"].push_back("CD");
	degOfFreedomlabels["GLU"]["chi3"].push_back("OE1");

	// ---------------  GLY  ---------------------- 
	// phi psi
	degOfFreedomlabels["GLY"]["phi"].push_back("-C");
	degOfFreedomlabels["GLY"]["phi"].push_back("N");
	degOfFreedomlabels["GLY"]["phi"].push_back("CA");
	degOfFreedomlabels["GLY"]["phi"].push_back("C");

	degOfFreedomlabels["GLY"]["psi"].push_back("N");
	degOfFreedomlabels["GLY"]["psi"].push_back("CA");
	degOfFreedomlabels["GLY"]["psi"].push_back("C");
	degOfFreedomlabels["GLY"]["psi"].push_back("+N");

	// ---------------  HIS  ---------------------- 
	// phi psi
	degOfFreedomlabels["HIS"]["phi"].push_back("-C");
	degOfFreedomlabels["HIS"]["phi"].push_back("N");
	degOfFreedomlabels["HIS"]["phi"].push_back("CA");
	degOfFreedomlabels["HIS"]["phi"].push_back("C");

	degOfFreedomlabels["HIS"]["psi"].push_back("N");
	degOfFreedomlabels["HIS"]["psi"].push_back("CA");
	degOfFreedomlabels["HIS"]["psi"].push_back("C");
	degOfFreedomlabels["HIS"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["HIS"]["chi1"].push_back("N");
	degOfFreedomlabels["HIS"]["chi1"].push_back("CA");
	degOfFreedomlabels["HIS"]["chi1"].push_back("CB");
	degOfFreedomlabels["HIS"]["chi1"].push_back("CG");

	degOfFreedomlabels["HIS"]["chi2"].push_back("CA");
	degOfFreedomlabels["HIS"]["chi2"].push_back("CB");
	degOfFreedomlabels["HIS"]["chi2"].push_back("CG");
	degOfFreedomlabels["HIS"]["chi2"].push_back("ND1");

	// ---------------  ILE  ---------------------- 
	// phi psi
	degOfFreedomlabels["ILE"]["phi"].push_back("-C");
	degOfFreedomlabels["ILE"]["phi"].push_back("N");
	degOfFreedomlabels["ILE"]["phi"].push_back("CA");
	degOfFreedomlabels["ILE"]["phi"].push_back("C");

	degOfFreedomlabels["ILE"]["psi"].push_back("N");
	degOfFreedomlabels["ILE"]["psi"].push_back("CA");
	degOfFreedomlabels["ILE"]["psi"].push_back("C");
	degOfFreedomlabels["ILE"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["ILE"]["chi1"].push_back("N");
	degOfFreedomlabels["ILE"]["chi1"].push_back("CA");
	degOfFreedomlabels["ILE"]["chi1"].push_back("CB");
	degOfFreedomlabels["ILE"]["chi1"].push_back("CG1");

	degOfFreedomlabels["ILE"]["chi2"].push_back("CA");
	degOfFreedomlabels["ILE"]["chi2"].push_back("CB");
	degOfFreedomlabels["ILE"]["chi2"].push_back("CG1");
	degOfFreedomlabels["ILE"]["chi2"].push_back("CD1");

	// ---------------  LEU  ---------------------- 
	// phi psi
	degOfFreedomlabels["LEU"]["phi"].push_back("-C");
	degOfFreedomlabels["LEU"]["phi"].push_back("N");
	degOfFreedomlabels["LEU"]["phi"].push_back("CA");
	degOfFreedomlabels["LEU"]["phi"].push_back("C");

	degOfFreedomlabels["LEU"]["psi"].push_back("N");
	degOfFreedomlabels["LEU"]["psi"].push_back("CA");
	degOfFreedomlabels["LEU"]["psi"].push_back("C");
	degOfFreedomlabels["LEU"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["LEU"]["chi1"].push_back("N");
	degOfFreedomlabels["LEU"]["chi1"].push_back("CA");
	degOfFreedomlabels["LEU"]["chi1"].push_back("CB");
	degOfFreedomlabels["LEU"]["chi1"].push_back("CG");

	degOfFreedomlabels["LEU"]["chi2"].push_back("CA");
	degOfFreedomlabels["LEU"]["chi2"].push_back("CB");
	degOfFreedomlabels["LEU"]["chi2"].push_back("CG");
	degOfFreedomlabels["LEU"]["chi2"].push_back("CD1");

	// ---------------  LYS  ---------------------- 
	// phi psi
	degOfFreedomlabels["LYS"]["phi"].push_back("-C");
	degOfFreedomlabels["LYS"]["phi"].push_back("N");
	degOfFreedomlabels["LYS"]["phi"].push_back("CA");
	degOfFreedomlabels["LYS"]["phi"].push_back("C");

	degOfFreedomlabels["LYS"]["psi"].push_back("N");
	degOfFreedomlabels["LYS"]["psi"].push_back("CA");
	degOfFreedomlabels["LYS"]["psi"].push_back("C");
	degOfFreedomlabels["LYS"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["LYS"]["chi1"].push_back("N");
	degOfFreedomlabels["LYS"]["chi1"].push_back("CA");
	degOfFreedomlabels["LYS"]["chi1"].push_back("CB");
	degOfFreedomlabels["LYS"]["chi1"].push_back("CG");

	degOfFreedomlabels["LYS"]["chi2"].push_back("CA");
	degOfFreedomlabels["LYS"]["chi2"].push_back("CB");
	degOfFreedomlabels["LYS"]["chi2"].push_back("CG");
	degOfFreedomlabels["LYS"]["chi2"].push_back("CD");

	degOfFreedomlabels["LYS"]["chi3"].push_back("CB");
	degOfFreedomlabels["LYS"]["chi3"].push_back("CG");
	degOfFreedomlabels["LYS"]["chi3"].push_back("CD");
	degOfFreedomlabels["LYS"]["chi3"].push_back("CE");

	degOfFreedomlabels["LYS"]["chi4"].push_back("CG");
	degOfFreedomlabels["LYS"]["chi4"].push_back("CD");
	degOfFreedomlabels["LYS"]["chi4"].push_back("CE");
	degOfFreedomlabels["LYS"]["chi4"].push_back("NZ");

	// ---------------  MET  ---------------------- 
	// phi psi
	degOfFreedomlabels["MET"]["phi"].push_back("-C");
	degOfFreedomlabels["MET"]["phi"].push_back("N");
	degOfFreedomlabels["MET"]["phi"].push_back("CA");
	degOfFreedomlabels["MET"]["phi"].push_back("C");

	degOfFreedomlabels["MET"]["psi"].push_back("N");
	degOfFreedomlabels["MET"]["psi"].push_back("CA");
	degOfFreedomlabels["MET"]["psi"].push_back("C");
	degOfFreedomlabels["MET"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["MET"]["chi1"].push_back("N");
	degOfFreedomlabels["MET"]["chi1"].push_back("CA");
	degOfFreedomlabels["MET"]["chi1"].push_back("CB");
	degOfFreedomlabels["MET"]["chi1"].push_back("CG");

	degOfFreedomlabels["MET"]["chi2"].push_back("CA");
	degOfFreedomlabels["MET"]["chi2"].push_back("CB");
	degOfFreedomlabels["MET"]["chi2"].push_back("CG");
	degOfFreedomlabels["MET"]["chi2"].push_back("SD");

	degOfFreedomlabels["MET"]["chi3"].push_back("CB");
	degOfFreedomlabels["MET"]["chi3"].push_back("CG");
	degOfFreedomlabels["MET"]["chi3"].push_back("SD");
	degOfFreedomlabels["MET"]["chi3"].push_back("CE");

	// ---------------  PHE  ---------------------- 
	// phi psi
	degOfFreedomlabels["PHE"]["phi"].push_back("-C");
	degOfFreedomlabels["PHE"]["phi"].push_back("N");
	degOfFreedomlabels["PHE"]["phi"].push_back("CA");
	degOfFreedomlabels["PHE"]["phi"].push_back("C");

	degOfFreedomlabels["PHE"]["psi"].push_back("N");
	degOfFreedomlabels["PHE"]["psi"].push_back("CA");
	degOfFreedomlabels["PHE"]["psi"].push_back("C");
	degOfFreedomlabels["PHE"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["PHE"]["chi1"].push_back("N");
	degOfFreedomlabels["PHE"]["chi1"].push_back("CA");
	degOfFreedomlabels["PHE"]["chi1"].push_back("CB");
	degOfFreedomlabels["PHE"]["chi1"].push_back("CG");

	degOfFreedomlabels["PHE"]["chi2"].push_back("CA");
	degOfFreedomlabels["PHE"]["chi2"].push_back("CB");
	degOfFreedomlabels["PHE"]["chi2"].push_back("CG");
	degOfFreedomlabels["PHE"]["chi2"].push_back("CD1");

	// ---------------  PRO  ---------------------- 
	// phi psi
	degOfFreedomlabels["PRO"]["phi"].push_back("-C");
	degOfFreedomlabels["PRO"]["phi"].push_back("N");
	degOfFreedomlabels["PRO"]["phi"].push_back("CA");
	degOfFreedomlabels["PRO"]["phi"].push_back("C");

	degOfFreedomlabels["PRO"]["psi"].push_back("N");
	degOfFreedomlabels["PRO"]["psi"].push_back("CA");
	degOfFreedomlabels["PRO"]["psi"].push_back("C");
	degOfFreedomlabels["PRO"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["PRO"]["chi1"].push_back("N");
	degOfFreedomlabels["PRO"]["chi1"].push_back("CA");
	degOfFreedomlabels["PRO"]["chi1"].push_back("CB");
	degOfFreedomlabels["PRO"]["chi1"].push_back("CG");

	degOfFreedomlabels["PRO"]["chi2"].push_back("CA");
	degOfFreedomlabels["PRO"]["chi2"].push_back("CB");
	degOfFreedomlabels["PRO"]["chi2"].push_back("CG");
	degOfFreedomlabels["PRO"]["chi2"].push_back("CD");

	// ---------------  SER  ---------------------- 
	// phi psi
	degOfFreedomlabels["SER"]["phi"].push_back("-C");
	degOfFreedomlabels["SER"]["phi"].push_back("N");
	degOfFreedomlabels["SER"]["phi"].push_back("CA");
	degOfFreedomlabels["SER"]["phi"].push_back("C");

	degOfFreedomlabels["SER"]["psi"].push_back("N");
	degOfFreedomlabels["SER"]["psi"].push_back("CA");
	degOfFreedomlabels["SER"]["psi"].push_back("C");
	degOfFreedomlabels["SER"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["SER"]["chi1"].push_back("N");
	degOfFreedomlabels["SER"]["chi1"].push_back("CA");
	degOfFreedomlabels["SER"]["chi1"].push_back("CB");
	degOfFreedomlabels["SER"]["chi1"].push_back("OG");

	degOfFreedomlabels["SER"]["chi2"].push_back("CA");
	degOfFreedomlabels["SER"]["chi2"].push_back("CB");
	degOfFreedomlabels["SER"]["chi2"].push_back("OG");
	degOfFreedomlabels["SER"]["chi2"].push_back("HG");

	// ---------------  THR  ---------------------- 
	// phi psi
	degOfFreedomlabels["THR"]["phi"].push_back("-C");
	degOfFreedomlabels["THR"]["phi"].push_back("N");
	degOfFreedomlabels["THR"]["phi"].push_back("CA");
	degOfFreedomlabels["THR"]["phi"].push_back("C");

	degOfFreedomlabels["THR"]["psi"].push_back("N");
	degOfFreedomlabels["THR"]["psi"].push_back("CA");
	degOfFreedomlabels["THR"]["psi"].push_back("C");
	degOfFreedomlabels["THR"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["THR"]["chi1"].push_back("N");
	degOfFreedomlabels["THR"]["chi1"].push_back("CA");
	degOfFreedomlabels["THR"]["chi1"].push_back("CB");
	degOfFreedomlabels["THR"]["chi1"].push_back("OG1");

	degOfFreedomlabels["THR"]["chi2"].push_back("CA");
	degOfFreedomlabels["THR"]["chi2"].push_back("CB");
	degOfFreedomlabels["THR"]["chi2"].push_back("OG1");
	degOfFreedomlabels["THR"]["chi2"].push_back("HG1");

	// ---------------  TRP  ---------------------- 
	// phi psi
	degOfFreedomlabels["TRP"]["phi"].push_back("-C");
	degOfFreedomlabels["TRP"]["phi"].push_back("N");
	degOfFreedomlabels["TRP"]["phi"].push_back("CA");
	degOfFreedomlabels["TRP"]["phi"].push_back("C");

	degOfFreedomlabels["TRP"]["psi"].push_back("N");
	degOfFreedomlabels["TRP"]["psi"].push_back("CA");
	degOfFreedomlabels["TRP"]["psi"].push_back("C");
	degOfFreedomlabels["TRP"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["TRP"]["chi1"].push_back("N");
	degOfFreedomlabels["TRP"]["chi1"].push_back("CA");
	degOfFreedomlabels["TRP"]["chi1"].push_back("CB");
	degOfFreedomlabels["TRP"]["chi1"].push_back("CG");

	degOfFreedomlabels["TRP"]["chi2"].push_back("CA");
	degOfFreedomlabels["TRP"]["chi2"].push_back("CB");
	degOfFreedomlabels["TRP"]["chi2"].push_back("CG");
	degOfFreedomlabels["TRP"]["chi2"].push_back("CD1");

	// ---------------  TYR  ---------------------- 
	// phi psi
	degOfFreedomlabels["TYR"]["phi"].push_back("-C");
	degOfFreedomlabels["TYR"]["phi"].push_back("N");
	degOfFreedomlabels["TYR"]["phi"].push_back("CA");
	degOfFreedomlabels["TYR"]["phi"].push_back("C");

	degOfFreedomlabels["TYR"]["psi"].push_back("N");
	degOfFreedomlabels["TYR"]["psi"].push_back("CA");
	degOfFreedomlabels["TYR"]["psi"].push_back("C");
	degOfFreedomlabels["TYR"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["TYR"]["chi1"].push_back("N");
	degOfFreedomlabels["TYR"]["chi1"].push_back("CA");
	degOfFreedomlabels["TYR"]["chi1"].push_back("CB");
	degOfFreedomlabels["TYR"]["chi1"].push_back("CG");

	degOfFreedomlabels["TYR"]["chi2"].push_back("CA");
	degOfFreedomlabels["TYR"]["chi2"].push_back("CB");
	degOfFreedomlabels["TYR"]["chi2"].push_back("CG");
	degOfFreedomlabels["TYR"]["chi2"].push_back("CD1");

	degOfFreedomlabels["TYR"]["chi3"].push_back("CE1");
	degOfFreedomlabels["TYR"]["chi3"].push_back("CZ");
	degOfFreedomlabels["TYR"]["chi3"].push_back("OH");
	degOfFreedomlabels["TYR"]["chi3"].push_back("HH");

	// ---------------  VAL  ---------------------- 
	// phi psi
	degOfFreedomlabels["VAL"]["phi"].push_back("-C");
	degOfFreedomlabels["VAL"]["phi"].push_back("N");
	degOfFreedomlabels["VAL"]["phi"].push_back("CA");
	degOfFreedomlabels["VAL"]["phi"].push_back("C");

	degOfFreedomlabels["VAL"]["psi"].push_back("N");
	degOfFreedomlabels["VAL"]["psi"].push_back("CA");
	degOfFreedomlabels["VAL"]["psi"].push_back("C");
	degOfFreedomlabels["VAL"]["psi"].push_back("+N");

	// chis
	degOfFreedomlabels["VAL"]["chi1"].push_back("N");
	degOfFreedomlabels["VAL"]["chi1"].push_back("CA");
	degOfFreedomlabels["VAL"]["chi1"].push_back("CB");
	degOfFreedomlabels["VAL"]["chi1"].push_back("CG1");

}


bool ConformationEditor::readDefinitionFile(string _defiFile) { 
	/**************************************************
	 *  This overwrites the definitions from file.
	 *
	 *    File format:
	 *
	 * ALA phi -C N CA C 
	 * ALA psi N CA C +N 
	 * ARG phi -C N CA C 
	 * ARG psi N CA C +N 
	 * ARG chi1 N CA CB CG 
	 * ARG chi2 CA CB CG CD 
	 * ARG chi3 CB CG CD NE 
	 * ARG chi4 CG CD NE CZ 
	 * ASN phi -C N CA C 
	 * ASN psi N CA C +N 
	 * ASN chi1 N CA CB CG 
	 * ASN chi2 CA CB CG OD1 
	 * ...
	 **************************************************/

	degOfFreedomlabels.clear();
	//cout << "UUU cleared" << endl;
	DegreeOfFreedomReader reader;
	if (!reader.read(_defiFile)) {
		return false;
	}

	degOfFreedomlabels = reader.getDegreesOfFreedom();
	vector<string> v = reader.getSingleDegreeOfFreedom("LEU", "chi2");
	//cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
/*
	Reader reader;
	if (!reader.open(_defiFile)) {
		cerr << "WARNING 34842: cannot read definition file " << _defiFile << endl;
		return false;
	}
	vector<string> lines = reader.getAllLines();
	reader.close();

	degOfFreedomlabels.clear();

	for (unsigned int i=0; i<lines.size(); i++) {
		string line = MslTools::uncomment(lines[i]);
		vector<string> tokens = MslTools::tokenizeAndTrim(line, " ");
		if (tokens.size() == 0) {
			continue;
		}
		//for (unsigned int j=0; j<tokens.size(); j++) {
		//	cout << tokens[j] << " ";
		//}
		//cout << endl;
		if (tokens.size() < 4) {
			// not enough elements
			cerr << "WARNING 34847: syntax error in definition line " << line << " in file " << _defiFile << endl;
			continue;
		}
		for (unsigned int j=2; j<tokens.size(); j++) {
			// 0: residue   1: label   2-5: atoms
			degOfFreedomlabels[tokens[0]][tokens[1]].push_back(tokens[j]);
		}
	}
*/

	return true;

}

