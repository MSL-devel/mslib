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

#include "PDBFormat.h"



PDBFormat::AtomData PDBFormat::parseAtomLine(const string &_pdbAtomLine){
	AtomData atom;

	int lineLength = _pdbAtomLine.size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_RECORD_NAME)    strcpy(atom.D_RECORD_NAME,    MslTools::trim(_pdbAtomLine.substr(S_RECORD_NAME, L_RECORD_NAME)).c_str());
		if (lineLength >= E_ATOM_NAME)      strcpy(atom.D_ATOM_NAME,      MslTools::trim(_pdbAtomLine.substr(S_ATOM_NAME, L_ATOM_NAME)).c_str());
		if (lineLength >= E_ALT_LOC)        strcpy(atom.D_ALT_LOC,        MslTools::trim(_pdbAtomLine.substr(S_ALT_LOC,L_ALT_LOC)).c_str());
		if (lineLength >= E_RES_NAME)       strcpy(atom.D_RES_NAME,       MslTools::trim(_pdbAtomLine.substr(S_RES_NAME, L_RES_NAME)).c_str());
		if (lineLength >= E_CHAIN_ID)       strcpy(atom.D_CHAIN_ID,       MslTools::trim(_pdbAtomLine.substr(S_CHAIN_ID, L_CHAIN_ID)).c_str());
		if (lineLength >= E_I_CODE)         strcpy(atom.D_I_CODE,         MslTools::trim(_pdbAtomLine.substr(S_I_CODE, L_I_CODE)).c_str());
		if (lineLength >= E_SEG_ID)         strcpy(atom.D_SEG_ID,         MslTools::trim(_pdbAtomLine.substr(S_SEG_ID, L_SEG_ID)).c_str());
		if (lineLength >= E_ELEMENT_SYMBOL) strcpy(atom.D_ELEMENT_SYMBOL, MslTools::trim(_pdbAtomLine.substr(S_ELEMENT_SYMBOL, L_ELEMENT_SYMBOL)).c_str());

		if (lineLength >= E_SERIAL)    atom.D_SERIAL         = MslTools::toInt(MslTools::trim(_pdbAtomLine.substr(S_SERIAL,L_SERIAL)), "Atom number(Serial) problem");
		if (lineLength >= E_RES_SEQ)   atom.D_RES_SEQ        = MslTools::toInt(MslTools::trim(_pdbAtomLine.substr(S_RES_SEQ,L_RES_SEQ)), "Residue Number problem");

		if (lineLength >= E_X)         atom.D_X              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_X,L_X)), "X-coord problem");
		if (lineLength >= E_Y)         atom.D_Y              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_Y,L_Y)), "Y-coord problem");
		if (lineLength >= E_Z)         atom.D_Z              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_Z,L_Z)), "Z-coord problem");
		if (lineLength >= E_OCCUP)     atom.D_OCCUP          = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_OCCUP, L_OCCUP)), "Occupation problem");
		if (lineLength >= E_TEMP_FACT) atom.D_TEMP_FACT      = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_TEMP_FACT, L_TEMP_FACT)), "Temp. Factor problem");

		// Charge is  a problem
		if (lineLength >= E_CHARGE && MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)) != "")    atom.D_CHARGE         = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)), "CHARGE problem");
			


	} catch(exception &e){
		cerr << "ERROR 34918 PDBFormat parseAtomLine "<<e.what()<<endl;
		exit(34918);
	}

	return atom;
}


PDBFormat::AtomData PDBFormat::createAtomData(string _resName, Real &_x, Real &_y, Real &_z, string _element){
	Atom a(_resName, _x,_y,_z,_element);
	return createAtomData(a);
		
}
PDBFormat::AtomData PDBFormat::createAtomData(const Atom &_at){
	PDBFormat::AtomData atom;

	atom.D_X = _at.getX();
	atom.D_Y = _at.getY();
	atom.D_Z = _at.getZ();

	// TODO: Need to check getName size < L_ATOM_NAME ???, does it truncate automatically
	strncpy(atom.D_ATOM_NAME, _at.getName().c_str(), PDBFormat::L_ATOM_NAME);


	strncpy(atom.D_RECORD_NAME , "ATOM  ",PDBFormat::L_RECORD_NAME);
	strncpy(atom.D_RES_NAME, _at.getResidueName().c_str(), PDBFormat::L_RES_NAME);
		// 
	strncpy(atom.D_CHAIN_ID, _at.getChainId().c_str(), PDBFormat::L_CHAIN_ID);
	strncpy(atom.D_I_CODE, _at.getResidueIcode().c_str(), PDBFormat::L_I_CODE);
	strncpy(atom.D_ELEMENT_SYMBOL, _at.getElement().c_str(), PDBFormat::L_ELEMENT_SYMBOL);
	strncpy(atom.D_SEG_ID,_at.getSegID().c_str(),PDBFormat::L_SEG_ID);

	atom.D_SERIAL  = 1;
	atom.D_RES_SEQ = _at.getResidueNumber();
	atom.D_TEMP_FACT = _at.getTempFactor();
	atom.D_OCCUP   = 1.0;

	return atom;

}

string PDBFormat::createAtomLine(const PDBFormat::AtomData &ad){

	// the first column of the atom name is blank unless it is a number (used by some hydrogen atoms)
	bool startWithNum = false;
	string atomName = ad.D_ATOM_NAME;
	if (atomName.length() > 0) {
		int asciiCode = atomName[0];
		if (asciiCode >= 48 && asciiCode <= 57) {
			startWithNum = true;
		}
	}

	string atomNamePdb = atomName;
	if (!startWithNum && atomName.length() < 4) {
		atomNamePdb = " " + atomName;
	}

	/*
		 1         2         3         4         5         6         7         8
	----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
	<----><--->x<-->|<->x|<-->|xxx<------><------><------><----><---->xxxxxx<--><><>
	ATOM    105  HD2 HIS A   9       8.129   4.666  -7.324  1.00  0.00           H  
	ATOM    106  HE1 HIS A   9      11.912   3.392  -8.763  1.00  0.00           H  
	ATOM    107  HE2 HIS A   9      10.415   5.367  -8.312  1.00  0.00           H  
	ATOM    108  N   VAL A  10       5.773  -0.015  -6.169  1.00  0.00           N  
	ATOM    109  CA  VAL A  10       4.784  -0.600  -5.271  1.00  0.00           C
	*/
	char c [1000];
	sprintf(c, "%6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s", ad.D_RECORD_NAME, ad.D_SERIAL, atomNamePdb.c_str(), ad.D_ALT_LOC,ad.D_RES_NAME,ad.D_CHAIN_ID,ad.D_RES_SEQ, ad.D_I_CODE, ad.D_X, ad.D_Y, ad.D_Z,ad.D_OCCUP,ad.D_TEMP_FACT,ad.D_SEG_ID,ad.D_ELEMENT_SYMBOL);
	
	return (string)c;

}

string PDBFormat::createTerLine(const PDBFormat::AtomData &ad){

	// the first column of the atom name is blank unless it is a number (used by some hydrogen atoms)
	bool startWithNum = false;
	string atomName = ad.D_ATOM_NAME;
	if (atomName.length() > 0) {
		int asciiCode = atomName[0];
		if (asciiCode >= 48 && asciiCode <= 57) {
			startWithNum = true;
		}
	}

	string atomNamePdb = atomName;
	if (!startWithNum && atomName.length() < 4) {
		atomNamePdb = " " + atomName;
	}

	/*
		 1         2         3         4         5         6         7         8
	----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
	<----><--->x<-->|<->x|<-->|xxx<------><------><------><----><---->xxxxxx<--><><>
	ATOM    105  HD2 HIS A   9       8.129   4.666  -7.324  1.00  0.00           H  
	ATOM    106  HE1 HIS A   9      11.912   3.392  -8.763  1.00  0.00           H  
	ATOM    107  HE2 HIS A   9      10.415   5.367  -8.312  1.00  0.00           H  
	ATOM    108  N   VAL A  10       5.773  -0.015  -6.169  1.00  0.00           N  
	ATOM    109  CA  VAL A  10       4.784  -0.600  -5.271  1.00  0.00           C
	*/
	char c [1000];
	sprintf(c, "%6s%5d     %1s%-3s %1s%4d%1s                                                   ", ad.D_RECORD_NAME, ad.D_SERIAL, ad.D_ALT_LOC,ad.D_RES_NAME,ad.D_CHAIN_ID,ad.D_RES_SEQ, ad.D_I_CODE);
	
	return (string)c;

}
