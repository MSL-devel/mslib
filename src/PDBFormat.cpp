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

#include "PDBFormat.h"

using namespace MSL;
using namespace std;



PDBFormat::CrystData PDBFormat::parseCrystLine(const string &_pdbCrystLine){
		
	CrystData cryst;

	int lineLength = MslTools::trim(_pdbCrystLine).size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_CRYSTRECORD)         strcpy(cryst.D_CRYSTRECORD,    MslTools::trim(_pdbCrystLine.substr(S_CRYSTRECORD, L_CRYSTRECORD)).c_str());
		if (lineLength >= E_CRYSTSPACEGROUP)     strcpy(cryst.D_CRYSTSPACEGROUP,    MslTools::trim(_pdbCrystLine.substr(S_CRYSTSPACEGROUP, L_CRYSTSPACEGROUP)).c_str());

		if (lineLength >= E_CRYSTINDEX)   cryst.D_CRYSTINDEX       = MslTools::toInt(MslTools::trim(_pdbCrystLine.substr(S_CRYSTINDEX,L_CRYSTINDEX)), "CRYST INDEX problem");
		/* Causing problems, some PDBs not to spec? */
		//if (lineLength >= E_CRYSTZ)       cryst.D_CRYSTZ           = MslTools::toInt(MslTools::trim(_pdbCrystLine.substr(S_CRYSTZ,L_CRYSTZ)), "CRYST Z problem");

		if (lineLength >= E_CRYSTA)       cryst.D_CRYSTA               = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTA,L_CRYSTA)), "CRYST A problem");
		if (lineLength >= E_CRYSTB)       cryst.D_CRYSTB               = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTB,L_CRYSTB)), "CRYST B problem");
		if (lineLength >= E_CRYSTC)       cryst.D_CRYSTC               = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTC,L_CRYSTC)), "CRYST C problem");
		if (lineLength >= E_CRYSTALPHA)   cryst.D_CRYSTALPHA           = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTALPHA,L_CRYSTALPHA)), "CRYST ALPHA problem");
		if (lineLength >= E_CRYSTBETA)    cryst.D_CRYSTBETA            = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTBETA,L_CRYSTBETA)), "CRYST BETA problem");
		if (lineLength >= E_CRYSTGAMMA)    cryst.D_CRYSTGAMMA           = MslTools::toDouble(MslTools::trim(_pdbCrystLine.substr(S_CRYSTGAMMA,L_CRYSTGAMMA)), "CRYST GAMMA problem");


	} catch(exception &e){
		cerr << "ERROR 34918 PDBFormat parseCrystLine "<<e.what()<<endl;
		exit(34918);
	}

	return cryst;
}

PDBFormat::ScaleData PDBFormat::parseScaleLine(const string &_pdbScaleLine){
	
	ScaleData scale;

	int lineLength = MslTools::trim(_pdbScaleLine).size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_SCALERECORD)     strcpy(scale.D_SCALERECORD,    MslTools::trim(_pdbScaleLine.substr(S_SCALERECORD, L_SCALERECORD)).c_str());

		if (lineLength >= E_SCALELINE)   scale.D_SCALELINE       = MslTools::toInt(MslTools::trim(_pdbScaleLine.substr(S_SCALELINE,L_SCALELINE)), "SCALE LINE problem");
		if (lineLength >= E_SCALEX)       scale.D_SCALEX           = MslTools::toDouble(MslTools::trim(_pdbScaleLine.substr(S_SCALEX,L_SCALEX)), "SCALE X problem");
		if (lineLength >= E_SCALEY)       scale.D_SCALEY           = MslTools::toDouble(MslTools::trim(_pdbScaleLine.substr(S_SCALEY,L_SCALEY)), "SCALE Y problem");
		if (lineLength >= E_SCALEZ)       scale.D_SCALEZ           = MslTools::toDouble(MslTools::trim(_pdbScaleLine.substr(S_SCALEZ,L_SCALEZ)), "SCALE Z problem");
		if (lineLength >= E_SCALETRANS)   scale.D_SCALETRANS       = MslTools::toDouble(MslTools::trim(_pdbScaleLine.substr(S_SCALETRANS,L_SCALETRANS)), "SCALE trans problem");


	} catch(exception &e){
		cerr << "ERROR 34918 PDBFormat parseScaleLine "<<e.what()<<endl;
		exit(34918);
	}

	return scale;
}

PDBFormat::ModelData PDBFormat::parseModelLine(const string &_pdbModelLine){
	
	ModelData model;

	int lineLength = MslTools::trim(_pdbModelLine).size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_MODEL_RECORD) {
			if (_pdbModelLine.substr(S_MODEL_RECORD, L_MODEL_RECORD) == "MODEL ") {
				model.D_ENDMODEL_FLAG = false;
				if(lineLength >= E_MODEL_NUMBER) {
					model.D_MODEL_NUMBER = MslTools::toInt(MslTools::trim(_pdbModelLine.substr(S_MODEL_NUMBER,L_MODEL_NUMBER)), "Model line number problem");
				}
			}
		}

	} catch(exception &e){
		cerr << "ERROR 34922 PDBFormat parseModelLine "<<e.what()<<endl;
		exit(34922);
	}

	return model;
}

PDBFormat::SymData PDBFormat::parseSymLine(const string &_pdbSymLine){

	SymData sym;

	int lineLength = MslTools::trim(_pdbSymLine).size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_SYMMRECORD)     strcpy(sym.D_SYMMRECORD,    MslTools::trim(_pdbSymLine.substr(S_SYMMRECORD, L_SYMMRECORD)).c_str());
		
		if (lineLength >= E_SYMMLINE)    sym.D_SYMMLINE        = MslTools::toInt(MslTools::trim(_pdbSymLine.substr(S_SYMMLINE,L_SYMMLINE)), "REMARK 290 Symmetry line number problem");
		if (lineLength >= E_SYMMINDEX)   sym.D_SYMMINDEX       = MslTools::toInt(MslTools::trim(_pdbSymLine.substr(S_SYMMINDEX,L_SYMMINDEX)), "REMARK 290 Symmetry matrix index problem");

		if (lineLength >= E_SYMMX)       sym.D_SYMMX           = MslTools::toDouble(MslTools::trim(_pdbSymLine.substr(S_SYMMX,L_SYMMX)), "REMARK 290 Symmetry X problem");
		if (lineLength >= E_SYMMY)       sym.D_SYMMY           = MslTools::toDouble(MslTools::trim(_pdbSymLine.substr(S_SYMMY,L_SYMMY)), "REMARK 290 Symmetry Y problem");
		if (lineLength >= E_SYMMZ)       sym.D_SYMMZ           = MslTools::toDouble(MslTools::trim(_pdbSymLine.substr(S_SYMMZ,L_SYMMZ)), "REMARK 290 Symmetry Z problem");
		if (lineLength >= E_SYMTRANS)    sym.D_SYMTRANS        = MslTools::toDouble(MslTools::trim(_pdbSymLine.substr(S_SYMTRANS,L_SYMTRANS)), "REMARK 290 Symmetry trans problem");


	} catch(exception &e){
		cerr << "ERROR 34918 PDBFormat parseSymLine "<<e.what()<<endl;
		exit(34918);
	}

	return sym;
}
PDBFormat::BioUData PDBFormat::parseBioULine(const string &_pdbBioULine){
	BioUData bio;

	int lineLength = _pdbBioULine.size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_BIOURECORD)     strcpy(bio.D_BIOURECORD,    MslTools::trim(_pdbBioULine.substr(S_BIOURECORD, L_BIOURECORD)).c_str());

		if (lineLength >= E_BIOULINE)    bio.D_BIOULINE        = MslTools::toInt(MslTools::trim(_pdbBioULine.substr(S_BIOULINE,L_BIOULINE)), "REMARK 350 BioUnit line number problem");
		if (lineLength >= E_BIOUINDEX)   bio.D_BIOUINDEX       = MslTools::toInt(MslTools::trim(_pdbBioULine.substr(S_BIOUINDEX,L_BIOUINDEX)), "REMARK 350 BioUnit matrix index problem");

		if (lineLength >= E_BIOUX)       bio.D_BIOUX           = MslTools::toDouble(MslTools::trim(_pdbBioULine.substr(S_BIOUX,L_BIOUX)), "REMARK 350 BioUnit X problem");
		if (lineLength >= E_BIOUY)       bio.D_BIOUY           = MslTools::toDouble(MslTools::trim(_pdbBioULine.substr(S_BIOUY,L_BIOUY)), "REMARK 350 BioUnit Y problem");
		if (lineLength >= E_BIOUZ)       bio.D_BIOUZ           = MslTools::toDouble(MslTools::trim(_pdbBioULine.substr(S_BIOUZ,L_BIOUZ)), "REMARK 350 BioUnit Z problem");
		if (lineLength >= E_BIOUTRANS)   bio.D_BIOUTRANS       = MslTools::toDouble(MslTools::trim(_pdbBioULine.substr(S_BIOUTRANS,L_BIOUTRANS)), "REMARK 350 BioUnit trans problem");


	} catch(exception &e){
		cerr << "ERROR 34918 PDBFormat parseBioULine "<<e.what()<<endl;
		exit(34918);
	}

	return bio;
}

PDBFormat::AtomData PDBFormat::parseAtomLine(const string &_pdbAtomLine){
	AtomData atom;

	int lineLength = _pdbAtomLine.size();
	try {
		// Make sure line is long enough for each field. 
		if (lineLength >= E_RECORD_NAME)    strcpy(atom.D_RECORD_NAME,    MslTools::trim(_pdbAtomLine.substr(S_RECORD_NAME, L_RECORD_NAME)).c_str());
		if (lineLength >= E_ATOM_NAME)      strcpy(atom.D_ATOM_NAME,      MslTools::trim(_pdbAtomLine.substr(S_ATOM_NAME, L_ATOM_NAME)).c_str());
		if (lineLength >= E_RES_NAME)       strcpy(atom.D_RES_NAME,       MslTools::trim(_pdbAtomLine.substr(S_RES_NAME, L_RES_NAME)).c_str());


		
		if (lineLength >= E_ALT_LOC)        strcpy(atom.D_ALT_LOC,        MslTools::trim(_pdbAtomLine.substr(S_ALT_LOC,L_ALT_LOC)).c_str());
		if (lineLength >= E_CHAIN_ID)       strcpy(atom.D_CHAIN_ID,       MslTools::trim(_pdbAtomLine.substr(S_CHAIN_ID, L_CHAIN_ID)).c_str());
		if (lineLength >= E_I_CODE)         strcpy(atom.D_I_CODE,         MslTools::trim(_pdbAtomLine.substr(S_I_CODE, L_I_CODE)).c_str());
		if (lineLength >= E_SEG_ID)         strcpy(atom.D_SEG_ID,         MslTools::trim(_pdbAtomLine.substr(S_SEG_ID, L_SEG_ID)).c_str());


		  /*
		    Extracting Element from PDBLINE: asenes email:
		    1. Use the element field from the input PDB file as the element.
		    2. Else, if the first character of the atom name is a space (as " CA ") or a number (as "1HG1"), use the second character of the atom name as the element
		    3. Else, if the first character is a H and the last two character are number ("HG21"), take the first character (this excludes a case like a mercury "HG  ")
		    4. Else (a calcium, "CA  "), use the first two characters of the atom name as the element. 
		  */
		 if (lineLength >= E_ELEMENT_SYMBOL) {

		   string element = _pdbAtomLine.substr(S_ELEMENT_SYMBOL, L_ELEMENT_SYMBOL);

		   if (MslTools::trim(element) == ""){
		     string atomname = _pdbAtomLine.substr(S_ATOM_NAME, L_ATOM_NAME);
		     if (isdigit(atomname[0]) ||  atomname[0] == ' '){
		       element = MslTools::stringf("%c",atomname[1]);
		     } else if (atomname[0] == 'H'  && isdigit(atomname[2]) && isdigit(atomname[3])){
		       element = MslTools::stringf("%c",atomname[0]);
		     } else {
		       element = atomname.substr(0,2);
		     }
		   } 

		   // Copy extracted element to D_ELEMENT_SYMBOL
		   strcpy(atom.D_ELEMENT_SYMBOL, MslTools::trim(element).c_str());
		 }


		if (lineLength >= E_SERIAL)    atom.D_SERIAL         = MslTools::toInt(MslTools::trim(_pdbAtomLine.substr(S_SERIAL,L_SERIAL)), "Atom number(Serial) problem");
		if (lineLength >= E_RES_SEQ)   atom.D_RES_SEQ        = MslTools::toInt(MslTools::trim(_pdbAtomLine.substr(S_RES_SEQ,L_RES_SEQ)), "Residue Number problem");

		if (lineLength >= E_X)         atom.D_X              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_X,L_X)), "X-coord problem");
		if (lineLength >= E_Y)         atom.D_Y              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_Y,L_Y)), "Y-coord problem");
		if (lineLength >= E_Z)         atom.D_Z              = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_Z,L_Z)), "Z-coord problem");
		if (lineLength >= E_OCCUP)     atom.D_OCCUP          = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_OCCUP, L_OCCUP)), "Occupation problem");

		


		// Charge is  a problem
	//	if (lineLength >= E_CHARGE && MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)) != "")    atom.D_CHARGE         = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)), "CHARGE problem");
		try {
			if (lineLength >= E_CHARGE && MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)) != "")    atom.D_CHARGE         = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_CHARGE, L_CHARGE)), "CHARGE problem");
			if (lineLength >= E_TEMP_FACT) atom.D_TEMP_FACT      = MslTools::toDouble(MslTools::trim(_pdbAtomLine.substr(S_TEMP_FACT, L_TEMP_FACT)), "Temp. Factor problem");
		} catch(exception &e){
			// perhaps we should let the charge fail silently since too many PDBs will have strange stuff there
			//cerr << "WARNING 34923 PDBFormat parseAtomLine "<<e.what()<<endl;
		}
			


	} catch(exception &e){

	        cerr << "ERROR 34918 PDBFormat parseAtomLine "<<e.what()<<" line is: "<<_pdbAtomLine<<endl;
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
