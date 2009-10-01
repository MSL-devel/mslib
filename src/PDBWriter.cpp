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

#include "PDBWriter.h"



bool PDBWriter::write(vector<CartesianPoint> &_cv){
	int atomCount = 1;
	int resCount  = 1;
	for (vector<CartesianPoint>::iterator it = _cv.begin(); it != _cv.end(); it++){
		PDBFormat::AtomData atom;

		atom.D_X = it->getX();
		atom.D_Y = it->getY();
		atom.D_Z = it->getZ();

		strncpy(atom.D_RECORD_NAME , "ATOM  ",PDBFormat::L_RECORD_NAME);
		strncpy(atom.D_ATOM_NAME, " DUM", PDBFormat::L_ATOM_NAME);
		strncpy(atom.D_RES_NAME, "DUM", PDBFormat::L_RES_NAME);
		strncpy(atom.D_CHAIN_ID, "A", PDBFormat::L_CHAIN_ID);
		strncpy(atom.D_ELEMENT_SYMBOL, "C", PDBFormat::L_ELEMENT_SYMBOL);

		// Allow for more than 10000 points, but no check for more than 10000 residues
		//   thus 100,000,000 atom limit!
		if ( atomCount / 10000 >= 1) { 
			resCount++;
			atomCount = 1;
		}		
		atom.D_SERIAL  = atomCount++;
		atom.D_RES_SEQ = resCount;
		atom.D_OCCUP   = 1.0;
		
		string pdbline = PDBFormat::createAtomLine(atom);

		writeln(pdbline);
		
	}

	return true;
}


bool PDBWriter::write(AtomVector &_av, bool _addTerm, bool _noHydrogens,bool _writeAsModel) {

	/******************************************************
	 *
	 *   TODO:  introduce support for alt location for pdb 
	 *   with multiple coord
	 *   :w
	 *
	 ******************************************************/

    if( is_open() == false )
       return false;


    if (_writeAsModel && _av.size() > 0){
	    string model = "MODEL";
	    writeln(model);
    }

	int atomCount = 1;
	for (AtomVector::iterator it = _av.begin(); it != _av.end(); it++){

		PDBFormat::AtomData atom;

		/*
		atom.D_X = (*it)->getX();
		atom.D_Y = (*it)->getY();
		atom.D_Z = (*it)->getZ();

		// TODO: Need to check getName size < L_ATOM_NAME ???, does it truncate automatically
		strncpy(atom.D_ATOM_NAME, (*it)->getName().c_str(), PDBFormat::L_ATOM_NAME);


		strncpy(atom.D_RECORD_NAME , "ATOM  ",PDBFormat::L_RECORD_NAME);
		strncpy(atom.D_RES_NAME, (*it)->getResidueName().c_str(), PDBFormat::L_RES_NAME);
		// 
		strncpy(atom.D_CHAIN_ID, (*it)->getChainId().c_str(), PDBFormat::L_CHAIN_ID);
		strncpy(atom.D_I_CODE, (*it)->getResidueIcode().c_str(), PDBFormat::L_I_CODE);
		strncpy(atom.D_ELEMENT_SYMBOL, (*it)->getElement().c_str(), PDBFormat::L_ELEMENT_SYMBOL);
		*/
		if ( _noHydrogens && (*it)->getElement() == "H") continue;

		const Atom *a = (*it);
		atom = PDBFormat::createAtomData(*a);
		// Allow for more than 10000 points, but no check for more than 10000 residues
		//   thus 100,000,000 atom limit!
		if ( atomCount / 10000 >= 1) { 
			atomCount = 1;
		}		
		atom.D_SERIAL  = atomCount++;
		
		
		string pdbline = PDBFormat::createAtomLine(atom);

		//writeln(pdbline);
		if (!writeln(pdbline)) {
			cerr << "WARNING 12491: cannot write atom line in bool PDBWriter::write(AtomVector &_av, bool _addTerm, bool _noHydrogens,bool _writeAsModel)" << endl;
			return false;
		}

		// add a TER line at the end of each chain
		if ( ((it+1 == _av.end()) && _addTerm) || (it+1 != _av.end() && (*(it+1))->getChainId() !=  atom.D_CHAIN_ID)) {
			PDBFormat::AtomData ter;
			ter.D_SERIAL  = atomCount++;
			ter.D_RES_SEQ = (*it)->getResidueNumber();
			strncpy(ter.D_RECORD_NAME , "TER   ",PDBFormat::L_RECORD_NAME);
			strncpy(ter.D_RES_NAME, (*it)->getResidueName().c_str(), PDBFormat::L_RES_NAME);
			strncpy(ter.D_CHAIN_ID, (*it)->getChainId().c_str(), PDBFormat::L_CHAIN_ID);
			strncpy(ter.D_I_CODE, (*it)->getResidueIcode().c_str(), PDBFormat::L_I_CODE);
			pdbline = PDBFormat::createTerLine(ter);
			//writeln(pdbline);
			if (!writeln(pdbline)) {
				cerr << "WARNING 12496: cannot write ter line in bool PDBWriter::write(AtomVector &_av, bool _addTerm, bool _noHydrogens,bool _writeAsModel)" << endl;
				return false;
			}
		}
		
	}
	if (_writeAsModel && _av.size() > 0){
		string endmdl = "ENDMDL";
		writeln(endmdl);
	}

	return true;
}

void PDBWriter::writeREMARKS(){

	if (fileHandler != stringstyle && !fileStream.is_open()){
		cout << "File : "<<fileName<< " is not open for writing\n";
		return;
	}

	// CORRECTED FORMAT.  DO WE NEED THIS CREDIT LINE?
	//string credit = "REMARK 000 File written by PDBWriter, which is part of the MSL libraries.";
        string credit = "REMARK   0 File written by PDBWriter, which is part of the MSL libraries.       ";
	writeln(credit);


	
	vector<string>::iterator it;
	for (it = remarks.begin(); it != remarks.end(); it++){
		vector<string> singleRemarks = MslTools::tokenize(*it, "\n");
		for (uint i = 0;i < singleRemarks.size();i++){
			if (singleRemarks[i].find("from pymol") == string::npos){
				// WHAT IS THIS "from pymol" FOR??????????
				string line = "REMARK "+singleRemarks[i];
				writeln(line);
			}

		}
	}
}
