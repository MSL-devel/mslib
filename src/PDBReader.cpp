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

#include "PDBReader.h"


/**
 * This method will delete all of the Atom pointers
 * held in this PDBReader.
 */
void PDBReader::deletePointers() {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
}



/// @todo Why are these two read methods commented out in PDBReader?
/*
bool PDBReader::read(){
	//cout << "Read function"<<endl;
	return true;
}

bool PDBReader::read(vector<CartesianPoint> &_cv){

	_cv.clear();

	try { 
		while (!endOfFileTest()){

			string line = Reader::getLine();
			// DO NOT TRIM: CHECK FOR "ATOM  " AND "HETATM" (THE WHOLE 6 CHARACTERS)
			//string header 	= MslTools::trim(line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME));
			string header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
			
			if (header == "ATOM  " || header == "HETATM"){
				//cout << "LINE: "<<line<<endl;

				PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line);


				CartesianPoint c;
				c.setX(atom.D_X);
				c.setY(atom.D_Y);
				c.setZ(atom.D_Z);

				_cv.push_back(c);

			}


		}

	} catch(...){
		cerr << "ERROR 5623 in PDBReader::read(vector<CartesianPoint> &_cv)\n";
		exit(5623);
	}


	//cout << "Number of Points Read in : "<<_cv.size()<<endl;
	return true;
}
*/


/*
bool PDBReader::read(AtomVector &_av){

	try { 
		while (!endOfFileTest()){

			string line = Reader::getLine();
			//string header 	= MslTools::trim(line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME));
			string header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
			
			if (header == "ATOM  " || header == "HETATM"){
				PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line);

				Atom *a = new Atom();
				a->setName(atom.D_ATOM_NAME);
				a->setCoor(atom.D_X,atom.D_Y, atom.D_Z);


				_av.push_back(a);
				a = NULL;

			}


		}

	} catch(...){
		cerr << "ERROR 5623 in PDBReader::read(vector<CartesianPoint> &_cv)\n";
		exit(5623);
	}

	return true;
}
*/

/**
 * This method will actually read the data in from the
 * given PDB file.  Currently the reader only looks at 
 * ATOM and HETATM information.  All other information 
 * (REMARKS, HEADER, TITLE, SEQRES, HET, etc.) is
 * simply ignored.
 */
bool PDBReader::read() {

	try { 


		string currentResidue = "";
		map<string, vector<PDBFormat::AtomData> > currentResidueAtoms;

		while (!endOfFileTest()){
			string line = Reader::getLine();
			string header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);


			
			if (header == "ATOM  " || header == "HETATM"){
				PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line);



				// NORMAL READING MODE = singleAltLocFlag is false;
				if (!singleAltLocFlag){
					atoms.push_back(new Atom);

					// atom name, residue name, residue icode, chain id, coor, element
					atoms.back()->setName(atom.D_ATOM_NAME);
					atoms.back()->setResidueName(atom.D_RES_NAME);
					atoms.back()->setResidueIcode(atom.D_I_CODE);
					atoms.back()->setResidueNumber(atom.D_RES_SEQ);
					atoms.back()->setChainId(atom.D_CHAIN_ID);
					atoms.back()->setCoor(atom.D_X,atom.D_Y, atom.D_Z);
					atoms.back()->setElement(atom.D_ELEMENT_SYMBOL);
					atoms.back()->setTempFactor(atom.D_TEMP_FACT);

				} else {

					
					/*
					  Try to determine a single Alternate Location for each residue.
					*/


					// Residue description string
					stringstream resDescription;
					resDescription << atom.D_CHAIN_ID <<":"<<atom.D_RES_NAME<<":"<<atom.D_RES_SEQ<<":"<<atom.D_I_CODE;

					// New residue means its time to decide which atoms will get added to "atoms"
					if (currentResidue != resDescription.str()){


						string altLoc = discoverProperResidueAltLoc(currentResidueAtoms);

						addAtoms(currentResidueAtoms, altLoc);
						
						currentResidueAtoms.clear();
					}

					currentResidueAtoms[atom.D_ATOM_NAME].push_back(atom);
					currentResidue = resDescription.str();
				}

			}

		}
		return true;

		// Add last residue..
		if (singleAltLocFlag){
			string altLoc = discoverProperResidueAltLoc(currentResidueAtoms);

			addAtoms(currentResidueAtoms, altLoc);
			
		}

		
	} catch(...){
		cerr << "ERROR 5623 in PDBReader::read(vector<CartesianPoint> &_cv)\n";
		return false;
	}

	return true;
}



void PDBReader::addAtoms(map<string, vector<PDBFormat::AtomData> > &_currentResidueAtoms, string _altLoc){

	// Convert AtomData to Atom and add to "atoms"
	map<string,vector<PDBFormat::AtomData> >::iterator it;
	for (it = _currentResidueAtoms.begin();it != _currentResidueAtoms.end();it++){

		int addIt = -1;
		if (it->second.size() == 1){
			addIt = 0;
		} else {
			for (uint i = 0; i < it->second.size();i++){	
				if (it->second[i].D_ALT_LOC == _altLoc){
					addIt = i;
					break;
				}
			}

		}


		if (addIt > -1){
			atoms.push_back(new Atom);

			// atom name, residue name, residue icode, chain id, coor, element
			atoms.back()->setName(it->second[addIt].D_ATOM_NAME);
			atoms.back()->setResidueName(it->second[addIt].D_RES_NAME);
			atoms.back()->setResidueIcode(it->second[addIt].D_I_CODE);
			atoms.back()->setResidueNumber(it->second[addIt].D_RES_SEQ);
			atoms.back()->setChainId(it->second[addIt].D_CHAIN_ID);
			atoms.back()->setCoor(it->second[addIt].D_X,it->second[addIt].D_Y, it->second[addIt].D_Z);
			atoms.back()->setElement(it->second[addIt].D_ELEMENT_SYMBOL);
			atoms.back()->setTempFactor(it->second[addIt].D_TEMP_FACT);
		}

	}
}
string PDBReader::discoverProperResidueAltLoc(map<string, vector<PDBFormat::AtomData> > &_currentResidueAtoms){

	// For each atom that has > 1 instances, add up occupancies
	map<string,double> occ;
	map<string,vector<PDBFormat::AtomData> >::iterator it;
	for (it = _currentResidueAtoms.begin();it != _currentResidueAtoms.end();it++){
		if (it->second.size() == 1){
			continue;
		}
							
		for (uint i = 0; i < it->second.size();i++){	
			occ[it->second[i].D_ALT_LOC] += it->second[i].D_OCCUP;
			//cout << "AltLoc Value: '"<<it->second[i].D_ALT_LOC<<"' is "<<it->second[i].D_OCCUP<<endl;
		}
	}
	string altLoc = "";
	double occupancy = 0.0;
	map<string,double>::iterator it2;
	for (it2 = occ.begin();it2 != occ.end();it2++){
		if (it2->second > occupancy){
			altLoc = it2->first;
			occupancy = it2->second;
		}
	}


	return altLoc;
}
