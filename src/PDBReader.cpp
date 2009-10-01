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


	for (uint i = 0; i < symmetryRotations.size();i++){
		delete(symmetryRotations[i]);
	}

	for (uint i = 0; i < symmetryTranslations.size();i++){
		delete(symmetryTranslations[i]);
	}

	for (uint i = 0; i < biounitRotations.size();i++){
		delete(biounitRotations[i]);
	}

	for (uint i = 0; i < biounitTranslations.size();i++){
		delete(biounitTranslations[i]);
	}

	delete(scaleRotation);
	delete(scaleTranslation);

}


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
	//		string header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
			// check the length
			string header = "";
			if (line.size() >= PDBFormat::S_RECORD_NAME + PDBFormat::L_RECORD_NAME) {
				header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
			}

		// Deal with remark parsing..
			if (header == "REMARK"){
				// REMEMBER TO VALIDATE THE SUBSTR!!! (most important, the skip cannot be more than the srting length)
				if (line.size() >= PDBFormat::S_SYMMRECORD + PDBFormat::L_SYMMRECORD && line.substr(7,3) == "290"){
					string symlinetype = line.substr(PDBFormat::S_SYMMRECORD,PDBFormat::L_SYMMRECORD);
					if (symlinetype != "SMTRY"){
						continue;
					}

					PDBFormat::SymData sym = PDBFormat::parseSymLine(line);

					if (symmetryRotations.size() < sym.D_SYMMINDEX){
						symmetryRotations.push_back(new Matrix(3,3,0.0));
						(*symmetryRotations.back())[0][0] = 1.0;
						(*symmetryRotations.back())[1][1] = 1.0;
						(*symmetryRotations.back())[2][2] = 1.0;
						symmetryTranslations.push_back(new CartesianPoint(0.0,0.0,0.0));
					}
					(*symmetryRotations[sym.D_SYMMINDEX-1])[sym.D_SYMMLINE-1][0] = sym.D_SYMMX;
					(*symmetryRotations[sym.D_SYMMINDEX-1])[sym.D_SYMMLINE-1][1] = sym.D_SYMMY;
					(*symmetryRotations[sym.D_SYMMINDEX-1])[sym.D_SYMMLINE-1][2] = sym.D_SYMMZ;
					(*symmetryTranslations[sym.D_SYMMINDEX-1])[sym.D_SYMMLINE-1] = sym.D_SYMTRANS;
					
				}
				if (line.size() >= PDBFormat::S_BIOURECORD + PDBFormat::L_BIOURECORD && line.substr(7,3) == "350"){
					/*
					  This does not handle multiple BIOMT sections for different chains.
					  Therefore BIO UNIT matrices are not stored properly and BIO UNITS will not be properly generated.
					  See 3DVH as an example (its in testData.h)
					 */
					string biolinetype = line.substr(PDBFormat::S_BIOURECORD,PDBFormat::L_BIOURECORD);
					if (biolinetype != "BIOMT"){
						continue;
					}

					PDBFormat::BioUData bio = PDBFormat::parseBioULine(line);

					if (biounitRotations.size() < bio.D_BIOUINDEX){
						biounitRotations.push_back(new Matrix(3,3,0.0));
						(*biounitRotations.back())[0][0] = 1.0;
						(*biounitRotations.back())[1][1] = 1.0;
						(*biounitRotations.back())[2][2] = 1.0;
						biounitTranslations.push_back(new CartesianPoint(0.0,0.0,0.0));
					}
					(*biounitRotations[bio.D_BIOUINDEX-1])[bio.D_BIOULINE-1][0] = bio.D_BIOUX;
					(*biounitRotations[bio.D_BIOUINDEX-1])[bio.D_BIOULINE-1][1] = bio.D_BIOUY;
					(*biounitRotations[bio.D_BIOUINDEX-1])[bio.D_BIOULINE-1][2] = bio.D_BIOUZ;

					(*biounitTranslations[bio.D_BIOUINDEX-1])[bio.D_BIOULINE-1] = bio.D_BIOUTRANS;
				}
				// NOTE: CHANGE TO USE PDBFormat FOR MISSING ATOMS AND RESIDUES!!!!
				if(line.size() >= 27 && line.substr(7,3) == "465") {
					// missing residues

					if(line.substr(11,1) == " " && line.substr(15,1) != " " && line.substr(13,1) != "M") {
						MissingResidue res;	
						res.model = line.substr(13,1) == " " ? 0 : MslTools::toInt(line.substr(13,1));
						res.resName = line.substr(15,3);
						res.chainId = line.substr(19,1);
						res.resNum = MslTools::toInt(MslTools::trim(line.substr(22,4)));
						res.resIcode = line.substr(26,1);
						misRes.push_back(res);
					} else {
						continue;
					}
				}

				if(line.size() >= 27 && line.substr(7,3) == "470") {
					// missing atoms
					if(line.substr(11,1) == " " && line.substr(15,1) != " " && line.substr(13,1) != "M") {
						MissingAtoms mAtoms;	
						mAtoms.model = line.substr(13,1) == " " ? 0 : MslTools::toInt(line.substr(13,1));
						mAtoms.resName = line.substr(15,3);
						mAtoms.chainId = line.substr(19,1);
						mAtoms.resNum = MslTools::toInt(line.substr(20,4));
						mAtoms.resIcode = line.substr(24,1);
						string temp = line.substr(27);
						mAtoms.atoms = MslTools::tokenize(MslTools::trim(temp));
						misAtoms.push_back(mAtoms);

					} else {
						continue;
					}
				}
			}

			if (header == "SCALE1" || header == "SCALE2" || header == "SCALE3"){
					PDBFormat::ScaleData scale = PDBFormat::parseScaleLine(line);
					
					(*scaleRotation)[scale.D_SCALELINE-1][0] = scale.D_SCALEX;
					(*scaleRotation)[scale.D_SCALELINE-1][1] = scale.D_SCALEY;
					(*scaleRotation)[scale.D_SCALELINE-1][2] = scale.D_SCALEZ;

					(*scaleTranslation)[scale.D_SCALELINE-1] = scale.D_SCALETRANS;
					
			}
			if (header == "CRYST1"){
				PDBFormat::CrystData cryst = PDBFormat::parseCrystLine(line);

				unitCellParams.push_back(cryst.D_CRYSTA);
				unitCellParams.push_back(cryst.D_CRYSTB);
				unitCellParams.push_back(cryst.D_CRYSTC);
				unitCellParams.push_back(cryst.D_CRYSTALPHA);
				unitCellParams.push_back(cryst.D_CRYSTBETA);
				unitCellParams.push_back(cryst.D_CRYSTGAMMA);
				
			}
			
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


					if (atom.D_X < boundingCoords["minX"]){
						boundingCoords["minX"] = atom.D_X;
					}
					if (atom.D_X > boundingCoords["maxX"]){
						boundingCoords["maxX"] = atom.D_X;
					}

					if (atom.D_Y < boundingCoords["minY"]){
						boundingCoords["minY"] = atom.D_Y;
					}
					if (atom.D_Y > boundingCoords["maxY"]){
						boundingCoords["maxY"] = atom.D_Y;
					}

					if (atom.D_Z < boundingCoords["minZ"]){
						boundingCoords["minZ"] = atom.D_Z;
					}
					if (atom.D_Z > boundingCoords["maxZ"]){
						boundingCoords["maxZ"] = atom.D_Z;
					}

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

		boundingCoords["deltaX"] = boundingCoords["maxX"] - boundingCoords["minX"];
		boundingCoords["deltaY"] = boundingCoords["maxY"] - boundingCoords["minY"];
		boundingCoords["deltaZ"] = boundingCoords["maxZ"] - boundingCoords["minZ"];
		boundingCoords["maxDelta"] = boundingCoords["deltaX"];
		if (boundingCoords["deltaY"] > boundingCoords["deltaX"]){
			boundingCoords["maxDelta"] = boundingCoords["deltaY"];
		} else if (boundingCoords["deltaZ"] > boundingCoords["deltaY"]){
			boundingCoords["maxDelta"] = boundingCoords["deltaZ"];
		}

		return true;

		// Add last residue..
		if (singleAltLocFlag){
			string altLoc = discoverProperResidueAltLoc(currentResidueAtoms);

			addAtoms(currentResidueAtoms, altLoc);
			
		}

		
	} catch(...){
		cerr << "ERROR 5623 in PDBReader::read()\n";
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
