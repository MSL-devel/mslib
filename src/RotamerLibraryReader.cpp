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


#include "RotamerLibraryReader.h"

using namespace MSL;
using namespace std;


RotamerLibraryReader::RotamerLibraryReader() : Reader() {
	setup(NULL);
}

RotamerLibraryReader::RotamerLibraryReader(const string &_filename, RotamerLibrary * _rotlib) : Reader(_filename) {
	setup(_rotlib);
}

RotamerLibraryReader::RotamerLibraryReader(RotamerLibrary * _rotlib) : Reader() {
	setup(_rotlib);
}

RotamerLibraryReader::RotamerLibraryReader(const RotamerLibraryReader & _reader) {
	copy(_reader);
}

RotamerLibraryReader::~RotamerLibraryReader() {
	close();
}

void RotamerLibraryReader::setup(RotamerLibrary * _pRotlib) {
	pRotLib = _pRotlib;
}

void RotamerLibraryReader::copy(const RotamerLibraryReader & _reader) {
	pRotLib = _reader.pRotLib;
}


bool RotamerLibraryReader::read() {
	if (!is_open()) {
		return false;
	}


	if (pRotLib == NULL) {
		return false;
	}

	try { 
		string currentLib = "";
		string currentRes = "";
		bool foundLib = false;
		bool foundRes = false;
		unsigned int lineCounter = 0;

		vector<string> levRes; // store the order of residues in LEVRES line
		bool foundLevRes = false; // Have we seen a LEVRES line ?
	
		while (!endOfFileTest()){
			lineCounter++;

			string line = Reader::getLine();
			//cout << line << endl;
			line = MslTools::uncomment(line);
			vector<string> tokens = MslTools::tokenize(line);
			if (tokens.size() == 0) {
				continue;
			}

			// found the library name
			if (tokens[0] == "LIBRARY") {
				//cout << "UUU in LIBRARY line" << endl;
				if (tokens.size() > 1 &&  tokens[1].size() > 0) {
					currentLib = tokens[1];
					pRotLib->addLibrary(tokens[1]);
					//libraries[currentLib].clear();
					foundLib = true;
					foundRes = false;
					//if (defaultLibrary == "") {
					//	defaultLibrary = currentLib;
					//}
					continue;
				} else {
					cerr << "ERROR 5010: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
					return false;
				}
			}
			
			if(tokens[0] == "LEVRES") {
				foundLevRes = true;
				tokens.erase(tokens.begin());
				levRes = tokens;
				continue;
			}

			if(tokens[0] == "LEVEL" ) {
				if(foundLevRes) {
					if(tokens.size() > 2) {
						string levelName = tokens[1];
						tokens.erase(tokens.begin(),tokens.begin()+2);
						for(int i = 0; i < tokens.size(); i++) {
							pRotLib->setLevel(levelName,levRes[i],MslTools::toUnsignedInt(tokens[i]));
						}
					}
				} else {
					cerr << "ERROR 5014: A \"LEVRES\" line needs to appear before any \"LEVEL\" line in the rotamer library File." << endl;
					return false;
				}
				continue;
			}
			
			// found the residue name
			if (tokens[0] == "RESI") {
				//cout << "UUU in RESI line" << endl;
				if (!foundLib) {
					cerr << "ERROR 5015: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
				return false;
				}
				if (tokens.size() > 1 &&  tokens[1].size() > 0) {
					// add the residue and the library to the object
					//addResidue(tokens[1]);
					currentRes = tokens[1];
					//libraries[currentLib][currentRes].mobileAtoms.clear();
					//libraries[currentLib][currentRes].defi.clear();
					//libraries[currentLib][currentRes].internalCoor.clear();
					pRotLib->addResidue(currentLib, tokens[1]);
					foundRes = true;
				} else {
					cerr << "ERROR 5016: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
					return false;
				}
				continue;
			}

			// list of atoms to be initialized
			if (foundRes && (tokens[0] == "MOBI" || tokens[0] == "INIT")) {
				//cout << "UUU in INIT line" << endl;
				if (!foundRes) {
					cerr << "ERROR 5017: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
					return false;
				}
				tokens.erase(tokens.begin());
				//libraries[currentLib][currentRes].mobileAtoms = tokens;
				pRotLib->addMobileAtoms(currentLib, currentRes, tokens);
			}

			// found a bond, angle, improper or dihedral definiton
			if (foundRes && tokens[0] == "DEFI") {
				//cout << "UUU in DEFI line" << endl;
				if (tokens.size() < 3 || tokens.size() > 5) {
					cerr << "ERROR 5018: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
					return false;
				}
				tokens.erase(tokens.begin()); // remove the DEFI token
				while (tokens.size() < 4) {
					// add one or two blanks if it was a bond or angle defintion
					tokens.push_back("");
				}
				pRotLib->addInternalCoorDefinition(currentLib, currentRes, tokens);
			//	cout << "UUU " << libraries[currentLib][currentRes].defi.back().type << endl;
			//	for (unsigned int i=0; i<libraries[currentLib][currentRes].defi.back().atomNames.size(); i++) {
			//		cout << " UUU - " << libraries[currentLib][currentRes].defi.back().atomNames[i] << endl;
			//	}
			//	for (unsigned int i=0; i<libraries[currentLib][currentRes].defi.back().resnumCorrectors.size(); i++) {
			//		cout << " UUU -> " << libraries[currentLib][currentRes].defi.back().resnumCorrectors[i] << endl;
			//	}
			}
			// found a bond, angle, improper or dihedral definiton using ICDEF
			if (foundRes && tokens[0] == "ICDEF") {
				//cout << "UUU in DEFI line" << endl;
				if (tokens.size() < 5) {
					cerr << "ERROR 5019: syntax error in rotamer library " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
					return false;
				}
				tokens.erase(tokens.begin()); // remove the ICDEF token

				// Format ICDEF A1 A2 A3 A4 => <Dist A3-A4> <ANGLE A2-A3-A4> <DIHEDRAL/IMPROPER A1-A2-A3-A4>
				vector<string> tmpTokens;
				if(tokens[2].substr(0,1) == "*") {
					// removing the * from distance,angle internalCoorDefinitions
					tmpTokens.push_back(tokens[2].substr(1));
				} else {
					tmpTokens.push_back(tokens[2]);
				}
				tmpTokens.push_back(tokens[3]);

				tmpTokens.push_back("");
				tmpTokens.push_back("");
				pRotLib->addInternalCoorDefinition(currentLib, currentRes, tmpTokens); // bond distance
				
				tmpTokens.insert(tmpTokens.begin(),tokens[1]);
				tmpTokens.erase(tmpTokens.end());
				pRotLib->addInternalCoorDefinition(currentLib, currentRes, tmpTokens); // bond angle

				pRotLib->addInternalCoorDefinition(currentLib, currentRes, tokens); // dihedral
			//	cout << "UUU " << libraries[currentLib][currentRes].defi.back().type << endl;
			//	for (unsigned int i=0; i<libraries[currentLib][currentRes].defi.back().atomNames.size(); i++) {
			//		cout << " UUU - " << libraries[currentLib][currentRes].defi.back().atomNames[i] << endl;
			//	}
			//	for (unsigned int i=0; i<libraries[currentLib][currentRes].defi.back().resnumCorrectors.size(); i++) {
			//		cout << " UUU -> " << libraries[currentLib][currentRes].defi.back().resnumCorrectors[i] << endl;
			//	}
			}

			// found a conformation
			if (foundRes && tokens[0] == "CONF") {
				//cout << "UUU in CONF line" << endl;
				tokens.erase(tokens.begin());
				vector<double> doubleTokens;
				unsigned int rotamerBin = 0;
				int numInternalCoors = (pRotLib->getInternalCoorDefinition(currentLib,currentRes)).size();
				for (int i = 0; i < numInternalCoors; i++) {
					doubleTokens.push_back(MslTools::toDouble(tokens[i]));
				}
				if(tokens.size() == numInternalCoors + 1) {
					rotamerBin = MslTools::toUnsignedInt(tokens[numInternalCoors]); 
				}
				//libraries[currentLib][currentRes].internalCoor.push_back(doubleTokens);
				pRotLib->addConformation(currentLib, currentRes, doubleTokens,rotamerBin);
			//	for (unsigned int i=0; i<libraries[currentLib][currentRes].internalCoor.back().size(); i++) {
			//		cout << " UUU - " << libraries[currentLib][currentRes].internalCoor.back()[i] << endl;
			//	}
			}
		}

		if (!pRotLib->calculateBuildingICentries()) {
			cerr << "ERROR 5023: error in crating the IC building instructions " << fileName << " at line " << lineCounter << ", in bool RotamerLibraryReader::read()" << endl;
			return false;
		}


	} catch(...){
		cerr << "ERROR 5623 in RotamerLibraryReader::read(vector<CartesianPoint> &_cv)\n";
		return false;
	}

	return true;
}



