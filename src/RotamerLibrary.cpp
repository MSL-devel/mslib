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

#include "RotamerLibrary.h"

RotamerLibrary::RotamerLibrary() {
	setup();
}


RotamerLibrary::RotamerLibrary(const RotamerLibrary & _rotlib) {
}

RotamerLibrary::~RotamerLibrary() {
}

void RotamerLibrary::copy(const RotamerLibrary & _rotlib) {
	defaultLibrary = _rotlib.defaultLibrary;
	libraries = _rotlib.libraries;
}

void RotamerLibrary::setup() {
	lastFoundRes = libraries.begin()->second.begin();
	defaultLibrary = "";
}

bool RotamerLibrary::addInternalCoorDefinition(string _libName, string _resName, const vector<string> & _atoms) {

	if (_libName == "") {
		_libName = defaultLibrary;
	}

	if (libraries.find(_libName) == libraries.end()) {
		cerr << "WARNING 48234: unknown library named " << _libName << " in bool RotamerLibrary::addDegreeOfFreedom(string _libName, string _resName, vector<string> _atoms)" << endl;
		return false;
	} else if (libraries[_libName].find(_resName) == libraries[_libName].end()) {
		cerr << "WARNING 48239: unknown residue " << _resName << " in library " << _libName << " in bool RotamerLibrary::addDegreeOfFreedom(string _libName, string _resName, vector<string> _atoms)" << endl;
		return false;
	} else if (_atoms.size() < 2 || _atoms.size() > 4) {
		cerr << "WARNING 48239: invalid number of atoms " << _atoms.size() << " in bool RotamerLibrary::addDegreeOfFreedom(string _libName, string _resName, vector<string> _atoms)" << endl;
		return false;
	}

	/************************************************************
	 *  The function ASSUMES receving _atoms of size 4
	 *
	 *  Atoms may be given as -C or +N meaning the C of the previous
	 *  residue or the N of /he next one.  The resnumCorrectors
	 *  are the correction to be made to the absres when an atom identifier
	 *  is generated, and the name of the atom is stored as C or N in
	 *  the atomNames vector
	 *
	 *  This loop also checks and clean the * of the 3rd atom of impropers
	 ************************************************************/
	InternalCoorDefi D;
	//D.resnumCorrectors = vector<int>(4, 0);
	D.atomNames = _atoms;
	D.type = 0; // bond
	bool improper_flag = false;

	for (unsigned int i=0; i<4; i++) {
		if (D.atomNames[i] == "") {
			if (i<2) {
				// the first two atoms must not be blank
				return false;
			} else {
				// it is a bond or angle
				break;
			}
		}
		if (i==2) {
			D.type = 1; // angle
		} else if (i==3) {
			if (improper_flag) { 
				D.type = 3; // improper
			} else {
				D.type = 2; // dihedral
			}
		}
		// check for modefiers (+ atom belongs to next residue, - atom belongs to prev res
//		if (D.atomNames[i].substr(0,1) == "+") {
//			D.atomNames[i] = D.atomNames[i].substr(1,D.atomNames[i].size()-1);
//			D.resnumCorrectors[i] = 1;
//		} else if (D.atomNames[i].substr(0,1) == "-") {
//			D.atomNames[i] = D.atomNames[i].substr(1,D.atomNames[i].size()-1);
//			D.resnumCorrectors[i] = -1;
//		}
		if (i==2 && D.atomNames[i].size() > 0) {
			if (D.atomNames[i].substr(0,1) == "*") {
				// this is an improper dihedral
				improper_flag = true;
				D.atomNames[i] = D.atomNames[i].substr(1,D.atomNames[i].size()-1);
				continue;
			}
		}
	}
	if (D.type == 0) {
		D.atomNames.erase(D.atomNames.begin()+2, D.atomNames.end());
//		D.resnumCorrectors.erase(D.resnumCorrectors.begin()+2, D.resnumCorrectors.end());
	} else if (D.type == 1) {
		D.atomNames.erase(D.atomNames.begin()+3, D.atomNames.end());
//		D.resnumCorrectors.erase(D.resnumCorrectors.begin()+3, D.resnumCorrectors.end());
	}

	libraries[_libName][_resName].defi.push_back(D);
	return true;
}

bool RotamerLibrary::calculateBuildingICentries() {

	/*********************************************************************
	 *  This function calculates from the DEFI a series of IC entries
	 *  that can build the atoms that are moved by the rotamer
	 *
	 *  The IC entries are in the form of 4 atom names (the last being the one
	 *  of the atoms to be rebuilt) and the three values (distance, angle and
	 *  dihedral) that are needed to rebuild
	 *
	 *  The values vary for each conf and what it is stored is the index
	 *  of the value within each defi so that the conf can be looked up
	 *
	 *  For example:
	 *
	 *       RESI ALA
	 *       INIT CB HA HB1 HB2 HB3
	 *     0 DEFI N C *CA CB
	 *     1 DEFI N C *CA HA
	 *     2 DEFI C CA CB HB1
	 *     3 DEFI HB1 CA *CB HB2
	 *     4 DEFI HB1 CA *CB HB3
	 *     5 DEFI C CA CB
	 *     6 DEFI C CA HA
	 *     7 DEFI CA CB HB1
	 *     8 DEFI CA CB HB2
	 *     9 DEFI CA CB HB3
	 *    10 DEFI CA CB
	 *    11 DEFI CA HA
	 *    12 DEFI CB HB1
	 *    13 DEFI CB HB2
	 *    14 DEFI CB HB3
	 *       CONF   120.21 -120.45   60.00  119.13 -119.58  114.29  106.39  109.60  111.05  111.61    1.55    1.08    1.11    1.11    1.11
	 *       CONF   122.35 -120.45   60.00  119.13 -119.58  109.36  106.39  109.60  111.05  111.61    1.55    1.08    1.11    1.11    1.11
	 *       ...
	 *
	 *  Will produce
	 *       atom names        indices
	 *       -----------       --------
	 *    1) N   C   CA  CB    10  5  0 
	 *    2) N   C   CA  HA    11  6  1
	 *    3) C   CA  CB  HB1   12  7  2
	 *    4) HB1 CA  CB  HB2   13  8  3
	 *    5) HB1 CA  CB  HB3   14  9  4
	 *                   |      |  |  |
	 *       atom to build   bond ang dihe
	 *********************************************************************/

	for (map<string, map<string, Res> >::iterator libItr=libraries.begin(); libItr!=libraries.end(); libItr++) {
		// for each library
		for (map<string, Res>::iterator resItr=libItr->second.begin(); resItr!=libItr->second.end(); resItr++) {
			// for each residue

			// get the atoms that will be re-built
			vector<string> initAtoms = resItr->second.initAtoms;
			vector<InternalCoorDefi> defi = resItr->second.defi;
			for (unsigned int i=0; i<initAtoms.size(); i++) {
				// for each atom to rebuild
				vector<string> bond;
				vector<string> angle;
				vector<string> dihe;
				vector<unsigned int> indeces(3, 0);
				bool improper_flag = false;
				for (unsigned int j=0; j<defi.size(); j++) {
					// for each internal coor definition
					if (defi[j].type == 0 && defi[j].atomNames[1] == initAtoms[i]) {
						// found a bond where the last atom is this atom
						bond =  defi[j].atomNames;
						indeces[0] = j;
					} else if (defi[j].type == 1 && defi[j].atomNames[2] == initAtoms[i]) {
						// found an angle where the last atom is this atom
						angle =  defi[j].atomNames;
						indeces[1] = j;
					} else if ((defi[j].type == 2 || defi[j].type == 3) && defi[j].atomNames[3] == initAtoms[i]) {
						// found a dihe or improper where the last atom is this atom
						dihe =  defi[j].atomNames;
						indeces[2] = j;
						if (defi[j].type == 3) {
							improper_flag = true;
						} else {
							improper_flag = false;
						}
					}
				}

				// make sure we found bond, angle and dihe and that the atoms are
				// in common
				if (bond.size() == 0) {
					cerr << "WARNING 95123: bond not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << initAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
					return false;
				} else if (angle.size() == 0) {
					cerr << "WARNING 95128: angle not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << initAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				} else if (dihe.size() == 0) {
					cerr << "WARNING 95133: dihedral/improper not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << initAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				} else if (bond[0] != dihe[2]) {
					cerr << "WARNING 95138: bond atoms do not match dihe/improper building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << initAtoms[i] << " (" << bond[0] << " " << bond[1] << " != " << dihe[0] << " " << dihe[1] << " " << dihe[2] << " " << dihe[3] << "), in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				}else if (angle[0] != dihe[1] || angle[1] != dihe[2]) {
					cerr << "WARNING 95143: angle atoms do not match dihe/improper building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << initAtoms[i] << " (" << angle[0] << " " << angle[1] << " " << angle[2] << " != " << dihe[0] << " " << dihe[1] << " " << dihe[2] << " " << dihe[3] << "), in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				}
				resItr->second.buildingInstructions[initAtoms[i]].atomNames = dihe;
				resItr->second.buildingInstructions[initAtoms[i]].defiIndeces = indeces;
				resItr->second.buildingInstructions[initAtoms[i]].improper_flag = improper_flag;
			}
			/*
			cout << "UUU lib " << libItr->first << ", res " << resItr->first << ", atom " << endl;
			for (map<string, RotamerBuildingIC>::iterator bldItr = resItr->second.buildingInstructions.begin(); bldItr!=resItr->second.buildingInstructions.end(); bldItr++) {
			cout << "  atom " << bldItr->first << " " << endl;
				cout << "   ";
				for (vector<string>::iterator strItr=bldItr->second.atomNames.begin(); strItr!=bldItr->second.atomNames.end(); strItr++) {
					cout << " " << *strItr;
				}
				cout << " : ";
				for (vector<unsigned int>::iterator uiItr=bldItr->second.defiIndeces.begin(); uiItr!=bldItr->second.defiIndeces.end(); uiItr++) {
					cout << " " << *uiItr;
				}
				cout << endl;
				//bldItr->second
			}
			*/
		}
	}

	


	return true;
}


string RotamerLibrary::toString(){
	string result;

	
	// for each library
	for (map<string, map<string, Res> >::iterator libItr=libraries.begin(); libItr!=libraries.end(); libItr++) {

		result += "LIBRARY "+libItr->first+"\n\n";
		result += "CHARMMPAR 22 27\n\n";

		// for each residue
		for (map<string, Res>::iterator resItr=libItr->second.begin(); resItr!=libItr->second.end(); resItr++) {

			result += "RESI "+resItr->first+"\n";
			// Print INIT atoms
			vector<string> initAtoms = resItr->second.initAtoms;
			result += "INIT ";
			for (unsigned int i=0; i<initAtoms.size(); i++) {
				result +=  initAtoms[i] + " ";
			}

			result += "\n";

			vector<InternalCoorDefi> defi = resItr->second.defi;

			for (unsigned int i=0; i<defi.size(); i++) {
				result += "DEFI ";
				for (unsigned int j = 0; j < defi[i].atomNames.size();j++){
					result += defi[i].atomNames[j] + " ";
				}
				result += "\n";
			}

			vector<vector<double> > iCoor = resItr->second.internalCoor;
			for (unsigned int i = 0; i < iCoor.size();i++){
				result += "CONF ";
				for (unsigned int j = 0; j < iCoor[i].size();j++){
					char c[10];
					sprintf(c,"%8.3f",iCoor[i][j]);
					result += (string)c + " ";
				}
				result += "\n";
			}
			

			
		}
	}

	return result;
	
}
