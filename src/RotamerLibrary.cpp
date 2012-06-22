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

#include "RotamerLibrary.h"
#include "RotamerLibraryWriter.h"
#include "RotamerLibraryReader.h"
#include <stdio.h>

using namespace MSL;
using namespace std;


RotamerLibrary::RotamerLibrary() {
	setup();
}


RotamerLibrary::RotamerLibrary(const RotamerLibrary & _rotlib) {
	setup();
	copy(_rotlib);
}

RotamerLibrary::~RotamerLibrary() {
	deletePointers();
}

void RotamerLibrary::deletePointers() {
	libraries.clear();
	delete rotReader;
	rotReader = NULL;
	delete rotWriter;
	rotWriter = NULL;
}

void RotamerLibrary::copy(const RotamerLibrary & _rotlib) {
	defaultLibrary = _rotlib.defaultLibrary;
	libraries = _rotlib.libraries;
	levels = _rotlib.levels;
}

void RotamerLibrary::reset() {
	libraries.clear();
	levels.clear();
	lastFoundRes = libraries.begin()->second.begin();
	defaultLibrary = "";
}

void RotamerLibrary::setup() {
	reset();
	rotReader = new RotamerLibraryReader();
	rotReader->setRotamerLibrary(this);
	rotWriter = new RotamerLibraryWriter();
}
void RotamerLibrary::removeAllConformations () {
	for (map<string, map<string, Res> > ::iterator lib = libraries.begin(); lib != libraries.end(); lib++) {
		for (map <string, Res>::iterator resName = lib->second.begin(); resName != lib->second.end(); resName++) {
			libraries[lib->first][resName->first].internalCoor.clear();
			libraries[lib->first][resName->first].rotamerBins.clear();
			//libraries[lib->first][resName->first].buildingInstructions.clear();
		}
	}
}


bool RotamerLibrary::removeRotamer(string _libName,string _resName,int _num) {
	if(_libName == "") {
		_libName = defaultLibrary;
	}

	if (residueExists(_libName, _resName)) {
	//if((libraries.find(_libName) != libraries.end()) && (libraries[_libName].find(_resName) != libraries[_libName].end()) ) {
			if(libraries[_libName][_resName].internalCoor.size() > _num) {
				libraries[_libName][_resName].internalCoor.erase(libraries[_libName][_resName].internalCoor.begin() + _num);
				libraries[_libName][_resName].rotamerBins.erase(libraries[_libName][_resName].rotamerBins.begin() + _num);
				return true;
			} else {
				//cerr << "Warning 2320: Not enough conformers for the residue" << endl;
				return false;
			}
	} else {

		//cerr << "Warning 2321: Unable to find conformer to remove" << endl;
		return false;

	}	

}	

bool RotamerLibrary::removeRotamers(std::string _libName,std::string _resName,unsigned _startIdx, unsigned _numOfRotamers) {
	if(_libName == "") {
		_libName = defaultLibrary;
	}

	if (residueExists(_libName, _resName)) {
			if(libraries[_libName][_resName].internalCoor.size() > _startIdx) {
				if(_numOfRotamers == -1) {
					libraries[_libName][_resName].internalCoor.erase(libraries[_libName][_resName].internalCoor.begin() + _startIdx,libraries[_libName][_resName].internalCoor.end());
					libraries[_libName][_resName].rotamerBins.erase(libraries[_libName][_resName].rotamerBins.begin() + _startIdx,libraries[_libName][_resName].rotamerBins.end());
				} else {
					libraries[_libName][_resName].internalCoor.erase(libraries[_libName][_resName].internalCoor.begin() + _startIdx,libraries[_libName][_resName].internalCoor.begin() + _startIdx + _numOfRotamers);
					libraries[_libName][_resName].rotamerBins.erase(libraries[_libName][_resName].rotamerBins.begin() + _startIdx,libraries[_libName][_resName].rotamerBins.begin() + _startIdx + _numOfRotamers);
				}
				return true;
			} else {
				//cerr << "Warning 2320: Not enough conformers for the residue" << endl;
				return false;
			}
	} else {

		//cerr << "Warning 2321: Unable to find conformer to remove" << endl;
		return false;

	}	
}

bool RotamerLibrary::trimToLevel(string _libName, string _level) {
	if(_libName == "") {
		_libName = defaultLibrary;
	}
	if(libraryExists(_libName)) {
		set<string> resList = getResList(_libName);
		bool foundLevel = false;
		for(set<string>::iterator it = resList.begin(); it != resList.end(); it++) {
			unsigned levelRots = getLevel(_level,*it);
			if(levelRots > 0) {
				foundLevel = true;
				if(!removeRotamers(_libName,*it,levelRots)) {
					return false;
				} 
			}
		}
		if(!foundLevel) {
			return false;
		}
	} else {
		return false;
	}
	
	return true;
}

vector<string> RotamerLibrary::getInternalCoorDefinitionLines( string _libName,string _resName) {
	vector<string> lines;
	vector<InternalCoorDefi> icDefis;
	if(residueExists(_libName,_resName)) {
		icDefis = lastFoundRes->second.defi;
	} else {
		//cerr << "(Res,Lib) Not found: " << _resName << "," << _libName << endl; 
		return lines;
	}
	for (vector<InternalCoorDefi>::iterator i = icDefis.begin(); i != icDefis.end(); i++) {
		lines.push_back("DEFI ");
		switch ((*i).type) {
			case 0:
				//bond
				lines.back() +=  (*i).atomNames[0] + " " +(*i).atomNames[1];
				break;
			case 1:
				// angle
				lines.back() += (*i).atomNames[0] + " " + (*i).atomNames[1] + " " + (*i).atomNames[2];
				break;
			
			case 2:
				// Dihedral 
				lines.back() +=  (*i).atomNames[0] + " " + (*i).atomNames[1] + " " + (*i).atomNames[2] + " " + (*i).atomNames[3];
				break;
                         
			case 3:
				// Improper 
				lines.back() +=  (*i).atomNames[0] + " " + (*i).atomNames[1] + " *" + (*i).atomNames[2] + " " + (*i).atomNames[3];
				break;
		}
	}
	return lines;
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
	} else if (D.type == 1) {
		D.atomNames.erase(D.atomNames.begin()+3, D.atomNames.end());
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
	 *       MOBI CB HA HB1 HB2 HB3
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
			vector<string> mobileAtoms = resItr->second.mobileAtoms;
			vector<InternalCoorDefi> defi = resItr->second.defi;
			for (unsigned int i=0; i<mobileAtoms.size(); i++) {
				// for each atom to rebuild
				vector<string> bond;
				vector<string> angle;
				vector<string> dihe;
				vector<unsigned int> indeces(3, 0);
				bool improper_flag = false;
				for (unsigned int j=0; j<defi.size(); j++) {
					// for each internal coor definition
					if (defi[j].type == 0 && defi[j].atomNames[1] == mobileAtoms[i]) {
						// found a bond where the last atom is this atom
						bond =  defi[j].atomNames;
						indeces[0] = j;
					} else if (defi[j].type == 1 && defi[j].atomNames[2] == mobileAtoms[i]) {
						// found an angle where the last atom is this atom
						angle =  defi[j].atomNames;
						indeces[1] = j;
					} else if ((defi[j].type == 2 || defi[j].type == 3) && defi[j].atomNames[3] == mobileAtoms[i]) {
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
					cerr << "WARNING 95123: bond not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << mobileAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
					return false;
				} else if (angle.size() == 0) {
					cerr << "WARNING 95128: angle not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << mobileAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				} else if (dihe.size() == 0) {
					cerr << "WARNING 95133: dihedral/improper not found for building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << mobileAtoms[i] << ", in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				} else if (bond[0] != dihe[2]) {
					cerr << "WARNING 95138: bond atoms do not match dihe/improper building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << mobileAtoms[i] << " (" << bond[0] << " " << bond[1] << " != " << dihe[0] << " " << dihe[1] << " " << dihe[2] << " " << dihe[3] << "), in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				}else if (angle[0] != dihe[1] || angle[1] != dihe[2]) {
					cerr << "WARNING 95143: angle atoms do not match dihe/improper building IC definition for library " << libItr->first << ", residue " << resItr->first << ", atom " << mobileAtoms[i] << " (" << angle[0] << " " << angle[1] << " " << angle[2] << " != " << dihe[0] << " " << dihe[1] << " " << dihe[2] << " " << dihe[3] << "), in bool RotamerLibrary::calculateBuildingICentries()" << endl;
				}
				resItr->second.buildingInstructions[mobileAtoms[i]].atomNames = dihe;
				resItr->second.buildingInstructions[mobileAtoms[i]].defiIndeces = indeces;
				resItr->second.buildingInstructions[mobileAtoms[i]].improper_flag = improper_flag;
			}
		}
	}

	


	return true;
}
string RotamerLibrary::getInitAtomsLine( string _libName, string _resName) {
	// DEPRECATED 
	cerr << "WARNING: using deprecated function string RotamerLibrary::getInitAtomsLine( string _libName, string _resName), use getMobileAtomsLine instead" << endl;
	return getMobileAtomsLine(_libName, _resName);
}
string RotamerLibrary::getMobileAtomsLine( string _libName, string _resName) {
	
	string line = "MOBI";
	vector<string> iAtoms;
	if(residueExists(_libName,_resName)) {
		iAtoms = lastFoundRes->second.mobileAtoms;
	} else {
		//cerr << "(Res,Lib) Not found: " << _resName << "," << _libName << endl; 
		return line;

	}
	for (int i = 0; i < iAtoms.size(); i++) {
		line += (" " + iAtoms[i]);
	}
	return line;
}

string RotamerLibrary::toString(){
	string result;

	// for each library
	for (map<string, map<string, Res> >::iterator libItr=libraries.begin(); libItr!=libraries.end(); libItr++) {
		result += "LIBRARY "+libItr->first+"\n\n";
		//result += "CHARMMPAR 22 27\n\n";
		// for each residue
		for (map<string, Res>::iterator resItr=libItr->second.begin(); resItr!=libItr->second.end(); resItr++) {
			result += "RESI "+resItr->first+"\n";
			// Print mobile atoms
			vector<string> mobileAtoms = resItr->second.mobileAtoms;
			result += "MOBI ";
			for (unsigned int i=0; i<mobileAtoms.size(); i++) {
				result +=  mobileAtoms[i] + " ";
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
					char c[100];
					sprintf(c,"%11.5f",iCoor[i][j]);
					result += (string)c + " ";
				}
				result += "\n";
			}
			
			
		}
	}
	return result;
	
}
vector<string> RotamerLibrary::getAllInternalCoorLines(string _libName, string _resName) {
	vector<string> lines;
	vector<vector<double> > coor;
	if(residueExists(_libName, _resName)) {
		coor = libraries[_libName][_resName].internalCoor;
	} else {
		//cerr << "(Res,Lib) Not found: " << _resName << "," << _libName << endl; 
		return lines;
	}

	for(int i=0; i< coor.size(); i++) {
		lines.push_back("CONF ");
		for (int j = 0; j < coor[i].size(); j++) {
			char coord[1000];
			sprintf(coord,"%11.5f",coor[i][j]);
			lines.back() += coord;
		}
		if(libraries[_libName][_resName].rotamerBins[i] != 0) {
			char bin[12];
			sprintf(bin,"%11u",libraries[_libName][_resName].rotamerBins[i]);
			lines.back() += bin;
		}
		//lines += "\n";
	}
	return lines;
}

string RotamerLibrary::getInternalCoorLine(string _libName, string _resName, unsigned int _num) {
	string line = "";
	vector<double> coor;
	if(residueExists(_libName, _resName) && _num < libraries[_libName][_resName].internalCoor.size()) {
		coor = libraries[_libName][_resName].internalCoor[_num];
	} else {
		//cerr << "(Res,Lib) Not found: " << _resName << "," << _libName << endl; 
		return line;
	}

	for (int j = 0; j < coor.size(); j++) {
		char coord[1000];
		sprintf(coord,"%11.5f",coor[j]);
		line += coord;
	}
	if(libraries[_libName][_resName].rotamerBins[_num] != 0) {
		char bin[12];
		sprintf(bin,"%11u",libraries[_libName][_resName].rotamerBins[_num]);
		line += bin;
	}

	return line;
}


bool RotamerLibrary::readFile(string _filename, bool _append) {
	if( _append == false) {
		reset();  // Remove this??
	}
	
	if (!rotReader->open(_filename) || !rotReader->read()) { 
		return false;
	 }
	rotReader->close();
	return true;
}

bool RotamerLibrary::writeFile(string _filename) {if (!rotWriter->open(_filename)) return false; bool result = rotWriter->write(this); rotWriter->close();return result;}


