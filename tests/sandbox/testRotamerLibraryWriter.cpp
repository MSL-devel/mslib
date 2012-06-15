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

#include <iostream>

#include "RotamerLibraryReader.h"
#include "RotamerLibraryWriter.h"

using namespace std;

using namespace MSL;


int main(int argc, char* argv[]) {

	if(argc != 3) {
		cerr << "Usage: testRotamerLibraryWriter <rotlib file> <outputFile> " << endl;
		exit(0);

	}

	RotamerLibrary rotlib;
	if(!rotlib.readFile(argv[1])) {
		cerr << "Unable to read " << argv[1] << endl;
		exit(0);
	}
	if(!rotlib.writeFile(argv[2])) {
		cerr << "Unable to write " << argv[2] << endl;
		exit(0);
	}
//	RotamerLibraryReader rotRead("/library/rotlib/balanced/rotlib-balanced.txt", rotlib);

	cout << "Read " << rotlib.getNumberOfLibraries() << " rotamer libraries" << endl;
/*	RotamerLibraryWriter rotWrite;
	rotWrite.open(argv[2]);
	rotWrite.write(&rotlib);
	rotWrite.close();
	
//	rotlib.reset();
	RotamerLibraryReader rotRead1("/library/rotlib/balanced/rotlib-balanced-200.txt", rotlib);
	rotRead1.close();
	rotWrite.open("rotWrite1.txt");
	rotWrite.write(rotlib);
	rotWrite.close();
	rotlib.removeAllConformations();
	rotWrite.open("rotWrite-clear.txt");
	rotWrite.write(rotlib);
	rotWrite.close();
*/	
}
/*

	vector<string> residues;

	residues.push_back("ALA");
	residues.push_back("LEU");
	residues.push_back("PIZZA");
	string lib = "BALANCED-200";
	cout << "Library " << lib;
	if (rotlib.libraryExists(lib)) {
		cout << " exists" << endl;
	} else {
		cout << " does not exists" << endl;
	}
	for (unsigned int i=0; i<residues.size(); i++) {
		cout << "    Residue " << residues[i];
		if (rotlib.residueExists(lib, residues[i])) {
			cout << " exists" << endl;
		} else {
			cout << " does not exists" << endl;
		}
	}
	lib = "PIZZA";
	cout << "Library " << lib;
	if (rotlib.libraryExists(lib)) {
		cout << " exists" << endl;
	} else {
		cout << " does not exists" << endl;
	}
	for (unsigned int i=0; i<residues.size(); i++) {
		cout << "    Residue " << residues[i];
		if (rotlib.residueExists(lib, residues[i])) {
			cout << " exists" << endl;
		} else {
			cout << " does not exists" << endl;
		}
	}
	rotRead.open("/tmp/pizza.txt");
	rotRead.read();
	rotRead.close();
	cout << "We have read " << rotlib.getNumberOfLibraries() << " rotamer libraries" << endl;
	lib = "BALANCED-200";
	cout << "Library " << lib;
	if (rotlib.libraryExists(lib)) {
		cout << " exists" << endl;
	} else {
		cout << " does not exists" << endl;
	}
	for (unsigned int i=0; i<residues.size(); i++) {
		cout << "    Residue " << residues[i];
		if (rotlib.residueExists(lib, residues[i])) {
			cout << " exists" << endl;
		} else {
			cout << " does not exists" << endl;
		}
	}
	lib = "PIZZA";
	cout << "Library " << lib;
	if (rotlib.libraryExists(lib)) {
		cout << " exists" << endl;
	} else {
		cout << " does not exists" << endl;
	}
	for (unsigned int i=0; i<residues.size(); i++) {
		cout << "    Residue " << residues[i];
		if (rotlib.residueExists(lib, residues[i])) {
			cout << " exists" << endl;
		} else {
			cout << " does not exists" << endl;
		}
	}

	vector<string> init = rotlib.getMobileAtoms("BALANCED-200", "VAL");
	cout << "BALANCED-200/VAL has " << init.size() << " init atoms" << endl;
	for (vector<string>::iterator k=init.begin(); k<init.end(); k++) {
		cout << *k << endl;
	}
	vector<unsigned int> type;
	vector<vector<string> > atomNames;
	vector<vector<int> > resnumCorrectors;
	vector<vector<double> > intCoor = rotlib.getInternalCoor("BALANCED-200", "VAL");
	rotlib.getInternalCoorDefinition("BALANCED-200", "VAL", type, atomNames, resnumCorrectors);
	cout << "BALANCED-200/VAL has " << type.size() << " types. " << atomNames.size() << " atoms names and " << resnumCorrectors.size() << " resnum correctors" << endl;
	for (unsigned int i=0; i<type.size(); i++) {
		cout << i << " type=" << type[i];
		for (unsigned int j=0; j<atomNames[i].size(); j++) {
			cout << " " << atomNames[i][j] << " (" << resnumCorrectors[i][j] << ")";
		}
		cout << endl;
	}
	for (unsigned int i=0; i<intCoor.size(); i++) {
		for (unsigned int j=0; j<intCoor[i].size(); j++) {
			cout << " " << intCoor[i][j] << " ";
		}
		cout << endl;
	}

*/
		/*
		vector<string> getMobileAtoms(string _libName, string _resName);
		void getInternalCoorDefinition(string _libName, string _resName, vector<unsigned int> & _type, vector<vector<string> > & _atomNames, vector<vector<int> > & _resnumCorrectors);
		vector<vector<double> > getInternalCoor(string _libName, string _resName);
		*/

/*
	return 0;
}
*/
