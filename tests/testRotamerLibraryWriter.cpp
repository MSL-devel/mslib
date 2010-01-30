
#include <iostream>

#include "RotamerLibraryReader.h"
#include "RotamerLibraryWriter.h"

using namespace std;

int main(int argc, char* argv[]) {

	if(argc != 3) {
		cerr << "Usage: testRotamerLibraryWriter <rotlib file> <outputFile> " << endl;
		exit(0);

	}

	RotamerLibrary rotlib;
//	RotamerLibraryReader rotRead("/library/rotlib/balanced/rotlib-balanced.txt", rotlib);
	RotamerLibraryReader rotRead(argv[1], rotlib);
	rotRead.close();

	cout << "Read " << rotlib.getNumberOfLibraries() << " rotamer libraries" << endl;
	RotamerLibraryWriter rotWrite;
	rotWrite.open(argv[2]);
	rotWrite.write(rotlib);
	rotWrite.close();
/*	
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

	vector<string> init = rotlib.getInitAtoms("BALANCED-200", "VAL");
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
		vector<string> getInitAtoms(string _libName, string _resName);
		void getInternalCoorDefinition(string _libName, string _resName, vector<unsigned int> & _type, vector<vector<string> > & _atomNames, vector<vector<int> > & _resnumCorrectors);
		vector<vector<double> > getInternalCoor(string _libName, string _resName);
		*/

/*
	return 0;
}
*/
