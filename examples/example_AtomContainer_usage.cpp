#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "MslTools.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to use the AtomContainer
 * 
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_AtomContainer_usage <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     Example on the usage of AtomContainer (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	cout << "Read example0001.pdb into the AtomContainer" << endl;
	string file = "example0001.pdb";
	file = (string)argv[1] + "/" + file;

	// read the PDB into a new AtomContainer
	AtomContainer atoms;
	if (!atoms.readPdb(file)) {
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	// print all atoms
	cout << "Print the AtomContainer" << endl;
	cout << atoms << endl;
	cout << endl;
	cout << "=============================" << endl;

	// cycles among atoms with the [] operator
	cout << "Cycle over the atoms of the AtomContainer" << endl;
	for (unsigned int i=0; i<atoms.size(); i++) {
		cout << "Atom " << i << " is " << atoms[i] << endl;
	}
	cout << endl;

	cout << endl;
	cout << "=============================" << endl;
	// check that a specific atom (chain resnum and atom name) exists and get it
	if (atoms.exists("A 2 CD1")) {
		cout << "Atom A 2 CD1 exists" << endl;
		Atom & a = atoms.getLastFoundAtom();
		cout << a << endl;
	} else {
		cout << "Atom A 2 CD1 does no exist" << endl;
	}

	cout << endl;
	cout << "=============================" << endl;
	// use the () operator to get an atom that surely exists (or you'll get an error
	cout << "Get the same atom with the () operator" << endl;
	cout << atoms("A 2 CD1") << endl;

	cout << endl;
	cout << "=============================" << endl;
	// let's do something with these atoms: pass all the atoms to a Transform object, apply a tranlation and write a new PDB
	Transforms T;
	T.Xrotate(atoms.getAtoms(), 90);

	cout << "Apply a 90 degree rotation to all atoms over the X axis using the Transform object and write a new PDB file" << endl;
	cout << endl;
	if (atoms.writePdb("/tmp/example0001_rotated.pdb")) {
		cout << "/tmp/example0001_rotated.pdb written" << endl;
	} else {
		cout << "Error writing /tmp/example0001_rotated.pdb" << endl;
	}

	return 0;
}
