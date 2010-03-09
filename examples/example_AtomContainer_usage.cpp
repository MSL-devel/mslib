/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
	AtomContainer mol;
	if (!mol.readPdb(file)) {
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	// print all atoms
	cout << "Print the AtomContainer" << endl;
	cout << mol << endl;
	cout << endl;
	cout << "=============================" << endl;

	// cycles among atoms with the [] operator
	cout << "Cycle over the atoms of the AtomContainer" << endl;
	for (unsigned int i=0; i<mol.size(); i++) {
		cout << "Atom " << i << " is " << mol[i] << endl;
	}
	cout << endl;

	cout << endl;
	cout << "=============================" << endl;
	// check that a specific atom (chain resnum and atom name) exists and get it
	cout << "Checking existange and getting an atom" << endl;
	if (mol.atomExists("A,2,CD1")) {
		cout << "Atom A 2 CD1 exists" << endl;
		Atom & a = mol.getLastFoundAtom();
		cout << a << endl;
	} else {
		cout << "Atom A 2 CD1 does no exist" << endl;
	}
	cout << endl;
	cout << "=============================" << endl;
	// or if you are really sure it exist
	cout << "Using the getAtom(\"A,2,CD1\") function" << endl;
	Atom & a = mol.getAtom("A,2,CD1");
	cout << a << endl;

	cout << endl;
	cout << "=============================" << endl;
	// use the () operator to get an atom that surely exists (or you'll get an error
	cout << "Get the same atom with the (\"A,2,CD1\") operator" << endl;
	cout << mol("A,2,CD1") << endl;

	cout << endl;
	cout << "=============================" << endl;
	// for a residue
	cout << "Using the getResidue(\"A,2\") function" << endl;
	AtomPointerVector av = mol.getResidue("A,2");
	cout << av << endl;

	cout << endl;
	cout << "=============================" << endl;
	// let's do something with these atoms: pass all the atoms to a Transform object, apply a tranlation and write a new PDB
	Transforms T;
	T.Xrotate(mol.getAtoms(), 90);

	cout << "Apply a 90 degree rotation to all atoms over the X axis using the Transform object and write a new PDB file" << endl;
	cout << endl;
	if (mol.writePdb("/tmp/example0001_rotated.pdb")) {
		cout << "/tmp/example0001_rotated.pdb written" << endl;
	} else {
		cout << "Error writing /tmp/example0001_rotated.pdb" << endl;
	}

	return 0;
}
