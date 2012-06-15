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
#include <cstdlib>

#include "AtomContainer.h"
#include "MslTools.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to read and write PDBs with the
 *  AtomContainer
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_read_write_PDBs_with_the_AtomContainer <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to read and write PDBs with the AtomContainer (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string file = "example0001.pdb";
	file = (string)argv[1] + "/" + file;
	cout << "Create an AtomContainer and read the atoms from " << file << endl;

	// read the PDB into a new AtomContainer "container"
	AtomContainer container;
	if (!container.readPdb(file)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	// all atoms can be conveniently be printed using the << operator
	cout << "Print the AtomContainer" << endl;
	cout << container << endl;
	cout << endl;
	cout << "=============================" << endl;

	// cycle among atoms with the [] operator
	cout << "Cycle over the atoms of the AtomContainer" << endl;
	for (unsigned int i=0; i<container.size(); i++) {
		cout << "Atom " << i << " is " << container[i] << endl;
	}
	cout << endl;

	cout << endl;
	cout << "=============================" << endl;

	// get an atom by identifier ("chain,resnum,name"): no error handling here, one needs to be sure that the atoms exists
	cout << "Using the getAtom(\"A,2,CD1\") function" << endl;
	Atom & a = container.getAtom("A,2,CD1");
	cout << a << endl;

	cout << endl;
	cout << "=============================" << endl;
	// or just the same with the () operator
	cout << "Get the same atom with the [\"A,2,CD1\"] operator" << endl;
	Atom & b = container["A,2,CD1"];
	cout << b << endl;

	cout << "=============================" << endl;
	// One can check if an atom exists and use the getLastFoundAtom to get it 
	cout << "Check if atom chain A, resnum 2 and atom name CD1 exists (\"A,2,CD1\") and print it" << endl;
	if (container.atomExists("A,2,CD1")) {
		cout << "Atom A 2 CD1 exists" << endl;
		// get the atom without a second lookup
		Atom & c = container.getLastFoundAtom();
		cout << c << endl;
	} else {
		cout << "Atom A 2 CD1 does no exist" << endl;
	}
	cout << endl;
	cout << endl;
	cout << "=============================" << endl;

	// write the coordinates to PDB
	cout << "Write the coordinates to file: /tmp/example0001_out.pdb" << endl;
	if (!container.writePdb("/tmp/example0001_out.pdb")) {
		cout << "Error writing pdb file" << endl;
	} else {
		cout << "OK" << endl;
	}

	return 0;
}
