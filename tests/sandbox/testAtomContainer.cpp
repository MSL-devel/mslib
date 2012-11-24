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
#include "AtomContainer.h"
#include "Transforms.h"
#include "MslTools.h"

using namespace MSL;
using namespace std;

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cout << "USAGE:\nexample_AtomContainer_usage <path_of_exampleFiles_directory> <atom identifier>\n\nAtom identifier: either A,37,CB or A,37,ILE,CB (with insertion code A,37A,CB)" << endl;
		exit(0);
	}

	AtomContainer molecule;
	molecule.setAddAtomsAsAltCoors(true);
	molecule.readPdb(argv[1]);
	cout << molecule;
	cout << "The PDB contains " << molecule.size() << " atoms" << endl;
	cout << endl;
	cout << "Check the alt coordinates of the first atom" << endl;

	for (unsigned int i=0; i<molecule[0].getNumberOfAltConformations(); i++) {
		molecule[0].setActiveConformation(i);
		cout << molecule[0] << endl;
	}
	cout << "Set the molecule in the second alt conformation" << endl;
	if (!molecule.setActiveConformation(1)) {
		cerr << "WARNING: at least one atom did not have 2 alt confs" << endl;
	}
	cout << molecule << endl;
	exit(0);

	if (argc < 3) {
		cout << endl;
		cout << "MISSING ATOM IDENTIFIER!\n\nUSAGE:\nexample_AtomContainer_usage <path_of_exampleFiles_directory> <atom identifier>\n\nAtom identifier: either A,37,CB or A,37,ILE,CB (with insertion code A,37A,CB)" << endl;
		exit(0);
	} else {
		if (molecule.atomExists(argv[2])) {
			cout << "Atom " << argv[2] << " exists" << endl;
			cout << "Print it using getLastFoundAtom()" << endl;
			cout << "   " << molecule.getLastFoundAtom() << endl;
			cout << "Print it using molecule(" << argv[2] << ")" << endl;
			cout << "   " << molecule(argv[2]) << endl;
		} else {
			cout << "Atom " << argv[2] << " does NOT exists" << endl;
		}
	}

	return 0;
}

