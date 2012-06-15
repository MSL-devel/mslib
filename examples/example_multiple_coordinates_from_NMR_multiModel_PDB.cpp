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

#include "System.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to add alternative coordinates to
 *  an Atom and how to switch the active one
 * 
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_multiple_coordinates_from_NMR_multiModel_PDB <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     Example on reading multi-model NMR style PDBs and changing the active model (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	cout << "Read example0007.pdb into the System" << endl;
	string file = "example0007.pdb";
	file = (string)argv[1] + "/" + file;

	System sys;
	sys.readPdb(file); // example0007.pdb contains 3 NMR style models

	cout << "=============================" << endl;
	cout << endl;
	cout << "The file example0007.pdb contains " << sys.getNumberOfModels() << " NMR style models. Let's write them into independent files" << endl;
	cout << endl;
	for (unsigned int i=0; i<sys.getNumberOfModels(); i++) {
		sys.setActiveModel(i);
		char c[1000];
		sprintf(c, "/tmp/example0007_model_%02u.pdb", i);
		sys.writePdb(c);
		cout << " - Written model " << i << " to pdb file " << c << endl;
	}
	cout << endl;

	return 0;
}
