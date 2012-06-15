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
#include "PDBTopologyBuilder.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to change the active identity of
 *  a system that has multiple identities at one position.
 * 
 *  It reads a non-canonical PDB file that has the coordinates of two
 *  different amino acid types at position 2 (ILE and LEU)
 *
 *  The program prints the sequence, writes a PDB with the ILE as the
 *  active identity, then it switches the position to LEU and writes
 *  a second PDB
 * 
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 3) {
		cerr << "USAGE:\nexample_add_identity_to_position <path_of_exampleFiles_directory> <path_of_toppar_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "    Example on adding a second residue identity at a position (" << MslTools::getMSLversion() << ")  " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;
	
	cout << "Read example0001.pdb, a 3 amino acid PDB file with ILE at position 2" << endl;
	string file = "example0001.pdb";
	file = (string)argv[1] + "/" + file;

	// read the PDB into a new System
	System sys;
	if (!sys.readPdb(file)) {
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}

	// print the sequence
	cout << endl;
	cout << "Print the sequence" << endl;
	cout << sys << endl;
	cout << endl;
	cout << "=============================" << endl;

	// Get the variable position (as a reference) using chain and residue number and print it
	Position & pos = sys.getPosition("A,2");

	cout << "The default identity at position 2 is the first one (marked with *) when the position is printed" << endl;
	cout << pos << endl;

	string file2 = "top_pdb2.3_noH.inp";
	file2 = (string)argv[2] + "/" + file2;
	PDBTopologyBuilder builder(sys, file2);

	if (builder.fail()) {
		cerr << "Error reading topology file top_pdb2.3.inp, please speficy the path of the \"toppar\" directory as an argument" << endl;
		exit(1);
	}
	builder.addIdentity("A,2", "LEU", "N H CA C O");

	cout << "The default identity at position 2 is the first one (marked with *) when the position is printed" << endl;
	cout << pos << endl;

	// write a PDB file with ILE as the active identity at position 2
	if (!sys.writePdb("/tmp/example0001_ile.pdb")) {
		cerr << "Error writing PDB /tmp/example0000_ile.pdb" << endl;
		exit(1);
	}
	cout << endl;
	cout << "Written a PDB file with ILE as the active identity at position 2: /tmp/example0001_ile.pdb" << endl;
	cout << endl;
	cout << "=============================" << endl;

	// change the active identity to LEU
	pos.setActiveIdentity("LEU");
	cout << "Set the active identity to LEU at position 2" << endl;
	cout << pos << endl;

	// write a PDB file with LEU as the active identity at position 2
	if (!sys.writePdb("/tmp/example0001_leu.pdb")) {
		cerr << "Error writing PDB /tmp/example0000_leu.pdb" << endl;
		exit(1);
	}
	cout << "Written a PDB file with ILE as the active identity at position 2: /tmp/example0001_leu.pdb" << endl;
	cout << endl;


	return 0;
}
