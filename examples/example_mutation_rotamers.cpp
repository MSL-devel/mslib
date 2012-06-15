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
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"

using namespace std;
using namespace MSL;

// MslOut can suppress output, to make example output clean
static MslOut MSLOUT("example");

int main(int argc, char *argv[]) {

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnAllOff();

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 3) {
		cerr << "USAGE:\nexample_mutation_rotamers <path_of_exampleFiles_directory> <full path of rotamer library file>" << endl;
		cerr << "\n\n./bin/example_mutation_rotamers exampleFiles /Users/dwkulp/software/mslib/library/rotlib/balanced/rotlib-balanced-200.txt\n\n";
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to mutate a position to any identity and get multiple conformations (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string pdbFile = "example0005.pdb";
	pdbFile = (string)argv[1] + "/" + pdbFile;
	cout << "Create an AtomContainer and read the atoms from " << pdbFile << endl;


	// Read pdb file into a System object
	System sys;
	sys.readPdb(pdbFile);

	// Get one of the residues : chain B position 4
	Residue &res = sys.getIdentity("B,1");

	// Create a PDBTopology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary((string)argv[2]);
	pdbTop.setAddAtomsFromRotLib(true);

	// Get backbone atoms  
	AtomPointerVector backboneAtoms = pdbTop.getBackboneAtoms(res);
	
	// *********  MUTATE *************** //
	//    B,2,ASP to B,2,ASN and get 3 ASN rotamers
	AtomContainer newAtoms = pdbTop.getResidue("B,1,ASN",backboneAtoms,3);

	// Add atoms to the system
	sys.getPosition("B,1").addIdentity(newAtoms.getAtomPointers(),"ASN");
	sys.getPosition("B,1").setActiveIdentity("ASN");

	// Write out new pdb with an Asn at position 1, chain B.  instead of Asp.
	sys.writePdb("/tmp/example00005_withAsn.pdb");
	cout << "Wrote file /tmp/example00005_withAsn.pdb , has Asn instead of Asp\n";

	// Change conformation , First conformation is 2.  Write out conformation 0 and 1.
	sys.getResidue("B,1,ASN").setActiveConformation(0);
	sys.writePdb("/tmp/example00005_withAsn_rotamer0.pdb");
	cout << "Wrote file /tmp/example00005_withAsn_rotamer0.pdb , has Asn (alternate rotamer) instead of Asp\n";

	sys.getResidue("B,1,ASN").setActiveConformation(1);
	sys.writePdb("/tmp/example00005_withAsn_rotamer1.pdb");
	cout << "Wrote file /tmp/example00005_withAsn_rotamer1.pdb , has Asn (alternate rotamer) instead of Asp\n";

	
}


