#include <iostream>
#include <cstdlib>

#include "System.h"
#include "MslTools.h"

using namespace std;

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
	if (argc < 2) {
		cerr << "USAGE:\nexample_multipleResidueIdentities <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "    Example on setting multiple residue identities (" << MslTools::getMSLversion() << ")  " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;
	
	cout << "Read example0000.pdb, a non-canonical PDB file that has two different identities (ILE and LEU) at position 2" << endl;
	string file = "example0000.pdb";
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
	Position & pos = sys.getPosition("A", 2);

	cout << "The default identity at position 2 is the first one (marked with *) when the position is printed" << endl;
	cout << pos << endl;

	// write a PDB file with ILE as the active identity at position 2
	if (!sys.writePdb("/tmp/example0000_ile.pdb")) {
		cerr << "Error writing PDB /tmp/example0000_ile.pdb" << endl;
		exit(1);
	}
	cout << endl;
	cout << "Written a PDB file with ILE as the active identity at position 2: /tmp/example0000_ile.pdb" << endl;
	cout << endl;
	cout << "=============================" << endl;

	// change the active identity to LEU
	pos.setActiveIdentity("LEU");
	cout << "Set the active identity to LEU at position 2" << endl;
	cout << pos << endl;

	// write a PDB file with LEU as the active identity at position 2
	if (!sys.writePdb("/tmp/example0000_leu.pdb")) {
		cerr << "Error writing PDB /tmp/example0000_leu.pdb" << endl;
		exit(1);
	}
	cout << "Written a PDB file with ILE as the active identity at position 2: /tmp/example0000_leu.pdb" << endl;
	cout << endl;


	return 0;
}
