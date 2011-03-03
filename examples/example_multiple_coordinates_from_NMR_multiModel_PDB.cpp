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
