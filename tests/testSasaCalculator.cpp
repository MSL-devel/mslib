#include <iostream>

#include "System.h"
#include "SasaCalculator.h"

using namespace std;

using namespace MSL;


int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\ntestSasaCalculator <file.pdb>" << endl;
		exit(0);
	}

	System sys;
	if (!sys.readPdb(argv[1])) {
		cerr << "Cannot read pdb " << argv[0] << endl;
	}

	SasaCalculator sas(sys.getAtoms());
	sas.calcSasa();
	sas.printSasaTable();
	cout << endl;
	sas.printResidueSasaTable();

	

	return 0;
}

