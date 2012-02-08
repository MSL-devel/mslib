#include "System.h"
#include "PDBReader.h"
#include "FormatConverter.h"

using namespace MSL;
using namespace std;

int main(int argc, char* argv[]) {
	
	if(argc < 4) {
		cerr << "Usage: testFormatConverter <inputFile> <outputFile> <conversion(p2c | c2p)> [<charmmVersion=22>] " << endl;
		exit(0);
	}

	string version = "22";
	if(argc == 5) {
		version = string(argv[4]);
	}
	string conv = MslTools::toUpper(argv[3]);
	PDBReader pRead;
	pRead.open(string(argv[1]));
	if(!pRead.read()) {
		cerr << "Unable to read " << argv[1] << endl;
		exit(0);
	}

	System sys;
	sys.addAtoms(pRead.getAtomPointers());

	FormatConverter fc(&sys);
	if(conv == "P2C") {
		cout << "Converting from Pdb to Charmm Names" << endl;
		fc.setCharmmFromPdb(version);
	} else {
		cout << "Converting from Charmm to Pdb Names" << endl;
		fc.setPdbFromCharmm(version);
	}
	sys.setActiveModel(2);
	if(sys.writePdb(string(argv[2]))) {
		cout << "Output written to " << argv[2] << endl;
	} else {
		cout << "Unable to write Output to " << argv[2] << endl;
	}
}

