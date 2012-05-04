#include "System.h"
#include "PDBReader.h"
#include "FormatConverter.h"

using namespace MSL;
using namespace std;

int main(int argc, char* argv[]) {
	
	if(argc < 5) {
		cerr << "Usage: testFormatConverter <inputFile> <outputFile> <origVersion> <targetVersion>" << endl;
		exit(0);
	}

	PDBReader pRead;
	pRead.open(string(argv[1]));
	if(!pRead.read()) {
		cerr << "Unable to read " << argv[1] << endl;
		exit(0);
	}


	FormatConverter fc;
	if(!fc.setNamespaces(string(argv[3]),string(argv[4]))) {
			cerr << "Unable to convert from " << argv[3] << " to " << argv[4] << endl;
			exit(0);
	}
	fc.convert(pRead.getAtomPointers());

	PDBWriter pWrite;
	pWrite.open(argv[2]);
	if(pWrite.write(pRead.getAtomPointers())) {
		cout << "Output written to " << argv[2] << endl;
	} else {
		cout << "Unable to write Output to " << argv[2] << endl;
	}
}

