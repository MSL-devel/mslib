#include "System.h"
#include "FormatConverter.h"

using namespace MSL;
using namespace std;

int main(int argc, char* argv[]) {
	if(argc != 3) {
		cerr << "usage: convertToPdbNames <inp> <out>" << endl;
		exit(0);
	}

	System sys;
	if(!sys.readPdb(string(argv[1]))) {
		cerr << "Unable to read " << argv[1] << endl;
		exit(0);
	}
	
	FormatConverter fc;

	fc.setPdbFromCharmm(&sys,"22");

	if(!sys.writePdb(string(argv[2]))) {
		cerr << "Unable to write " << argv[2] << endl;
		exit(0);
	}


}


