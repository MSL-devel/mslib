#include "AtomPointerVector.h"
#include "AtomBondBuilder.h"
#include "PDBReader.h"
#include "System.h"

using namespace std;

using namespace MSL;


int main(int argc,char *argv[]) {

	/*
	PDBResysr rAv;
	rAv.open(argv[1]);
	rAv.read();
	*/
	if (argc == 1) {
		cout << "Please specify a PDB file" << endl;
		exit(0);
	}
	System sys;
	if (!sys.readPdb(argv[1])) {
		cerr << "Cannot read pdb file " << argv[1] << endl;
		exit(1);
	}
	AtomPointerVector av = sys.getAtoms();
	//cout << av;

	AtomBondBuilder abb;
	abb.buildConnections(av);

	for (unsigned int i=0; i<av.size(); i++) {
		cout << *av[i] << endl;
		vector<Atom*> bonded = av[i]->getBonds();
		for (unsigned int j=0; j<bonded.size(); j++) {
			cout << "   " << *bonded[j] << endl;
		}
	}


}

