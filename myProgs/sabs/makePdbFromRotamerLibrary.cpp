#include "MslTools.h"
#include "SystemRotamerLoader.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"
#include <set>
#include <vector>

using namespace MSL;
using namespace std;

int main(int argc, char* argv[]) {
	if(argc != 4 ) {
		cout << "Usage: makePdbFromRotamerLibrary <RotamerLib> <number of rotamers> <outputFile> " << endl;
		exit(0);
	}
	
	string inputFile = string(argv[1]);
	int numRots = MslTools::toInt(string(argv[2]));

	PDBWriter pWriter;
	if(!pWriter.open(string(argv[3]))) {
		cerr << "unable to open " << argv[3] << endl;
		exit(0);
	}

	System sys;

	string seq = "A: ";
	
	SystemRotamerLoader rotLoader;
	rotLoader.setSystem(sys);

	if(!rotLoader.readRotamerLibraryFile(inputFile)) {
		cerr << "Unable to read " << inputFile << endl;
		exit(0);
	}

	RotamerLibrary *rotlib = rotLoader.getRotamerLibrary();
	set<string> resList = rotlib->getAllResList();

	for(set<string>::iterator it = resList.begin(); it != resList.end(); it++) {
		seq += *it + " ";
	}
	
	cout << seq << endl;

	CharmmSystemBuilder csb(sys,"/data00/sabs/pdb/top_all22_prot.inp","/data00/sabs/pdb/par_all22_prot.inp");
	if(!csb.buildSystem(seq)) {
		cerr << "Unable to build System" << endl;
		exit(0);
	}
	cout << sys << endl;
	if(!sys.seed("A 1 C","A 1 CA","A 1 N")) {
		cerr << "Unable to seed " << endl;
		exit(0);
	}
	sys.buildAllAtoms();
	cout << sys.getAtomPointers() << endl;
	vector<Position*> positions = sys.getPositions();
	for(vector<Position*>::iterator it = positions.begin(); it != positions.end(); it++) {
		string resName = (*it)->getResidueName();
		if(resName != "VAL") {
			continue;
		}
			
		int tempSize = numRots;
		if(rotlib->size("",resName) < numRots) {
			tempSize = rotlib->size("",resName);
		}
		if(tempSize == 0) {
			continue;
		}
		if(!rotLoader.loadRotamers(*it,resName,tempSize)) {
			cerr << "Unable to load rotamers " << resName << endl;
			continue;
		}
		cout << resName << " " << rotlib->size("",resName) << endl;
		for(int i = 0; i < tempSize; i++) {
			((*it)->getCurrentIdentity()).setActiveConformation(i);
			AtomPointerVector& atoms = (*it)->getAtomPointers();
			PDBWriter w;
			if(!w.open(string(argv[3]) + "_" + MslTools::intToString(i))) {
				cerr << "Unable to write " << i << endl;
				exit(0);
			}
			if(!w.write(atoms)) {
				cerr << "unable to write " << endl;
				exit(0);
			}
			w.close();
			//cout << atoms << endl;
			if(!pWriter.write(atoms,true,false,true)) {
				cerr << "unable to write " << endl;
				exit(0);
			}
		}
	}
	pWriter.close();

}
