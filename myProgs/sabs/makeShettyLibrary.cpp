#include "PDBReader.h"
#include "System.h"
#include "RotamerLibrary.h"
#include "RotamerLibraryBuilder.h"
#include "CharmmSystemBuilder.h"
#include "AtomSelection.h"

#include <fstream>

using namespace std;
using namespace MSL;

int main(int argc, char* argv[]) {
	if(argc != 3 ) {
		cout << "Usage: makeShettyLibrary <shettyFile> <outputFile>" << endl;
		exit(0);
	}
	
	string shettyPDB = string(argv[1]);
	string outfilename = string(argv[2]);
	cout << "Write to " << outfilename << endl;
	
	System sys;
	//CharmmSystemBuilder csb(sys,"/data01/sabs/pdb/top_all22_prot_eef1.1.inp","/data01/sabs/pdb/par_all22_prot.inp");
	CharmmSystemBuilder csb(sys,"/data00/sabs/pdb/top_all22_prot.inp","/data00/sabs/pdb/par_all22_prot.inp");
	
	PDBReader pRead(shettyPDB);
	pRead.open();
	if(!pRead.read()) {
		cerr << "Unable to read Shetty PDB" << endl;
		exit(0);
	}
	PolymerSequence seq(pRead.getAtomPointers());
	int tot = 0;
	for(unsigned int i = 0; i < seq.size(); i++) { // all chains
		tot += seq.chainSize(i);
//		cout << seq.chainSize(i) << endl;
		for(unsigned int j = 0; j < seq.chainSize(i); j++) { // all positions
			
			if(seq.getPositionIdentity(i,j,0) == "HSD") {
				seq.addPositionIdentity(i,j,"HSE");
				seq.addPositionIdentity(i,j,"HSP");
			}
		}
	}
//	cout << argv[1] << " " << tot << endl;
//	exit(0);

	csb.setBuildNonBondedInteractions(false);
	if(!csb.buildSystem(seq)) {
		cout << "Unable to build " << seq << endl;
		return 0;
	}
	
	
	AtomPointerVector& atoms = pRead.getAtomPointers();
	string selection = "his, resn HSD";
	AtomSelection sel(atoms);	
	AtomPointerVector his = sel.select(selection);
	
	sys.assignCoordinates(pRead.getAtomPointers());
	for(int i = 0; i < his.size(); i++) {
		his[i]->setResidueName("HSE");
	}
	
	sys.assignCoordinates(pRead.getAtomPointers());
	for(int i = 0; i < his.size(); i++) {
		his[i]->setResidueName("HSP");
	}
	
	sys.assignCoordinates(pRead.getAtomPointers());
	RotamerLibrary rotlib;
	if(!rotlib.readFile("/library/rotlib/balanced/rotlib-balanced.txt")) {
		cerr << "Unable to read /library/rotlib/balanced/rotlib-balanced.txt " << endl;
		exit(0);
	}

	rotlib.removeAllConformations();

	string library = rotlib.getDefaultLibrary();

	RotamerLibraryBuilder rBuild(&rotlib);
	sys.buildAllAtoms();

	cout <<  sys.positionSize() << endl;
	for(unsigned i = 0; i < sys.positionSize(); i++) {
		cout << i << endl;
		Position & p = sys.getPosition(i);
		for(unsigned j = 0; j < p.residueSize();j++) {
			Residue & r = p.getIdentity(j);
			if(!rBuild.addRotamer(r,library)) {
				cerr << "Unable to add " << r.getResidueName() << endl;
				continue;
			}
		}
	}

	if(!rotlib.writeFile(outfilename)) {
		cerr << "Unable to write to " << outfilename << endl;
	}
}
