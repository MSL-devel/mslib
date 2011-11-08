#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "PDBWriter.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
//#include "VectorHashing.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "discoverMotif.h"
#include "VectorHashing.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("alignByFrame");

int main(int argc, char *argv[]) {

	std::vector<std::string> filelist;
	MslTools::readTextFile(filelist,opt.pdblist);
	
	AtomContainer ref;
	ref.readPdb(filelist[0]);

	AtomSelection refSel(ref.getAtomPointers());
	
	Frame refFrame;
	refFrame.computeFrameFromPCA(refSel.select(opt.sel));

	for (uint i = 1; i < filelist.size();i++){
		AtomContainer pdb;
		pdb.readPdb(filelist[i]);
		
		AtomSelection pdbSel(pdb.getAtomPointers());

		Frame pdbFrame;
		pdbFrame.computeFrameFromPCA(pdbSel.select(opt.sel));

		
		

	}

}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "discoverMotif --pdb PDB --rotlib ROTLIB\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.sel = OP.getString("sel");
	if (OP.fail()){
		cerr << "ERROR 1111 sel not specified.\n";
		exit(1111);
	}

	return opt;
}
