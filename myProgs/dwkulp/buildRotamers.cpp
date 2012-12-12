#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "buildRotamers.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("buildRotamers");

// SysEnv
static SysEnv SYSENV;

int main(int argc, char *argv[]) {

        // MslOut can suppress output, to make example output clean
        //MSLOUT.turnAllOff();

	Options opt = setupOptions(argc,argv);

	Timer t;
	double start = t.getWallTime();

	// Read in the pdb with a set of functional groups defined
	System sys;
	sys.readPdb(opt.pdb);

	// Create a PDBTopology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);

	// Get one of the residues : chain B position 4
	Residue &res = sys.getIdentity(opt.position);

	// Get backbone atoms  
	AtomPointerVector backboneAtoms = pdbTop.getBackboneAtoms(res);

	for (uint i = 0 ; i < opt.identity.size(); i++){
	    string residueId = opt.position+","+opt.identity[i];
	    AtomContainer newAtoms = pdbTop.getResidue(residueId,backboneAtoms,opt.numRotamers);

	    for (uint j = 0; j < opt.numRotamers;j++){
	      
	      for (uint a = 0; a < newAtoms.size();a++){
		newAtoms(a).setActiveConformation(j);
	      }
	      string fname = MslTools::stringf("%s_%3s_%04d.pdb",MslTools::getFileName(opt.pdb).c_str(),opt.identity[i].c_str(),j);
	      newAtoms.writePdb(fname);
	    }


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
		cout << "buildRotamers --pdb PDB --rotlib ROTLIB --numRotamers NUM --identity ASN --identity ARG --position POSID\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
	  opt.rotlib =	SYSENV.getEnv("MSL_ROTLIB");
	  cout << "MSL_ROTLIB is set to: "<<opt.rotlib<<endl;
	}
	opt.numRotamers = OP.getInt("numRotamers");
	if (OP.fail()){
	  opt.numRotamers = 10;
	}
	opt.identity = OP.getMultiString("identity");
	if (OP.fail()){
	  cerr << "ERROR 1111 identity not defined\n";
	  exit(1111);
	}
	opt.position = OP.getString("position");
	if (OP.fail()){
	  cerr << "ERROR 1111 position not defined\n";
	  exit(1111);
	}

	return opt;
}
