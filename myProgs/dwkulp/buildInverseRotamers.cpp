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
#include "buildInverseRotamers.h"

using namespace std;
using namespace MSL;

// Build Inverse Rotamers from input PDB and collect the bounding boxes
void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers, vector<map<string,double> > &_boundingBoxes);


// MslOut 
static MslOut MSLOUT("buildInverseRotamers");

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

	// For each residue, build inverse rotamers
	vector<AtomContainer *> inverseRotamers;
	vector<map<string,double> > boundingBoxes;

	// Build inverse rotamers, get bounding boxes
	buildInverseRotamers(opt,sys,inverseRotamers,boundingBoxes);
	
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
		cout << "buildInverseRotamers --pdb PDB --rotlib ROTLIB\n";
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
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}

	opt.modelOut = OP.getBool("multiModelPdb");
	if (OP.fail()){
	  opt.modelOut = false;
	}
	
	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  vector<string> nameToks = MslTools::tokenize(opt.pdb,".");
	  opt.outPdb = nameToks[0]; 
	  cerr << "WARNING outPdb is defaulted to "<<opt.outPdb<<endl;
	}
	return opt;
}

void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers, vector<map<string,double> > &_boundingBoxes){
  
	// Create a pdb topology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(_opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);

	// Each position build opt.numRotamer rotamers
	for (uint i = 0; i < _sys.positionSize();i++){

	  cout << "I: "<<i<<endl;
	  Residue &res = _sys.getResidue(i);
	  if (res.getResidueName() == "HIS"){ res.setResidueName("HSD");}
	  string id = res.getIdentityId();
	 

	  // Build using inverse-rotamer IC building 
	  //AtomContainer ac = pdbTop.getResidue(id,res.getAtomPointers(),_opt.numRotamers);
	  
	  // Build a generic N-CA-C, then build rotmaers onto it
	  MSLOUT.stream() << "BUILD GENERICS"<<endl;
	  AtomContainer ac = pdbTop.getGenericResidue(id,_opt.numRotamers);
	  MSLOUT.stream() << "Done getting generic residues"<<endl;
	    

	  _inverseRotamers.push_back(new AtomContainer(ac.getAtomPointers()));
	  AtomContainer *lastRes = _inverseRotamers.back();

	  string outFile = _opt.outPdb;
	  if (_opt.modelOut){
	    outFile = "MODELS_"+_opt.outPdb;
	  }
 	  // write it out.
	  char name[80];
	  sprintf(name,"ROTAMERS_%s-%04d.pdb",outFile.c_str(),i);
	  PDBWriter pout;

	  pout.open((string)name);
	  for (uint j = 0; j < (*lastRes)(0).getNumberOfAltConformations();j++){

	    // Change each atoms active conformation
	    for (uint k =0;k <(*lastRes).size();k++){
	      (*lastRes)(k).setActiveConformation(j);
	      (*lastRes)(k).setResidueNumber((*lastRes)(k).getResidueNumber()+1);

	    }

	    // Need smartRmsdAlignment to take a filter for atom types to align.
	    Transforms t;
	    if (!t.smartRmsdAlignment((*lastRes).getAtomPointers(),res.getAtomPointers(),Transforms::MT_ATOMNAME)){
	      MSLOUT.stream() << "ERROR WITH ALIGNING lastRes and res"<<endl;
	    }

	    AtomSelection sel((*lastRes).getAtomPointers());
	    //AtomPointerVector vec = sel.select("not name N+C+O+HA");
	    AtomPointerVector vec = sel.select("not name HA");
	    pout.write(vec,false,true,_opt.modelOut);


	  }
	}
}

