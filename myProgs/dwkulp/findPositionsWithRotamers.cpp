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
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "findPositionsWithRotamers.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("findPositionsWithRotamers");

// SysEnv
static SysEnv SYSENV;

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc,argv);

	// Read in the pdb with epitope
	System scaffold;
	scaffold.readPdb(opt.scaffold);

	// Use a part of the scaffold structure for alignment, rest for finding positions
	System nonEpitopeScaffold;
	AtomPointerVector epitopeAts;
	if (opt.epitopeSelect != ""){
	  AtomSelection sel(scaffold.getAtomPointers());
	  epitopeAts = sel.select(opt.epitopeSelect);
	  nonEpitopeScaffold.addAtoms(sel.inverseSelect(opt.epitopeSelect));
	}

	nonEpitopeScaffold.writePdb("foo.pdb");
	exit(0);

	// Reference structure, holds epitope.
	System ref;
	ref.readPdb(opt.ref);

	// Use a part of the reference structure for alignment
	AtomPointerVector refAts;
	if (opt.refSelect != ""){
	  AtomSelection sel(ref.getAtomPointers());
	  refAts = sel.select(opt.refSelect);
	}

	// Make sure selections were the same size..
	if (epitopeAts.size() != refAts.size()){
	  MSLOUT.stream() << "ERROR epitope selection and reference selection do not contain same number of atoms!"<<endl;
	  MSLOUT.stream() << "ERROR epitope selection ("<<opt.epitopeSelect<<") got "<<epitopeAts.size()<<" atoms."<<endl;
	  MSLOUT.stream() << "ERROR ref     selection ("<<opt.refSelect     <<") got "<<refAts.size()<<" atoms."<<endl;
	  exit(2222);
	}

	// Align 
	Transforms tr;
	if (!tr.rmsdAlignment(epitopeAts,refAts,scaffold.getAtomPointers())){
	  MSLOUT.stream() << "ERROR aligning epitope("<<opt.epitopeSelect") and ref("<<opt.refSelect<<")"<<endl;
	  exit(3333);
	}

	// For debugging store RMSD of initial alignment
	double alignRMSD = epitopeAts.rmsd(refAts);
	
	// For each residue, build inverse rotamers
	vector<AtomContainer *> inverseRotamers;

	// Build inverse rotamers, get bounding boxes
	buildInverseRotamers(opt,ref,inverseRotamers);


	// RMS check..
	for (uint i = 0; i < inverseRotamers.size();i++){

	  AtomPointerVector inverseCaCb;
	  inverseCaCb.push_back(inverseRotamers[i]("CA"));
	  inverseCaCb.push_back(inverseRotamers[i]("CB"));

	  int minAltConf = inverseRotamers[i]("CA").getNumberOfAltConformations();
	  for (int c = 0; c < minAltConf;c++){
	    inverseRotamers[i]("CA").setActiveConformation(c);
	    inverseRotamers[i]("CB").setActiveConformation(c);

	    
	    for (uint p = 0; p < nonEpitopeScaffold.positionSize();p++){
	      AtomPointerVector scaffoldCaCb;
	      scaffoldCaCb.push_back(nonEpitopeScaffold.getPosition(p).getAtom("CA"));
	      scaffoldCaCb.push_back(nonEpitopeScaffold.getPosition(p).getAtom("CB"));

	      double rmsd = inverseCaCb.rmsd(scaffoldCaCb);
	      if (rmsd < opt.rmsd){
		MSLOUT.fprintf(stdout,"HIT %1s,%3d%1s to %-10s %1s,%3d%1s %8.3f\n",
			       inverseCaCb[0].getChainId().c_str(),
			       inverseCaCb[0].getResidueNumber(),
			       inverseCaCb[0].getResidueICode().c_str(),
			       MslTools::getFileName(opt.scaffold).c_str(),
			       scaffoldCaCb[0].getChainId().c_str(),
			       scaffoldCaCb[0].getResidueNumber(),
			       scaffoldCaCb[0].getResidueICode().c_str(),
			       rmsd
			       );

		// Check for clashes with glycan?
		

	      } // IF RMSD
	      
	    } // END nonEptiopeScaffold Positions

	    
	  } // END Alternative Conformation/Rotamer

	} // END InverseRotamer Positions	
	MSLOUT.stream() << "Done."<<endl;
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
		cout << "findPositionsWithRotamers --scaffold PDB --ref PDB --rotlib ROTLIB\n";
		exit(0);
	}

	opt.scaffold = OP.getString("scaffold");
	if (OP.fail()){
		cerr << "ERROR 1111 scaffold not specified.\n";
		exit(1111);
	}
	opt.ref = OP.getString("ref");
	if (OP.fail()){
		cerr << "ERROR 1111 ref not specified.\n";
		exit(1111);
	}
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
	  opt.rotlib =	SYSENV.getEnv("MSL_ROTLIB");
	  MSLOUT.stream() << "MSL_ROTLIB is set to: "<<opt.rotlib<<endl;
	}

	opt.numRotamers = OP.getInt("numRotamers");
	if (OP.fail()){
	  opt.numRotamers = 10;
	}
	
	opt.refSelect = OP.getString("refSelect");
	if (OP.fail()){
	  opt.refSelect = "chain Z";
	  MSLOUT.stream() << "WARNING refSelect defaulted to "<<opt.refSelect<<endl;
	}

	opt.epitopeSelect = OP.getString("epitopeSelect");
	if (OP.fail()){
	  opt.epitopeSelect = "chain Z and name CA";
	  MSLOUT.stream() << "WARNING epitopeSelect defaulted to "<<opt.epitopeSelect<<endl;
	}
	
	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  opt.rmsd = 0.5;
	  MSLOUT.stream() << "WARNING rmsd defaulted to "<<opt.rmsd<<endl;
	  
	}

	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}

	
	return opt;
}

void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers){
  
	// Create a pdb topology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(_opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);


	AtomSelection sel(_sys.getAtomPointers());
	
	System local;
	local.addAtoms(sel.select(_opt.rotamerSelect));

	// Each position build opt.numRotamer rotamers
	for (uint i = 0; i < local.positionSize();i++){

	  cout << "I: "<<i<<endl;
	  Residue &res = local.getResidue(i);
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

	  }
	}
}

