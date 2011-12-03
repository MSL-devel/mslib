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

        // Parse command-line options
	Options opt = setupOptions(argc,argv);

	ofstream fout;
	fout.open(opt.outFile.c_str());

	// Read in the pdb with epitope
	System scaffold;
	scaffold.readPdb(opt.scaffold);

	// Use a part of the scaffold structure for alignment, rest for finding positions
	AtomPointerVector epitopeAts;
	AtomSelection scaffoldSel(scaffold.getAtomPointers());
	if (opt.epitopeSelect != ""){
	  epitopeAts = scaffoldSel.select(opt.epitopeSelect);
	}

	// Reference structure, holds epitope.
	System ref;
	ref.readPdb(opt.ref);

	// Use a part of the reference structure for alignment
	// -- use residue numbers from epitopeAts , they will correspond.
	AtomSelection refSelect(ref.getAtomPointers());
	string refEpitopeSelectStr = MslTools::stringf("chain A and resi %03d-%03d and name CA",
								    epitopeAts(0).getResidueNumber(),
								    epitopeAts(epitopeAts.size()-1).getResidueNumber());
	AtomPointerVector refAts  = refSelect.select(refEpitopeSelectStr);
	AtomPointerVector glycans = refSelect.select(opt.glycanSelect);

	// Make sure selections were the same size..
	if (epitopeAts.size() != refAts.size()){
	  MSLOUT.stream() << "ERROR epitope selection and reference selection do not contain same number of atoms!"<<endl;
	  MSLOUT.stream() << "ERROR epitope selection ("<<opt.epitopeSelect<<") got "<<epitopeAts.size()<<" atoms."<<endl;
	  MSLOUT.stream() << "ERROR ref               ("<<refEpitopeSelectStr<<") got "<<refAts.size()    <<" atoms."<<endl;
	  exit(2222);
	}

	// Align 
	Transforms tr;
	if (!tr.rmsdAlignment(epitopeAts,refAts,scaffold.getAtomPointers())){
	  MSLOUT.stream() << "ERROR aligning epitope("<<opt.epitopeSelect<<") and ref("<<refEpitopeSelectStr<<")"<<endl;
	  exit(3333);
	}

	// For debugging store RMSD of initial alignment
	scaffold.writePdb("scaffold_aligned.pdb");

	// Create a pdb topology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);

	AtomSelection sel(ref.getAtomPointers());
	
	System local;
	local.addAtoms(sel.select(opt.rotamerSelect));

	// Each position build opt.numRotamer rotamers
	for (uint i = 0; i < local.positionSize();i++){

	  
	  Residue &res = local.getResidue(i);
	  if (res.getResidueName() == "HIS"){ res.setResidueName("HSD");}
	  string id = res.getIdentityId();
	 

	  // Build using inverse-rotamer IC building 
	  //AtomContainer ac = pdbTop.getResidue(id,res.getAtomPointers(),opt.numRotamers);
	  
	  // Build a generic N-CA-C, then build rotmaers onto it
	  MSLOUT.stream() << "BUILD GENERICS"<<endl;
	  AtomContainer ac = pdbTop.getGenericResidue(id,opt.numRotamers);

	  MSLOUT.stream() << "Done getting generic residues: "<<ac.getAtomPointers().getMaxAltConf()<<endl;
	    


	 
	  AtomContainer *lastRes = new AtomContainer();
	  lastRes->addAtoms(ac.getAtomPointers());


	  for (uint j = 0; j < (*lastRes)(0).getNumberOfAltConformations();j++){

	    // Change each atoms active conformation
	    for (uint k =0;k <(*lastRes).size();k++){
	      (*lastRes)(k).setActiveConformation(j);
	    }

	    // Need smartRmsdAlignment to take a filter for atom types to align.
	    Transforms t;
	    if (!t.smartRmsdAlignment((*lastRes).getAtomPointers(),res.getAtomPointers(),Transforms::MT_ATOMNAME)){
	      MSLOUT.stream() << "ERROR WITH ALIGNING lastRes and res"<<endl;
	    }


	    // Write out inverse rotamers to check that they were done properly
	    if (opt.debug){

	      PDBWriter pout;
	      pout.open(MslTools::stringf("ir_%02d_%02d.pdb",i,j));
	      pout.write((*lastRes).getAtomPointers());
	      pout.close();
	    }

	    // RMS check..
	    AtomPointerVector inverseCaCb;
	    for (uint a = 0; a < (*lastRes).size();a++){
	      if ((*lastRes)(a).getName() == "CA" ||
		  (*lastRes)(a).getName() == "CB") {
		inverseCaCb.push_back(&(*lastRes)[a]);
	      }
	    }


	    for (uint p = 0; p < scaffold.positionSize();p++){
	      Position &pos = scaffold.getPosition(p);
	      
	      if (!(pos.atomExists("CA"))) continue;
	       
	      if (pos.getCurrentIdentity().getResidueName() != "GLY" && !pos.atomExists("CB")) continue;

	    



	      //MSLOUT.stream() << "Testing "<<pos.toString()<<endl;
	      AtomContainer scaffoldCaCb;
	      scaffoldCaCb.addAtom(pos.getAtom("CA"));

	      if (scaffold.getPosition(p).getCurrentIdentity().getResidueName() == "GLY"){
		if (pos.atomExists("N") && pos.atomExists("C")){
		  CartesianPoint CBcoor = CartesianGeometry::build(pos.getAtom("CA").getCoor(), pos.getAtom("N").getCoor(), pos.getAtom("C").getCoor(), 1.521, 110.5, -122.5);
		  Atom cb(pos.getCurrentIdentity().getIdentityId()+",CB",CBcoor);
		  scaffoldCaCb.addAtom(cb);
		} else {
		  MSLOUT.stream() << "Postion "<<scaffold.getPosition(p)<<" does not have a C,N atoms and is a GLY residue."<<endl;
		  continue;
		}
		
	      } else {


		if (pos.atomExists("CB")){
		  scaffoldCaCb.addAtom(pos.getAtom("CB"));
		} else {
		  MSLOUT.stream() << "Postion "<<scaffold.getPosition(p)<<" does not have a CB and is not a GLY residue."<<endl;
		  continue;
		}
	      }


	      
	      double rmsd = scaffoldCaCb.getAtomPointers().rmsd(inverseCaCb);
	      if (rmsd < opt.rmsd){


		// Scan Hit for Glycan clash
		int numClashes = 0;
		for (uint g = 0; g < glycans.size();g++){
		  for (uint s = 0; s < scaffold.atomSize();s++){
		    double g_dist = glycans(g).distance2(scaffold.getAtom(s));
		    if (g_dist < 4.0){
		      numClashes++;
		    }
		  } // END SCAFFOLD ATOMS
		} // END GLYCAN ATOMS


		if (numClashes < opt.numGlycanClashesAllowed){
		  
		  fout << MslTools::stringf("HIT %1s,%3d%1s to %-10s %1s,%3d%1s %8.3f\n",
			       (*lastRes)(0).getChainId().c_str(),
			       (*lastRes)(0).getResidueNumber(),
			       (*lastRes)(0).getResidueIcode().c_str(),
			       MslTools::getFileName(opt.scaffold).c_str(),
			       scaffoldCaCb[0].getChainId().c_str(),
			       scaffoldCaCb[0].getResidueNumber(),
			       scaffoldCaCb[0].getResidueIcode().c_str(),
			       rmsd
			       );

		  
		}

	      } // IF RMSD else 
	    } // END nonEptiopeScaffold Positions

	  } // END Alternative Conformation/Rotamer

	} // END InverseRotamer Positions	
	fout.close();
	MSLOUT.stream() << "Done."<<endl;
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile
	OP.readArgv(theArgc, theArgv);


	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
	  OP.readFile(opt.configfile);
	  if (OP.fail()) {
	    cerr << "ERROR couldn't read : "<<opt.configfile<<endl;
	    exit(1);
	  }
	}




	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "findPositionsWithRotamers --scaffold PDB --ref PDB --rotlib ROTLIB\n";
		exit(0);
	}

	//OP.printConfFile();
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
	
	opt.epitopeSelect = OP.getString("epitopeSelect");
	if (OP.fail()){
	  opt.epitopeSelect = "chain Z and name CA";
	  MSLOUT.stream() << "WARNING epitopeSelect defaulted to "<<opt.epitopeSelect<<endl;
	}

	opt.rotamerSelect = OP.getString("rotamerSelect");
	if (OP.fail()){
	  opt.rotamerSelect = "chain Z";
	  MSLOUT.stream() << "WARNING rotamerSelect defaulted to "<<opt.rotamerSelect<<endl;
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

	opt.numGlycanClashesAllowed = OP.getInt("numGlycanClashesAllowed");
	if (OP.fail()){
	  opt.numGlycanClashesAllowed = 2;
	  MSLOUT.stream() << "WARNING numGlycanClashesAllowed defaulted to: "<<opt.numGlycanClashesAllowed<<endl;
	}

	opt.glycanSelect = OP.getString("glycanSelect");
	if (OP.fail()){
	  opt.glycanSelect = "not resn ALA+CYS+ASP+GLU+PHE+GLY+HIS+ILE+LYS+LEU+MET+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR";
	  MSLOUT.stream() << "WARNING glycanSelect defaulted to: '"<<opt.glycanSelect<<"'"<<endl;
	}

	opt.outFile = OP.getString("outFile");
	if (OP.fail()){
	  opt.outFile = "data.txt";
	  MSLOUT.stream() << "WARNING outFile defaulted to: "<<opt.outFile<<endl;
	}
	return opt;
}



