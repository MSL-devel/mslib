#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "RotamerLibrary.h"
#include "RotamerLibraryBuilder.h"
#include "RotamerLibraryWriter.h"
#include "generateRotamerLibrary.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("generateRotamerLibrary");

int main(int argc, char *argv[]) {

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnAllOff();

	// Read cmdline options
	Options opt = setupOptions(argc,argv);

	// Read in rotamer library for building rotamers on each position
	RotamerLibrary rotLib;
	rotLib.readFile(opt.rotlib);

	// Rotamer library to fill up with position specific rotamer libraries
	RotamerLibrary newRotLib;
	newRotLib.readFile(opt.rotlib);
	newRotLib.removeAllConformations();

	// A utility object for building rotamer libraries from residue conformations
	RotamerLibraryBuilder rBuild(&newRotLib);

	// Create a pdb topology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);
	

	// Read in the pdb with a set of functional groups defined
	System inSys;
	inSys.readPdb(opt.pdb);

	// Stip to poly-Alanine
	AtomPointerVector ala;
	for (uint i = 0; i < inSys.positionSize();i++){
	  Position &pos = inSys.getPosition(i);
	  if (!pos.atomExists("N") || !pos.atomExists("CA") || !pos.atomExists("C") || !pos.atomExists("O")) {
	    cout << "Skipping position: "<<pos.toString()<<endl;
	    continue;
	  }
	  ala.push_back(new Atom(pos.getAtom("N"))); ala.back()->setResidueName("XXX");
	  ala.push_back(new Atom(pos.getAtom("CA"))); ala.back()->setResidueName("XXX");
	  ala.push_back(new Atom(pos.getAtom("C"))); ala.back()->setResidueName("XXX");
	  ala.push_back(new Atom(pos.getAtom("O"))); ala.back()->setResidueName("XXX");

	  if (pos.atomExists("CB")){
	    ala.push_back(new Atom(pos.getAtom("CB"))); ala.back()->setResidueName("XXX");
	  } else {
	    ala.push_back(new Atom(*pdbTop.getPseudoCbeta(pos.getCurrentIdentity())));ala.back()->setResidueName("XXX");
	  }
	  
	}

	// Create the poly-alanine system
	System sys;
	sys.addAtoms(ala);

	char str[255];
	strcpy( str, opt.report.c_str());
	ofstream fout(str);
	fout << "# =================  POSITION ANALYSIS  ================="<<endl;
	stringstream goodRotamerReport;

	// For each residue
	for (uint i = 0; i < sys.positionSize();i++){

	  Position &pos = sys.getPosition(i);

	  // Get backbone atoms  
	  AtomPointerVector backboneAtoms = pdbTop.getBackboneAtoms(pos.getResidue("XXX"));

	  // For each amino acid in rotamer library
	  stringstream aaOutput;

	  set<string> AAs = rotLib.getAllResList();
	  set<string>::iterator it;
	  for (it = AAs.begin(); it != AAs.end();it++){
	    
	    MSLOUT.stream() << "Going to build: "<<*it<<" conformations:"<<rotLib.size(rotLib.getLibraryNames()[0],*it)<<" at position "<<pos.getPositionId()<<endl;

	    // Get all the rotamers from the rotamer library...
	    AtomContainer ac = pdbTop.getResidue(pos.getPositionId()+","+*it,backboneAtoms,rotLib.size(rotLib.getLibraryNames()[0],*it));

	    pos.addIdentity(ac.getAtomPointers(),*it);
	    pos.setActiveIdentity(*it);
	    Residue &res = pos.getCurrentIdentity();
	    goodRotamerReport << "GOOD_ROTAMERS: "<<res.getIdentityId()<<" ";

	    int numberRotamersAccepted = 0;
	    for (uint c = 0; c < res.getNumberOfAltConformations();c++){
	      res.setActiveConformation(c);

	      // Clash check
	      bool clash = false;
	      for (uint j = 0; j < sys.positionSize();j++){
		//if (j == pos.getIndexInSystem() || j == pos.getIndexInSystem()-1 || j == pos.getIndexInSystem()+1) continue; // skip this position
		if (j == pos.getIndexInSystem()) continue; // skip this position

		Residue &alaCheck = sys.getPosition(j).getResidue("XXX");

		// For each atom in res
		for (uint t = 0; t < res.size();t++){
		  if (res.getAtom(t).getName().substr(0,1) == "H") continue; // skip hydrogens
		  if (res.getAtom(t).getName() == "N" ||res.getAtom(t).getName() == "C" ||res.getAtom(t).getName() == "O" ||res.getAtom(t).getName() == "CA") continue; //skip backbone atoms

		  if (res.getAtom(t).distance(alaCheck.getAtom("CA")) < 2.5 || res.getAtom(t).distance(alaCheck.getAtom("CB")) < 2.5 ||
		      res.getAtom(t).distance(alaCheck.getAtom("C")) < 2.5  || res.getAtom(t).distance(alaCheck.getAtom("O")) < 2.5  ||
		      res.getAtom(t).distance(alaCheck.getAtom("N")) < 2.5){

		    if (opt.debug) {
		      MSLOUT.stream() << "RES CLASH: "<<res.getAtom(t).toString()<< " and " <<alaCheck.toString()<<" distC: "<<res.getAtom(t).distance(alaCheck.getAtom("C"))<< " distN: "<<res.getAtom(t).distance(alaCheck.getAtom("N"))<<endl;
		    // PDBWriter pout;
		    // pout.open("/tmp/test.pdb");
		    // pout.write(res.getAtomPointers());
		    // pout.close();

		    }
		    clash = true;

		    break;
		  }
		} // end for each atom in res
		if (clash) break;

	      } // end for position

	      if (clash){
		if (opt.debug && res.getResidueNumber() == 4){
		  MSLOUT.stream() << res.toString()<< " rotamer "<<c<<" has a clash"<<endl;
		}
	      } else {

		char tmp[80];
		sprintf(tmp,"%1s_%d",res.getChainId().c_str(),res.getResidueNumber());
		rBuild.addRotamer(res, rotLib.getLibraryNames()[0],(string)tmp);
		numberRotamersAccepted++;

		goodRotamerReport << c <<" ";

		if (opt.debug){
		  sprintf(tmp, "/tmp/rotamer-%1s_%04d_%3s_%03d.pdb", res.getChainId().c_str(),res.getResidueNumber(),res.getResidueName().c_str(),c);
		  PDBWriter pout;
		  pout.open((string)tmp);
		  pout.write(res.getAtomPointers());
		  pout.close();
		}
	      }
	    } // each rotamer

	    char tmp[80];
	    sprintf(tmp, "%3s_%03d ", (*it).c_str(),numberRotamersAccepted);
	    aaOutput << tmp;

	    goodRotamerReport << endl;
	  } // each amino acid	  

	  char tmp[1000];
	  sprintf(tmp, "%1s %4d%1s %s\n", pos.getChainId().c_str(), pos.getResidueNumber(),pos.getResidueIcode().c_str(), aaOutput.str().c_str());
	  fout << tmp;

	} // each position in the system

	fout <<endl<<endl;
	fout << "# ================= GOOD ROTAMER REPORT from ("<<opt.rotlib<<") ================="<<endl;
	fout << goodRotamerReport.str()<<endl;
	fout.close();

	newRotLib.writeFile(opt.outputLibrary);
	cout << "Done."<<endl;
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
		cout << "generateRotamerLibrary --pdb PDB --rotlib ROTLIB [ --outputLibrary PDB.rotlib ]\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "ERROR 1111 rotlib not specified.\n";
		exit(1111);
	}

	opt.outputLibrary = OP.getString("outputLibrary");
	if (OP.fail()){
	  opt.outputLibrary = MslTools::getFileName(opt.pdb)+".rotlib";
	  cout << "Output rotamer library will be saved as "<<opt.outputLibrary<<endl;
	}

	opt.report = OP.getString("report");
	if (OP.fail()){
	  opt.report = MslTools::getFileName(opt.pdb)+".report.txt";
	  cout << "Output report will be saved as "<<opt.report<<endl;
	}
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}
	return opt;
}
