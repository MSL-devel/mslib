#include <iostream>
#include <cstdlib>

#include "MslTools.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Timer.h"
#include "SasaCalculator.h"
#include "AtomPointerVector.h"
#include "patchAnalysis.h"


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("patchAnalysis");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnOn("patchAnalysis");

	Timer t;
	double start = t.getWallTime();
	
	vector<string> lines;
	if (opt.list != ""){
	  MslTools::readTextFile(lines,opt.list);
	} else {
	  lines.push_back(opt.pdb);
	}

	for (uint i = 0; i < lines.size();i++){
	  System pdb;
	  pdb.readPdb(lines[i]);
	  AtomSelection sel(pdb.getAtomPointers());

	  AtomPointerVector patch1  = sel.select(opt.patch1);
	  AtomPointerVector patch2  = sel.select(opt.patch2);
	  AtomPointerVector patches = patch1 + patch2;

	  AtomContainer patchesC;
	  patchesC.addAtoms(patches);
	  SasaCalculator sasaComplex(patchesC.getAtomPointers(),1.4,2000);
	  sasaComplex.calcSasa();
	  double boundSasa = sasaComplex.getTotalSasa();

	  AtomContainer patch1_new;
	  patch1_new.addAtoms(patch1);
	  SasaCalculator sasaPatch1(patch1_new.getAtomPointers(),1.4,2000);
	  sasaPatch1.calcSasa();
	  double patch1Sasa = sasaPatch1.getTotalSasa();

	  AtomContainer patch2_new;
	  patch2_new.addAtoms(patch2);
	  SasaCalculator sasaPatch2(patch2_new.getAtomPointers(),1.4,2000);
	  sasaPatch2.calcSasa();
	  double patch2Sasa = sasaPatch2.getTotalSasa();

	  //SasaCalculator sasaTotal(pdb.getAtomPointers(),1.4,2000);
	  //sasaTotal.calcSasa();
	  //double totalSasa = sasaTotal.getTotalSasa();

	  //patch1_new.writePdb("patch1.pdb");
	  //patch2_new.writePdb("patch2.pdb");
	  //patchesC.writePdb("patch_complex.pdb");

	  fprintf(stdout, "%35s %8.2f %10d %10d %10d %10d\n",lines[i].c_str(),(patch1Sasa+patch2Sasa)  - boundSasa, patch1.size(),patch2.size(),patches.size(),pdb.atomSize());
	  //fprintf(stdout, "Patches1 table: \n%s\n", sasaPatch1.getSasaTable().c_str());
	  //fprintf(stdout, "PatchesC table: \n%s\n", sasaComplex.getSasaTable().c_str());
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
		cout << "patchAnalysis --pdb foo.pdb --patch1 \"chain A and resi 5-150\" --patch2 \"chain B and resi 5-150\" \n";
		exit(0);
	}

	opt.list= OP.getString("list");
	if (OP.fail()){
	  opt.list = "";
	}
	opt.pdb= OP.getString("pdb");
	if (OP.fail() && opt.list == ""){
		cerr << "ERROR 1111 pdb or list not specified.\n";
		exit(1111);
	}
	opt.patch1 = OP.getString("patch1");
	if (OP.fail()){
		cerr << "ERROR 1111 patch1 not specified.\n";
		exit(1111);
	}
	opt.patch2 = OP.getString("patch2");
	if (OP.fail()){
		cerr << "ERROR 1111 patch2 not specified.\n";
		exit(1111);
	}
	return opt;
}
