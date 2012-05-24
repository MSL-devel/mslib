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
#include "FastaReader.h"
#include "Position.h"
#include "domainSasa.h"


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("domainSasa");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnOn("domainSasa");

	Timer t;
	double start = t.getWallTime();
	
	System pdb;
	pdb.readPdb(opt.pdb);
	pdb.getChain(0).setChainId("A");
	double dist = pdb.getAtom("A,112,CA").distance(pdb.getAtom("A,972,CA"));
	fprintf(stdout, "%s %8.3f\n", MslTools::getFileName(opt.pdb).c_str(), dist);
	// SasaCalculator sasaRef(pdb.getAtomPointers(), 1.4, 200);
	// sasaRef.calcSasa();
	// double totalSasa = sasaRef.getTotalSasa();

	// AtomSelection sel(pdb.getAtomPointers());	
	// for (uint i = 0; i < opt.domains.size();i++){


	//   AtomPointerVector domainAts    = sel.select(opt.domains[i]);
	//   AtomPointerVector notDomainAts = sel.select("not "+opt.domains[i]);
	//   AtomContainer dAts;
	//   dAts.addAtoms(domainAts);
	//   dAts.writePdb("domain.pdb");
	//   AtomContainer nAts;
	//   nAts.addAtoms(notDomainAts);
	//   nAts.writePdb("other.pdb");
	  
	//   SasaCalculator scDomain(dAts.getAtomPointers(), 1.4, 200);
	//   scDomain.calcSasa();

	//   SasaCalculator scOther(nAts.getAtomPointers(), 1.4, 200);
	//   scOther.calcSasa();

	//   double unBoundSasa = scDomain.getTotalSasa() + scOther.getTotalSasa();


	//   fprintf(stdout, "%s %s %8.3f %8.3f %8.3f %8.3f %8.3f\n", MslTools::getFileName(opt.pdb).c_str(), opt.domains[i].c_str(),scDomain.getTotalSasa(), scOther.getTotalSasa(),unBoundSasa, totalSasa, unBoundSasa - totalSasa);
	  
	  
	// } // FOR OPT.DOMAINS


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
		cout << "domainSasa --pdb foo.pdb --domains \"chain A and resi 5-150\" \n";
		exit(0);
	}
	opt.pdb= OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.domains = OP.getMultiString("domains");
	if (OP.fail()){
		cerr << "ERROR 1111 domains not specified.\n";
		exit(1111);
	}
	return opt;
}
