#include <iostream>
#include <cstdlib>

#include "ChiStatistics.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "PhiPsiReader.h"
#include "PhiPsiStatistics.h"
#include "SurfaceAreaAndVolume.h"
#include "SasaCalculator.h"
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
#include "glycineSearch.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("glycineSearch");

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

    PhiPsiStatistics pps;

    if (opt.phiPsiTable != ""){
	    PhiPsiReader ppr(opt.phiPsiTable);
	    ppr.open();
	    ppr.read();
	    ppr.close();
	    pps = ppr.getPhiPsiStatistics();
    }

	    for (uint i = 0 ; i < sys.positionSize();i++){

	    Residue & n   = sys.getResidue(i);

	    // Remove non-amino acids... bad way to do this, but should work.
	    string oneLetter = MslTools::getOneLetterCode(n.getResidueName());
	    if (oneLetter == "X") continue;

	    //fprintf(stdout, "%s %1s %3d%1s %3s ",filename.c_str(),n.getChainId().c_str(),n.getResidueNumber(),n.getResidueIcode().c_str(),n.getResidueName().c_str());


	    double phi           = 0.0;
	    double psi    	 = 0.0;
	    if (i > 0 && (i < sys.positionSize()-1  &&  
			  MslTools::getOneLetterCode(sys.getResidue(i+1).getResidueName()) != "X"  && 
			  MslTools::getOneLetterCode(sys.getResidue(i-1).getResidueName()) != "X")){
		    Residue & nm1 = sys.getResidue(i-1);
		    Residue & np1 = sys.getResidue(i+1);
		    if (nm1.getChainId() == n.getChainId() && np1.getChainId() == n.getChainId()){
			    phi     = PhiPsiStatistics::getPhi(nm1,n);
			    psi     = PhiPsiStatistics::getPsi(n,np1);

			    if (opt.phiPsiTable != ""){
			      int resCounts  = pps.getCounts(n.getResidueName(),phi,psi);
			      int glyCounts  = pps.getCounts("GLY",phi,psi);
			      double glyFreq = pps.getFreqInBin("GLY",phi,psi);
			      string extra="";
			      if (glyFreq > 0.25){
				extra=" **** ";
			      }
			      fprintf(stdout, "%10s %8d %8d %8.3f %s\n",n.getPositionId().c_str(),resCounts,glyCounts,glyFreq,extra.c_str());
				    
			    }
		    }
		    
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
		cout << "buildInverseRotamers --pdb PDB --rotlib ROTLIB\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.phiPsiTable = OP.getString("phiPsiTable");
	if (OP.fail()){
		cerr << "ERROR 1111 phiPsiTable not specified.\n";
		exit(1111);
	}
	return opt;
}

