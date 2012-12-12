
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "PhiPsiStatistics.h"
#include "getVectorPairs.h"
#include "VectorPair.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
#include "SasaCalculator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
using namespace MSL;
using namespace MslTools;

int main(int argc, char *argv[]) {

    Options opt = setupOptions(argc, argv);

    vector<string> lines;
    MslTools::readTextFile(lines,opt.pdblist);

    fprintf(stdout, "%20s %6s %5s %8s %8s %8s %6s %5s %8s %8s %8s %12s %6s %8s %8s %8s %8s %8s %8s %8s %8s\n","FILE","CHAIN1", "RESI1", "PHI1", "PSI1", "SASA1", "CHAIN2", "RESI2", "PHI2", "PSI2", "SASA2", "RESN1-RESN2", "SEQSEP", "D1", "D2", "A1", "A2", "A3", "A4", "T1", "T2");
    for (uint i = 0; i < lines.size();i++){

	System sys;
	sys.readPdb(lines[i]);

	SasaCalculator sasaCalc(sys.getAtomPointers(), 1.4, 200);
	sasaCalc.calcSasa();	  
	

	// Each chain
	for (uint c = 0; c < sys.chainSize();c++){

	    // Walk through residues
	    Chain &ch = sys.getChain(c);


	    for (uint r = 1; r < ch.positionSize()-1;r++){
	        
	        Residue &resm1   = ch.getResidue(r-1);
		Residue &res     = ch.getResidue(r);
		Residue &resp1   = ch.getResidue(r+1);

		if (!(resm1.atomExists("CA") && res.atomExists("CA") && resp1.atomExists("CA"))) {continue;}

		if (!(resm1.atomExists("N") && resm1.atomExists("C") &&
		      res.atomExists("N") && res.atomExists("C") &&
		      resp1.atomExists("N") && resp1.atomExists("C"))){ continue;}
		
		double phi1        = PhiPsiStatistics::getPhi(resm1,res);
		double psi1        = PhiPsiStatistics::getPsi(res,resp1);

	      for (uint c2 = c; c2 < sys.chainSize();c2++){
		Chain &ch2 = sys.getChain(c2);

		for (uint r2 = 1; r2 < ch2.positionSize()-1;r2++){
		  if (c == c2 && r >= r2) continue;

		  Residue &res2m1   = ch2.getResidue(r2-1);
		  Residue &res2     = ch2.getResidue(r2);
		  Residue &res2p1   = ch2.getResidue(r2+1);

		if (!(res2m1.atomExists("CA") && res2.atomExists("CA") && res2p1.atomExists("CA"))) {continue;}

		if (!(res2m1.atomExists("N") && res2m1.atomExists("C") &&
		      res2.atomExists("N") && res2.atomExists("C") &&
		      res2p1.atomExists("N") && res2p1.atomExists("C"))){ continue;}

		
		  double phi2        = PhiPsiStatistics::getPhi(res2m1,res2);
		  double psi2        = PhiPsiStatistics::getPsi(res2,res2p1);


		  // Skip if no 'CA' 'CB' atoms
		  //if (!(res.atomExists("CA") && res.atomExists("CB") && res2.atomExists("CA") && res2.atomExists("CB"))) continue;

		  Atom *CB1 = NULL;
		  if (res.getResidueName() == "GLY"){
		    CB1 = PDBTopology::getPseudoCbeta(res);
		  }
		  if (res.atomExists("CB")){
		    CB1 = &res.getAtom("CB");
		  }
		  if (CB1 == NULL) continue;

		  Atom *CB2 = NULL;
		  if (res2.getResidueName() == "GLY"){
		    CB2 = PDBTopology::getPseudoCbeta(res2);
		  }
		  if (res2.atomExists("CB")){
		    CB2 = &res2.getAtom("CB");
		  }
		  if (CB2 == NULL) continue;

		  VectorPair vp(res.getAtom("CA").getCoor(),
				CB1->getCoor(),
				res2.getAtom("CA").getCoor(),
				CB2->getCoor());


		  vp.calcAll();
		  if (res.getResidueName() == "GLY"){
		    delete(CB1);
		  }
		  if (res2.getResidueName() == "GLY"){
		    delete(CB2);
		  }

		  if (opt.maxCbCb != MslTools::doubleMax && vp.getDistance2() > opt.maxCbCb) continue; 

		  int seqsep = 0;
		  if (c == c2){
		    seqsep = r2-r;
		  }

		  string chainId1 = res.getChainId();
		  if (chainId1 == ""){
		    chainId1 = "_";
		  }
		  string chainId2 = res2.getChainId();
		  if (chainId2 == ""){
		    chainId2 = "_";
		  }

		  // Skip chain-chain pairs?
		  if (opt.interChainOnly && chainId1 == chainId2) continue;

		fprintf(stdout, "%20s %6s %4d%1s %8.2f %8.2f %8.2f %6s %4d%1s %8.2f %8.2f %8.2f  %5s-%-5s %6d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
			MslTools::getFileName(lines[i]).c_str(),
			chainId1.c_str(),
			res.getResidueNumber(),
			res.getResidueIcode().c_str(),
			phi1,psi1,
			sasaCalc.getResidueSasa(res.getPositionId()),
			chainId2.c_str(),
			res2.getResidueNumber(),
			res2.getResidueIcode().c_str(),
			phi2,psi2,
			sasaCalc.getResidueSasa(res2.getPositionId()),
			res.getResidueName().c_str(),
			res2.getResidueName().c_str(),
			seqsep,
			vp.getDistance1(),
			vp.getDistance2(),
			vp.getAngle1(),
			vp.getAngle2(),
			vp.getAngle3(),
			vp.getAngle4(),
			vp.getTorsion1(),
			vp.getTorsion2());


		} // END r2
	      } // END c2
	    } // END r
	} // END c
    } // END file
}

Options setupOptions(int theArgc, char * theArgv[]){
    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: getVectorPairs " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdblist PDB\n";
	cout << "maxCbCb DIST\n";
	cout << "interChainOnly\n";
	cout << endl;
	exit(0);
    }

    opt.pdblist = OP.getString("pdblist");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdblist specified."<<endl;
	exit(1111);
    }
    opt.maxCbCb = OP.getDouble("maxCbCb");
    if (OP.fail()){
      opt.maxCbCb = MslTools::doubleMax;
    }
    opt.interChainOnly = OP.getBool("interChainOnly");
    if (OP.fail()){
      opt.interChainOnly = false;
    }
    return opt;
}
