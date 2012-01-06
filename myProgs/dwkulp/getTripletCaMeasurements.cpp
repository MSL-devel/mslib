
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "PhiPsiStatistics.h"
#include "getTripletCaMeasurements.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
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

    for (uint i = 0; i < lines.size();i++){

	System sys;
	sys.readPdb(lines[i]);


	// Each chain
	int loopCount = 1;
	for (uint c = 0; c < sys.chainSize();c++){

	    int brokenChain = 0;

	    stringstream ss;

	    // Walk through residues
	    Chain &ch = sys.getChain(c);
	    int loopIndex = -1;

	    int startIndex = 0;
	    if (opt.startResidue != ""){
	        if (ch.positionExists(opt.startResidue)){
			Position &pos = ch.getPosition(opt.startResidue);
			startIndex = pos.getIndexInChain();
		}
	
	    }
            int endIndex = ch.positionSize();
	    if (opt.endResidue != ""){
	        if (ch.positionExists(opt.endResidue)){
			Position &pos = ch.getPosition(opt.endResidue);
			endIndex = pos.getIndexInChain();
		}
	    }

	    for (uint r = startIndex; r < endIndex;r++){

		// Skip first and last residues
		if (r == 0 || r == ch.positionSize()-1) {ss.str(""); continue;}

		Residue &resm1 = ch.getResidue(r-1);
		Residue &res   = ch.getResidue(r);
		Residue &resp1 = ch.getResidue(r+1);

		// Skip if no 'CA' atom
		if (!(resm1.atomExists("CA") && res.atomExists("CA") && resp1.atomExists("CA"))) {ss.str("");continue;}

		if (!(resm1.atomExists("N") && resm1.atomExists("C") &&
		      res.atomExists("N") && res.atomExists("C") &&
		      resp1.atomExists("N") && resp1.atomExists("C"))){ ss.str("");continue;}
		      

		// Only look for loop residues
		//if (!(res("CA").getSegID() == "LLLL"  || res("CA").getSegID() == "TTTT"))  {ss.str("");loopIndex=-1;continue;}


		// Set loop index to -1 for first and last residue of loop
		//if (!(resm1("CA").getSegID() == "LLLL"  || resm1("CA").getSegID() == "TTTT"))  {
		//loopIndex = -1;
		//}

		int localLoopIndex = loopIndex;
		if (!(resp1("CA").getSegID() == "LLLL"  || resp1("CA").getSegID() == "TTTT"))  {
		    localLoopIndex = -1;
		}


		string inLoop ="NO_LOOP";
		if (res("CA").getSegID() == "LLLL") { inLoop="YESLOOP";}




	      
		double cacacaAngle = resm1("CA").angle(res("CA"),resp1("CA"));
		double caca1       = resm1("CA").distance(res("CA"));
		double caca2       = res("CA").distance(resp1("CA"));

		double psi1        = PhiPsiStatistics::getPsi(resm1,res);
		double phi2        = PhiPsiStatistics::getPhi(resm1,res);
		double psi2        = PhiPsiStatistics::getPsi(res,resp1);
		double phi3        = PhiPsiStatistics::getPhi(res,resp1);


		char tmp[300];
		sprintf(tmp, "%10d %s %1s %5d %3s %3s %3s L%-3d, %8.2f %8.2f %8.2f , %8.2f %8.2f %8.2f %8.2f %s\n", 
			loopCount,
			MslTools::getFileName(lines[i]).c_str(),
			res.getChainId().c_str(),
			res.getResidueNumber(),
			resm1.getResidueName().c_str(),
			res.getResidueName().c_str(),
			resp1.getResidueName().c_str(),
			localLoopIndex,
			cacacaAngle,
			caca1,
			caca2,
			psi1,
			phi2,
			psi2,
			phi3,inLoop.c_str());

		ss << tmp;

		if (abs(psi1) > 180 || abs(phi2) > 180){
		  cout <<"BAD ANGLES: "<<ss.str()<<endl;
		  exit(2132);
		}

		if (!(resp1("CA").getSegID() == "LLLL"  || resp1("CA").getSegID() == "TTTT"))  {

		    if (loopIndex+3 > opt.loopMin && loopIndex+3 < opt.loopMax){
			cout << ss.str();
			loopCount++;
		    } 

		    ss.str("");


		}

		if (caca1 > 4.00) {
		    brokenChain++;
		    cout << "BROKEN! "<<brokenChain<<endl;
		}
		loopIndex++;
	    } // END for startIndex to endIndex

	    if (brokenChain <= opt.numCaBreaksAllowed){
		std::cout << "COMPLETE_CHAIN: "<<lines[i]<<" chain "<<ch.getChainId()<<endl;
	    } else {
	      std::cout << "CHAIN BROKEN: "<<brokenChain<<endl;
	    }
	} // END chains



    }

	

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
	cout << "Usage: getTripletCaMeasurements " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdblist PDB\n";
	cout << endl;
	exit(0);
    }

    opt.pdblist = OP.getString("pdblist");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdblist specified."<<endl;
	exit(1111);
    }
    opt.loopMin = OP.getInt("minLoopLength");
    if (OP.fail()){
	opt.loopMin = 0;
    }
    opt.loopMax = OP.getInt("maxLoopLength");
    if (OP.fail()){
	opt.loopMax = MslTools::intMax;
    }

    opt.startResidue = OP.getString("startResidue");
    if (OP.fail()){
	opt.startResidue = "";
    }

    opt.endResidue = OP.getString("endResidue");
    if (OP.fail()){
	opt.endResidue = "";
    }
    opt.numCaBreaksAllowed = OP.getInt("numCaBreaksAllowed");
    if (OP.fail()){
	opt.numCaBreaksAllowed = 0;
    }
    opt.naturalBreaks = OP.getMultiString("naturalBreaks");
    
    return opt;
}
