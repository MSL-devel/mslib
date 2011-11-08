
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "writeShiftedRamaStats.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
#include "PhiPsiStatistics.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
using namespace MSL;
using namespace MslTools;

int main(int argc, char *argv[]) {

    Options opt = setupOptions(argc, argv);
    cout << "READING: "<<opt.pdblist<<endl;

    vector<string> lines;
    MslTools::readTextFile(lines,opt.pdblist);

    map<string,bool> residueCombinations;
    PhiPsiStatistics pps;
    for (uint i = 0; i < lines.size();i++){

	System sys;
	sys.readPdb(lines[i]);

	bool brokenChain = false;
	// Each chain
	int loopCount = 1;
	for (uint c = 0; c < sys.chainSize();c++){

	    // Walk through residues
	    Chain &ch = sys.getChain(c);

	    int startIndex = 0;
            int endIndex = ch.positionSize();

	    for (uint r = startIndex; r < endIndex;r++){

		// Skip first and last residues
		if (r == 0 || r == ch.positionSize()-1) {continue;}

		Residue &resm1 = ch.getResidue(r-1);
		Residue &res   = ch.getResidue(r);
		Residue &resp1 = ch.getResidue(r+1);

		// Skip if no backbone atoms
		if (!(resm1.atomExists("N") && res.atomExists("N") && resp1.atomExists("N"))) {continue;}
		if (!(resm1.atomExists("CA") && res.atomExists("CA") && resp1.atomExists("CA"))) {continue;}
		if (!(resm1.atomExists("C") && res.atomExists("C") && resp1.atomExists("C"))) {continue;}


		// Only look for loop residues
		if (!(res("CA").getSegID() == "LLLL"  || res("CA").getSegID() == "TTTT"))  {continue;}

		double psi1        = PhiPsiStatistics::getPsi(resm1,res);
		double phi2        = PhiPsiStatistics::getPhi(resm1,res);
		double psi2        = PhiPsiStatistics::getPsi(res,resp1);
		double phi3        = PhiPsiStatistics::getPhi(res,resp1);

		stringstream ss1;
		ss1 << resm1.getResidueName()<<"-"<<res.getResidueName();

		double psi1round = MslTools::smartRound(psi1,20);
		if (abs(psi1round) > 180){
		  cerr << "ERROR psi1 = "<<psi1<<" rounded: "<<psi1round<<endl;
		  cerr << "\t"<<res.toString()<<endl;
		  exit(0);
		}

		stringstream psi1ss;
		psi1ss << psi1round;

		double phi2round = MslTools::smartRound(phi2,20);
		if (abs(phi2round) > 180){
		  cerr << "ERROR phi2 = "<<phi2<<" rounded: "<<phi2round<<endl;
		  cerr << "\t"<<res.toString()<<endl;
		  exit(0);
		}

		stringstream phi2ss;
		phi2ss << phi2round;
		
		//cout << "ADDING "<<ss1.str()<<" "<<psi1ss.str()<<" "<<phi2ss.str()<<endl;
		pps.addStatisitics(ss1.str(),psi1ss.str(),phi2ss.str(),1);

		residueCombinations[ss1.str()] = true;

		
	    }


	}



    }//lines[i]

    map<string,int> counts = pps.getPhiPsiCounts();
    RandomNumberGenerator rng;
    rng.setTimeBasedSeed();
    
    ofstream fout;
    fout.open("shiftedRamaStats.txt");
    stringstream ss;
    map<string,bool>::iterator it;
    for (it = residueCombinations.begin();it != residueCombinations.end();it++){
	    vector<string> toks = MslTools::tokenize(it->first,"-");
	    cout << "PAIR: "<<it->first<<" has "<<counts[it->first]<<" counts."<<endl;
	    for (uint i = 0; i < opt.numSamples;i++){


		    // Get a random pair of Phi/Psi values
		    std::pair<double,double> angles = pps.getRandomPhiPsi(it->first);


		    // Add some noise into the angles.. (random angle 1/2 the width of a bin)
		    angles.first  += (rng.getRandomDouble() * (20/2)  * (rng.getRandomDouble() < 0.5 ? -1 : 1));
		    angles.second += (rng.getRandomDouble() * (20/2)  * (rng.getRandomDouble() < 0.5 ? -1 : 1));

		    while (abs(angles.first) > 180){
		      angles.first  = pps.getRandomPhiPsi(it->first).first + (rng.getRandomDouble() * (20/2)  * (rng.getRandomDouble() < 0.5 ? -1 : 1));
		    }

		    while (abs(angles.second) > 180){
		      angles.second = pps.getRandomPhiPsi(it->first).second + (rng.getRandomDouble() * (20/2)  * (rng.getRandomDouble() < 0.5 ? -1 : 1));
		    }


		    char tmp[300];
		    sprintf(tmp, "%3s %8.3f %3s %8.3f\n",
			    toks[0].c_str(),
			    angles.first,
			    toks[1].c_str(),
			    angles.second);
		    fout << tmp;
		    
	    }
    }

    fout.close();
    //cout << ss.str()<<endl;

	

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

    opt.numSamples = OP.getInt("numSamples");
    if (OP.fail()){
      opt.numSamples = 10000;
      cerr << "WARNING numSamples set to "<<opt.numSamples<<endl;
    }
    return opt;
}

