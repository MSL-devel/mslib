
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
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

	bool brokenChain = false;

	
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
	cout << "Usage: getConsecutivePhiPsi " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdblist PDB\n";
	cout << "SSE L # L = Loop, S = Strand, H = Helix\n"
	cout << endl;
	exit(0);
    }

    opt.pdblist = OP.getString("pdblist");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdblist specified."<<endl;
	exit(1111);
    }

    opt.SSE = OP.getString("sse");
    if (OP.fail()){
	opt.SSE = "";
    }


    return opt;
}
