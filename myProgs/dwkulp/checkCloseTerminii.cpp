#include <iostream>
#include <cstdlib>
#include <fstream>


#include "MslTools.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Timer.h"
#include "SasaCalculator.h"
#include "Position.h"
#include "PDBTopology.h"
#include "VectorPair.h"
#include "checkCloseTerminii.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <numeric>
#include <vector>

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("checkCloseTerminii");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
	// MSLOUT.turnOff("calcSasaInterface");

	Timer t;
	double start = t.getWallTime();

	// Read in the pdb list
	vector<string> list;
	MslTools::readTextFile(list, opt.list);


	// first string is list[i] name, second is posId, double is tmpSasa - sasaRef.
	for (uint i = 0; i < list.size();i++){

	  System tmp;
	  tmp.readPdb(list[i]);

	  Chain &ch = tmp.getChain("X");
	  double dist = ch.getPosition(0).getAtom("CA").distance(ch.getPosition(ch.positionSize()-1).getAtom("CA"));
	  if (dist < 15){
	    cout << list[i]<<endl;
	  }
	  
	}// END list


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
		cout << "calcSasaInterface --list LIST --dist DIST\n";
		exit(0);
	}
	opt.list= OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}

	return opt;
}


