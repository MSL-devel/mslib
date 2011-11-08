
#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "PDBWriter.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "VectorHashing.h"

#include "createVectorHash.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("createVectorHash");

int main(int argc, char *argv[]) {


	// Parse commandline options
	Options opt = setupOptions(argc,argv);

	cout << "READ LIST"<<endl;
	vector<string> pdbs;  
	ifstream fs;

	fs.open(opt.list.c_str());
	if (fs.fail()){
		cerr<<"Cannot open file "<<opt.list<<endl;
		exit(1);
	}

	while(true){
		string line;
		getline(fs, line);

		if(fs.fail()){
			//no more lines to read, quite the while.
			break;
		}

		if(line==""){
			continue;
		}
		pdbs.push_back(line);
	}

	fs.close();


	// Read a list of PDBs into a single atom vector.
	Timer t;
	VectorHashing testHash;
	MSLOUT.turnOff("VectorHashing");
	for (uint i = 0; i < pdbs.size();i++){

	  double startHashTime = t.getWallTime();
	  System testSys;
	  //testSys.readPdb("/Users/dwkulp/work/DandD/test_section.pdb");
	  //testSys.readPdb("/Users/dwkulp/work/VaccineDesign_PGT128/tertFragSearch/pgt128_noGlycans.pdb");
	  testSys.readPdb(pdbs[i]);


	  MSLOUT.stream() << "Build VectorHash from PDB: "<<pdbs[i]<<"\n";
	  testHash.addToVectorHash(testSys,pdbs[i]);
	  MSLOUT.stream() << "done building hash\n";
	  double endTimeHash = t.getWallTime();
	  MSLOUT.fprintf(stdout,"Time %8.3f for building hash with test system\n",(endTimeHash - startHashTime));

	}


	testHash.save_checkpoint(opt.outfile);
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
		cout << "createVectorHash --list pdb.list --outfile foo.vh\n";
		exit(0);
	}
	opt.list = OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
	  cerr << "WARNING outfile not specified. Using: '"<<opt.list<<".vh'"<<endl;
	  opt.outfile = opt.list+".vh";

	}

	return opt;
}
