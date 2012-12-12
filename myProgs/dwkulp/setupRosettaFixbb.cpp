#include <iostream>
#include <cstdlib>
#include <fstream>
#include "OptionParser.h"
#include "FastaReader.h"
#include "System.h"
#include "MslOut.h"
#include "PolymerSequence.h"

#include "setupRosettaFixbb.h"


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("setupRosettaFixbb");



int main(int argc, char *argv[]) {
  Options opt = setupOptions(argc,argv);

  System sys;
  sys.readPdb(opt.pdb);

  // Create new files
  string basefile = MslTools::stringf("%s_%s", MslTools::getFileName(opt.pdb).c_str(),opt.output.c_str());
  string resfile = MslTools::stringf("%s.resfile", basefile.c_str());

  // Print out the files
  fstream fresfile;
  fresfile.open(resfile.c_str(),std::ios::out);
  fresfile << "NATRO" <<endl;
  fresfile << "start"<<endl;
  map<string,bool> uniquePositions;
  fresfile <<endl<< "# Mutated positions "<<endl;
  for (uint i = 0; i < opt.mutate.size();i++){

    // Parse mutate format L260A or L_260_A
    vector<string> mutation = parseMutateString(opt.mutate[i]);

    if (!sys.positionExists(mutation[0])){
      cerr << "ERROR 2222 Position: "<<mutation[0]<< " does not exist in "<<opt.pdb<<endl;
      exit(2222);
    }
    
    // Resfile PIKAA
    fresfile << MslTools::stringf("%d%1s %1s PIKAA %s EX 1 LEVEL 4 EX 2 LEVEL 4 EX 3 EX 4\n",sys.getPosition(mutation[0]).getResidueNumber(),sys.getPosition(mutation[0]).getResidueIcode().c_str(),sys.getPosition(mutation[0]).getChainId().c_str(),mutation[1].c_str());
    uniquePositions[sys.getPosition(mutation[0]).getPositionId()] = 1;
  }

  for (uint i = 0; i < opt.mutate.size();i++){

    // Parse mutate format L260A or L_260_A
    vector<string> mutation = parseMutateString(opt.mutate[i]);

    if (!sys.positionExists(mutation[0])){
      cerr << "ERROR 2222 Position: "<<mutation[0]<< " does not exist in "<<opt.pdb<<endl;
      exit(2222);
    }
    if (opt.neighbor_dist != 0.0){
      vector<int> neighbors = sys.getPosition(mutation[0]).getCurrentIdentity().findNeighbors(opt.neighbor_dist);
      fresfile <<endl<< "# Neighbors of "<<mutation[0]<<endl;
      for (uint n = 0; n < neighbors.size();n++){
	map<string,bool>::iterator it;
	it = uniquePositions.find(sys.getPosition(neighbors[n]).getPositionId());
	if (it == uniquePositions.end()){
	  fresfile << MslTools::stringf("%d%1s %1s NATAA EX 1 LEVEL 4 EX 2 LEVEL 4 EX 3 EX 4\n",sys.getPosition(neighbors[n]).getResidueNumber(),sys.getPosition(neighbors[n]).getResidueIcode().c_str(),sys.getPosition(neighbors[n]).getChainId().c_str());	
	}
	uniquePositions[sys.getPosition(neighbors[n]).getPositionId()] = true;
      }
      fresfile << endl << endl;
    }
  }
  fresfile.close();
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
		cout << "setupRosettaMSA --msa msa.fasta --pdb PDB\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.mutate = OP.getMultiString("mutate");
	if (OP.fail()){
	  cerr << "ERROR 111 mutate not specified\n";
	  exit(1111);
	}
	opt.neighbor_dist = OP.getDouble("neighbor_dist");
	if (OP.fail()){
	  opt.neighbor_dist = 0.0;
	}

	
	opt.output = OP.getString("output");
	if (OP.fail()){
	  opt.output = "mutate";
	}
	return opt;
}


vector<string> parseMutateString(string _mutation){

  vector<string> result;

  vector<string> parseUnderscore = MslTools::tokenize(_mutation,":");

  // FORMAT: A,260C:M   means chain A at position 260C should mutate to Met
  if (parseUnderscore.size() == 2){
    for (uint i = 0; i < parseUnderscore.size();i++){
      result.push_back(parseUnderscore[i]);
    }
  }


  return result;

}
