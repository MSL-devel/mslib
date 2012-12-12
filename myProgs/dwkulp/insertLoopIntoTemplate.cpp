#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "AtomContainer.h"
#include "FuseChains.h"
#include "insertLoopIntoTemplate.h"

using namespace MSL;
using namespace std;
using namespace MslTools;

int main(int argc, char *argv[]) {

        // Parse command line options into Options structure
	Options opt = setupOptions(argc, argv);

	// Template is the structure to insert into
	System templatePDB;
	templatePDB.readPdb(opt.templatePDB);
	if (templatePDB.size() > 1 ){
	  cerr << "ERROR 3333 currently insertLoopIntoTemplate assumes templatePDB has a single chain, I found "<<templatePDB.size()<<" chains."<<endl;
	  exit(3333);
	}

	// Fragment is the structure to take peice of structure from
	System fragmentPDB;
	fragmentPDB.readPdb(opt.fragmentPDB);

	// Find the shortest chain to use as a fragment
	int shortestChain = 0;
	for (uint c = 1; c < fragmentPDB.size();c++){
	  if (opt.fragmentChain == fragmentPDB.getChain(c).getChainId()){
	    shortestChain = c;
	    break;
	  } 

	  if (fragmentPDB.getChain(c).size() < fragmentPDB.getChain(shortestChain).size()){
	    shortestChain = c;
	  }

	}

	Chain &fragChain = fragmentPDB.getChain(shortestChain);
	cout << "Fragment Chain is "<<fragChain.getChainId()<<" from "<<opt.fragmentPDB<<endl;

	AtomContainer fusedProtein;
	FuseChains fuse;
	fusedProtein.addAtoms(fuse.fuseInsert(templatePDB.getChain(0), fragChain,opt.templateStem1,opt.templateStem2));
	
	
	// Add additional chains from the fragmentPDB
	for (uint i = 0; i < fragmentPDB.size();i++){
	    if (i == shortestChain) continue;
	    fprintf(stdout, "Adding chain %1s from fragmentPDB\n",fragmentPDB.getChain(i).getChainId().c_str());
	    fusedProtein.addAtoms(fragmentPDB.getChain(i).getAtomPointers());
	}

	int numClashes = 0;
	if (opt.clashCheck){
	  for (uint i = 0; i < fusedProtein.size();i++){
	    if (fusedProtein[i].getName() != "CA") continue;

	    for (uint j = i+1;j < fusedProtein.size();j++){
	      if (fusedProtein[j].getName() != "CA") continue;

	      if (fusedProtein[i].distance(fusedProtein[j]) < 2.5){
		numClashes++;
	      }

	    }
	  }

	}

	int numClashTolerance = 3;
	if (numClashes <= 3) {
	  char fname[200];
	  sprintf(fname, "%s_%s.pdb",MslTools::getFileName(opt.templatePDB).c_str(),MslTools::getFileName(opt.fragmentPDB).c_str());
	  cout << "Fused Protein: "<< fname << " number of clashes: "<<numClashes<<endl;
	  fusedProtein.writePdb(fname);
	} else {
	  fprintf(stdout, " ERROR fusion of %s with %s has %d clashes which is more than the tolerance of %d clashes, no file output\n",opt.templatePDB.c_str(),opt.fragmentPDB.c_str(),numClashes,numClashTolerance);
	}

}
	

Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;

	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage: insertSelectionIntoTemplate " << endl;
		cout << endl;
		cout << "\n";
		cout << "template PDB\n";
		cout << "templateStem1 positionId\n";
		cout << "templateStem2 positionId\n";
		cout << "fragment PDB\n";
		cout << endl;
		exit(0);
	}

	opt.templatePDB = OP.getString("template");
	if (OP.fail()){
		cerr << "ERROR 1111 no template specified."<<endl;
		exit(1111);
	}


	opt.fragmentPDB = OP.getString("fragment");
	if (OP.fail()){
		cerr << "ERROR 1111 no fragment specified."<<endl;
		exit(1111);
	}
	opt.fragmentChain = OP.getString("fragmentChain");
	if (OP.fail()){
	  opt.fragmentChain = "NO_CHAIN_INPUT_USE_SHORTEST_CHAIN_IN_FRAGMENT_PDB";
	}
	opt.templateStem1= OP.getString("templateStem1");
	opt.templateStem2 = OP.getString("templateStem2");


	opt.clashCheck = OP.getBool("clashCheck");
	
	
	return opt;
}
