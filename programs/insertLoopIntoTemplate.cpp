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
	if (templatePDB.chainSize() > 1 ){
	  cerr << "ERROR 3333 currently insertLoopIntoTemplate assumes templatePDB has a single chain, I found "<<templatePDB.chainSize()<<" chains."<<endl;
	  exit(3333);
	}

	cout << "Template chain "<<templatePDB.getChain(0).getChainId()<<" has "<<templatePDB.getChain(0).positionSize()<< " residues."<<endl;
 
	// Fragment is the structure to take peice of structure from
	System fragmentPDB;
	fragmentPDB.readPdb(opt.fragmentPDB);

	// Find the shortest chain to use as a fragment
	int shortestChain = 0;
	for (uint c = 1; c < fragmentPDB.chainSize();c++){
	  if (fragmentPDB.getChain(c).positionSize() < fragmentPDB.getChain(shortestChain).positionSize()){
	    shortestChain = c;
	  }
	}

	Chain &fragChain = fragmentPDB.getChain(shortestChain);
	cout << "Fragment Chain is "<<fragChain.getChainId()<<" from "<<opt.fragmentPDB<< " "<<fragChain.positionSize()<< " residues"<<endl;

	AtomContainer fusedProtein;
	FuseChains fuse;
	fusedProtein.addAtoms(fuse.fuseInsert(templatePDB.getChain(0), fragChain,opt.templateStem1,opt.templateStem2,opt.includeTemplateStems));
	

	if (opt.checkCaCaDistances){

	  bool chainOk = true;
	  int chainBreak1 = 0;
	  int chainBreak2 = 0;
	  for (uint i = 0; i < fusedProtein.size()-1;i++){
	    if (fusedProtein[i].getName() != "CA") continue;

	    int nextCa = 0;
	    for (uint j = i+1; j < fusedProtein.size()-1;j++){
	      if (fusedProtein[j].getName() == "CA") { nextCa = j; break; }
	    }

	    if (fusedProtein[i].distance(fusedProtein[nextCa]) > 4.0){
		chainOk = false;
		chainBreak1 = i;
		chainBreak2 = nextCa;
		break;
	      }
	  }

	  if (!chainOk){
	    fprintf(stdout, "ERROR fusion of %s with %s has chain break at %6s -> %6s, no file output\n",opt.templatePDB.c_str(),opt.fragmentPDB.c_str(),fusedProtein[chainBreak1].getPositionId().c_str(),fusedProtein[chainBreak2].getPositionId().c_str());
	    exit(2543);
	  }

	}

	// Add additional chains from the fragmentPDB
	int chainIdIndex = 0;
	for (uint i = 0; i < fragmentPDB.chainSize();i++){
	    if (i == shortestChain) continue;

	    string chains = "BCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	    string newChainId = chains.substr(chainIdIndex,1);	
	    while (fragmentPDB.chainExists(newChainId)){

	        chainIdIndex++;
	        newChainId = chains.substr(chainIdIndex,1);	

		if (chainIdIndex > chains.length()){
		  cerr << "ERROR 5555 insertLoopIntoTemplate couldn't find a free chain id to use for fragmentPDB chain "<< fragmentPDB.getChain(i).getChainId()<<endl;
		  exit(5555);
		}
	    } 

	    fprintf(stdout, "Adding chain %1s from fragmentPDB as chain %s\n",fragmentPDB.getChain(i).getChainId().c_str(),newChainId.c_str());

	    // Add each atom and change its chain id.
	    for (uint a = 0; a < fragmentPDB.getChain(i).getAtomPointers().size();a++){
		      fusedProtein.addAtom(fragmentPDB.getChain(i).getAtom(a));
		      fusedProtein.getAtom(fusedProtein.size()-1).setChainId(newChainId);
	    }

	    // Increment the chain index
	    chainIdIndex++;
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

	int numClashTolerance = opt.numClashes;
	if (numClashes <= opt.numClashes) {
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

	opt.templateStem1= OP.getString("templateStem1");
	opt.templateStem2 = OP.getString("templateStem2");

	opt.includeTemplateStems = OP.getBool("includeTemplateStems");

	opt.clashCheck = OP.getBool("clashCheck");
	opt.numClashes = OP.getInt("numClashes");	
	opt.checkCaCaDistances = OP.getBool("checkCaCaDistances");
	return opt;
}
