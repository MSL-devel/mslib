
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "PhiPsiStatistics.h"
#include "clashCheck.h"
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

        map<string,int> clashMap;
        int total_clashes = 0;

	// Each chain
	for (uint c1 = 0; c1 < sys.chainSize();c1++){

	  for (uint p1 = 0; p1 < sys.getChain(c1).positionSize();p1++){

	    int pos_clashes = 0;
	    for (uint p2 = 0; p2 < sys.getChain(c1).positionSize();p2++){
	      if (abs (p1 - p2) < 3) continue;

	      for (uint a1 = 0; a1 < sys.getChain(c1).getPosition(p1).atomSize();a1++){
		Atom &at1 = sys.getChain(c1).getPosition(p1).getAtom(a1);
		if (at1.getElement() == "H") continue;
		for (uint a2 = a1+1; a2 < sys.getChain(c1).getPosition(p2).atomSize();a2++){
		  Atom &at2 = sys.getChain(c1).getPosition(p2).getAtom(a2);
		  if (at2.getElement() == "H") continue;

		  double dist = at1.distance(at2);
		  if (dist < 1.5){
		    pos_clashes++;
		  }
		  
		  
		}
	      }


	    }

	    clashMap[sys.getChain(c1).getPosition(p1).getPositionId()] = pos_clashes;
	    total_clashes += pos_clashes;

	  }


	  
	}


    
	map<string,int>::iterator it;
	for (it = clashMap.begin(); it != clashMap.end();it++){
	  fprintf(stdout, "%10s %8d\n",it->first.c_str(),it->second);
	}
	fprintf(stdout, "TOTAL: %s %8d\n",lines[i].c_str(),total_clashes);

	if (total_clashes < opt.tooManyClashes) {
	  fprintf(stdout,"MIN CLASHES: %s",lines[i].c_str());
	}
	
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
	cout << "Usage: clashCheck " << endl;
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
    opt.tooManyClashes = OP.getInt("tooManyClashes");
    if (OP.fail()){
	opt.tooManyClashes = 0;
    }
    
    return opt;
}
