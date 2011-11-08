
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "getAAStats.h"
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

    map<string,int> aaMap;
    map<string,int> aaH;
    map<string,int> aaL;
    map<string,int> aaE;
    for (uint i = 0; i < lines.size();i++){

	System sys;
	sys.readPdb(lines[i]);

	
	for (uint p = 0; p < sys.positionSize();p++){
	  
	  Residue &r = sys.getPosition(p).getCurrentIdentity();
	  if (!r.atomExists("CA")) continue;

	  aaMap[r.getResidueName()]++;

	  if (r("CA").getSegID() == "HHHH") aaH[r.getResidueName()]++;
	  if (r("CA").getSegID() == "LLLL") aaL[r.getResidueName()]++;
	  if (r("CA").getSegID() == "EEEE") aaE[r.getResidueName()]++;

	}
    }

    map<string,int>::iterator it;
    fprintf(stdout,"PDB list=%-s\n",opt.pdblist.c_str());
    fprintf(stdout,"%3s %-20s %-20s %-20s %-15s\n","RES","HELIX","STRAND","LOOP","TOTAL");
    for (it = aaMap.begin();it != aaMap.end();it++){

      map<string,int>::iterator itH;
      map<string,int>::iterator itL;
      map<string,int>::iterator itE;

      int numH = 0;
      int numL = 0;
      int numE = 0;
      itH = aaH.find(it->first);
      if (itH != aaH.end()) numH = itH->second;
      itL = aaL.find(it->first);
      if (itL != aaL.end()) numL = itL->second;
      itE = aaE.find(it->first);
      if (itE != aaE.end()) numE = itE->second;

      fprintf(stdout,"%3s %-15d (%3.0f) %-15d (%3.0f) %-15d (%3.0f) %-15d\n",it->first.c_str(),numH,(double)numH/(double)it->second*100,numE,(double)numE/(double)it->second*100,numL,(double)numL/(double)it->second*100,it->second);
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
	cout << "Usage: getAAstats " << endl;
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
    return opt;
}
