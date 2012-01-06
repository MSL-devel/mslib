#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstdlib>

using namespace std;



struct DoFSpec {
	vector<string> atomNames;
        double minValue;
	double maxValue;
	double stepSize;
  DoFSpec(){
    minValue = 0.0;
    maxValue = 0.0;
    stepSize = 0.0;
    atomNames = vector<string>(4,"");
  }
        string toString(){
	  string outStr;
	  for (uint i = 0; i < atomNames.size();i++){
	    if (i == 0){
	      outStr = MslTools::stringf("%s",atomNames[i].c_str());
	    } else {
	      outStr = MslTools::stringf("%s %s",outStr.c_str(),atomNames[i].c_str());
	    }
	  }

	  outStr = MslTools::stringf("%s [ %8.3f - %8.3f , %8.3f ]",outStr.c_str(),minValue,maxValue,stepSize);
	  return outStr;
	}
};
struct DoF {
	vector<string> atomNames;
        double value;
};
// Define Input Options and Store them.

struct Options {

	// Set up options here...
	Options(){

		required.push_back("pdb");
		required.push_back("atomNames");
		required.push_back("range");
		optional.push_back("outPdb");
		optional.push_back("ref");
		optional.push_back("refSelect");
		optional.push_back("pdbSelect");
		optional.push_back("rmsd");
		optional.push_back("clashDist");


	}

	// Storage for the vales of each optional
	string pdb;
	string outPdb;	
        vector<DoFSpec> DoFs;
        string ref;
        string refSelect;
        string pdbSelect;
        double rmsd;
        double clashDist;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};



// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

