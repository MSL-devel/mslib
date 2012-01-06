#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdblist");
		required.push_back("minLoopLength");
		required.push_back("maxLoopLength");
		required.push_back("startResidue");
		required.push_back("endResidue");
		required.push_back("naturalBreaks");
		required.push_back("numCaBreaksAllowed");

	}




	// Storage for the vales of each option
	string pdblist;
        int loopMin;
        int loopMax;
	string startResidue;
	string endResidue;
        int numCaBreaksAllowed;
        vector<string> naturalBreaks;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

