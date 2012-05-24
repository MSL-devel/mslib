#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdblist");
		optional.push_back("maxCbCb");
	}




	// Storage for the vales of each option
	string pdblist;
        double maxCbCb;
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

