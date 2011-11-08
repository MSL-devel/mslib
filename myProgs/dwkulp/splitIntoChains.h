#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdblist");

	}




	// Storage for the vales of each option
	string pdblist;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

