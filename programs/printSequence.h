#include <vector>
struct Options {

	// Set up options here...
	Options(){

		// Energy Table
		required.push_back("pdb");

	}




	// Storage for the vales of each option
	string pdb;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

