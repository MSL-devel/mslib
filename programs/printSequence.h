#include <vector>
struct Options {

	// Set up options here...
	Options(){

		// Energy Table
		required.push_back("pdb");
		optional.push_back("fasta");

	}




	// Storage for the vales of each option
	string pdb;
        bool fasta;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

