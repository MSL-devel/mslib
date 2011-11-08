#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("list");
		required.push_back("residueType");
		required.push_back("separateByResidueType");

		// Number of rotamers
		optional.push_back("debug");
		optional.push_back("alignAtoms");

	}

	// Storage for the vales of each optional
	string list;
	vector<string> residueType;
        string alignAtoms;
	bool separateByResidueType;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

