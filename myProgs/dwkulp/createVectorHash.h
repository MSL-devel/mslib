#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// List of PDB files
		required.push_back("list");

		// Number of rotamers
		optional.push_back("outfile");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
        string list;
        string outfile;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

