#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("ref");
		
		// Rotamer library
		required.push_back("list");
		required.push_back("testPdb");

		// Max C-alpha clashes allowed between potential fusions
		optional.push_back("maxCaClashes");

		// Config file
	        defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string ref;
        string list;
        string testPdb;
        string configfile;
        double maxCaClashes;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

