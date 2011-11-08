#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb");
		required.push_back("outPdb");
		

		// Config file
	        defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdb;
        string outPdb;
        string configfile;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

