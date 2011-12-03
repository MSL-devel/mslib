#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

	
		required.push_back("complex");
		required.push_back("ref");
		required.push_back("chainInComplexToAlign");
		required.push_back("chainInComplexToCheck");
		required.push_back("chainInRef");

		optional.push_back("domain");
		optional.push_back("debug");

		// Config file
	        defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string complex;
	string ref;
        string chainToAlign;
        vector<string> chainsToCheck;
        string chainInRef;
        vector<vector<string> > domains;
	bool debug;
        string configfile;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

