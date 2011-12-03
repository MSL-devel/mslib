#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Scaffold with epitope
		required.push_back("scaffold");

		// Reference pdb
		required.push_back("ref");
		
		// Rotamer library
		required.push_back("rotlib");

		// Number of rotamers
		optional.push_back("numRotamers");
		optional.push_back("scaffoldSelect");
		optional.push_back("rotamerSelect");
		optional.push_back("glycanSelect");
		optional.push_back("numGlycanClashesAllowed");
		optional.push_back("rmsd");
		optional.push_back("outFile");
		optional.push_back("debug");


		// Config file
	        defaultArgs.push_back("configfile");
	}

	// Storage for the vales of each optional
	string scaffold;
	string epitopeSelect;

	string ref;
        string rotamerSelect;

	string rotlib;	
	int numRotamers;
	string glycanSelect;
	int numGlycanClashesAllowed;
        double rmsd;

	bool debug;

        string configfile;
        string outFile;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

