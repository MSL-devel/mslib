#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){


		// The input PDB file..
		required.push_back("pdb");

		// Optional arguments
		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("rotlib");
	}

	// Storage for the vales of each option
	string pdb;
	string topfile;
	string parfile;
	string rotlib;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
