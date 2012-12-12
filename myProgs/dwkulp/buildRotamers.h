#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb");
		
		// Rotamer library
		required.push_back("rotlib");
		required.push_back("identity");
		required.push_back("position");

		// Number of rotamers
		optional.push_back("numRotamers");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdb;
	string rotlib;	
	vector<string> identity;	
	string position;	
	int numRotamers;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

