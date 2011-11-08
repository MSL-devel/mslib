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

		// List of PDB files
		required.push_back("list");
		required.push_back("vh");

		// Number of rotamers
		optional.push_back("numRotamers");
		optional.push_back("rmsd");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdb;
	string rotlib;	
        string list;
        string vh;
	int numRotamers;
        double rmsd;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

