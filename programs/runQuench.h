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
		optional.push_back("outfile");
		optional.push_back("positions");
		optional.push_back("largeRotNum");
		optional.push_back("smallRotNum");
	}

	// Storage for the vales of each option
	string pdb;
	string topfile;
	string parfile;
	string rotlib;
	string outfile;
	vector<string> positions;
	int largeRotNum;
	int smallRotNum;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
