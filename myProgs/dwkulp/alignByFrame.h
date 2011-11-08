#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdblist");
		
		// Selection
		required.push_back("sel");

		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdblist;
	string sel;	
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

