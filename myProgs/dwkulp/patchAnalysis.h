#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input list
		required.push_back("list");
		required.push_back("pdb");
		required.push_back("patch1");
		required.push_back("patch2");

	}

	// Storage for the vales of each optional
	string list;
	string pdb;
	string patch1;
	string patch2;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

