#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input list
		required.push_back("pdb");
		required.push_back("domains");

	}

	// Storage for the vales of each optional
	string pdb;
        vector<string> domains;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

