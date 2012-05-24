#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		required.push_back("pdb");
		required.push_back("newStartRes");
		required.push_back("newEndRes");
		optional.push_back("outPdbp");

	}

	// Storage for the vales of each optional
        string pdb;
        string outPdb;
        string newStartRes;
        string newEndRes;


	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

