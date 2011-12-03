#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input list
		required.push_back("list");
		required.push_back("fasta");
		required.push_back("refSeqName");
		

	}

	// Storage for the vales of each optional
	string list;
        string refSeqName;
        string fasta;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

