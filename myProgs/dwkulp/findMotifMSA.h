#include <string>
#include <vector>
using namespace std;



// Define Input Options and Store them.
struct Options {

        enum vtype { freq=1,entropy=2,mostfreq=3};

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");

		// Input fasta
		required.push_back("msa");

		// Regex
		required.push_back("regex");

	}

	// Storage for the vales of each optional
	string pdb;
	string msa;	
        string regex;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

