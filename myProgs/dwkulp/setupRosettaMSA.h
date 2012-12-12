#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){


		required.push_back("pdb");
		required.push_back("msa");

		optional.push_back("select_seqs");
		optional.push_back("remodel_neighbors");
		optional.push_back("skipBlankEnds");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdb;
	string msa;	
        int remodel_neighbors;
        vector<string> select_seqs;
	bool skipBlankEnds;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

