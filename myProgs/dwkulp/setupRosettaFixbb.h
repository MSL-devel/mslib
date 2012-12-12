#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){


		required.push_back("pdb");
		required.push_back("mutate");


		optional.push_back("neighbor_dist");
		optional.push_back("debug");
		optional.push_back("output");
	}

	// Storage for the vales of each optional
	string pdb;
	vector<string> mutate;	
        double neighbor_dist;
        string output;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
vector<string> parseMutateString(string _mutation);
