#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("pdb");
		required.push_back("designOptions");
		optional.push_back("designResidues");

	}




	// Storage for the vales of each option
	string pdb;
        string designOptions;
        string designResidues;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
