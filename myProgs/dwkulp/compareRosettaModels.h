#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){


		required.push_back("pdb1");
		required.push_back("pdb2");
		required.push_back("resfile");
		required.push_back("out");

	}

	// Storage for the vales of each optional
	string pdb1;
	string pdb2;
        string resfile;
        string out;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
void parseResfile(vector<string> &_resfilePositions, string _filename);

