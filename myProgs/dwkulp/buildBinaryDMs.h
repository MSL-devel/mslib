#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// List of PDBs to search against
		required.push_back("pdbList");
		
		// Size of window
		required.push_back("windowSize");
		
		optional.push_back("allowIntraChainCompare");

		optional.push_back("likenessTolerance");
		optional.push_back("diagnolMWonly");
		optional.push_back("outName");
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdbList;
	int windowSize;	
	bool intraChainCompare;
	double likenessTolerance;
	string outName;
	bool diagnolMatrixWindowsOnly;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

