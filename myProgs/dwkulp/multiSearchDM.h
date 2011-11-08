#include "MatrixWindow.h"

#include <string>
#include <vector>
using namespace std;
using namespace MSL;
// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("inputPDB");
		
		// List of PDBs to search against or binary
		optional.push_back("pdbList");
		optional.push_back("binaryDMs");
		
		// Search Criteria
		required.push_back("searchCriteria");

		// Size of window
		required.push_back("windowSize");
		
		// Number of searches
		optional.push_back("numberOfIterations");

		optional.push_back("allowIntraChainCompare");

		optional.push_back("likenessTolerance");
		
		optional.push_back("alignPdbs");
		optional.push_back("rmsdTol");

		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string inputPDB;
	string pdbList;
        string dmd;
	string searchCriteria;
	int numberOfIterations;
	int windowSize;	
	bool intraChainCompare;
	double likenessTolerance;
	double rmsdTol;
	bool alignPdbs;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

//gets RMSD for MatrixWindow win2 aligned onto win1. System sys2 corresponds to win2. i and j index the two windows in our list of distance matrices
void getRMSD(MatrixWindow *_win1, MatrixWindow *_win2, System *_sys2);
