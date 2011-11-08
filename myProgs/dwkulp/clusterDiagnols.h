#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("pdb_db");
		
		// List of PDBs to search against or binary
		optional.push_back("pymol");
		optional.push_back("numDMs");
		
	}

	// Storage for the vales of each optional
        string pdb_binary_database;
        int numDMs;
        bool pymol;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
vector<int> getClusters(vector<vector<double> > &_distMatrix) ;
