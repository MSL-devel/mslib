#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb");
		
		// Rotamer library
		required.push_back("insertionPoint1domain1");
		required.push_back("insertionPoint1domain2");
		required.push_back("insertionPoint2domain1");
		required.push_back("insertionPoint2domain2");
		required.push_back("linker1lengths");
		required.push_back("linker2lengths");

		optional.push_back("designInterfaceDistance");
		optional.push_back("extraDesignPositions");
		optional.push_back("extraRemodelPositions");
		// Config file
	        defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdb;
        string insertionPoint1domain1;
        string insertionPoint1domain2;
        string insertionPoint2domain1;
        string insertionPoint2domain2;
        string linker1lengths;
        string linker2lengths;
        string configfile;
        double designInterfaceDistance;
        string extraDesignPositions;
        string extraRemodelPositions;
	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

