#include <string>
#include <vector>

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("pdb");
		required.push_back("x");
		required.push_back("y");
		required.push_back("z");
		optional.push_back("sele");
		optional.push_back("outfile");

	}

	// Storage for the vales of each optional
	std::string pdb;
        double x;
        double y;
        double z;
	std::string sele;
	std::string outfile;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

