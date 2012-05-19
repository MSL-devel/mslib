#include <string>
#include <vector>

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){


		// The input PDB file..
		required.push_back("pdb");

		// Optional arguments
		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("rotlib");
	}

	// Storage for the vales of each option
	std::string pdb;
	std::string topfile;
	std::string parfile;
	std::string rotlib;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);