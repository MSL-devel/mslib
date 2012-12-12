#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		required.push_back("pdb");
		required.push_back("rotlib");
		required.push_back("topfile");
		required.push_back("parfile");
		required.push_back("variablePosition");
		required.push_back("variableRotamer");
		required.push_back("fixedPosition");
		required.push_back("znPosition");
		optional.push_back("numModels");
		optional.push_back("configfile");


	}

	// Storage for the vales of each optional
        string configFile;
        vector<string> defaultArgs;
	string pdb;
	string rotlib;	
        string topfile;
        string parfile;
	std::vector<std::string> variable_positions;
        std::vector<std::vector<std::string> > variable_identities;
        std::vector<std::vector<int> > variable_rotamers;

	std::vector<std::string> fixed_positions;
	std::vector<std::string> zn_positions;

        int numModels;
  
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

bool checkMetalSiteGeometry(Position &_fix, Position &_var, Atom &_zn);
void cleanExit(int sig) ;
