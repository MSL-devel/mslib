#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("template");
		required.push_back("templateStem1");
		required.push_back("templateStem2");
		required.push_back("fragment");


	}




	// Storage for the vales of each option
	string templatePDB;
        string templateStem1;
        string templateStem2;
        string fragmentPDB;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
