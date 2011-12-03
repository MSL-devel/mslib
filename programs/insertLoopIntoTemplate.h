#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("template");
		required.push_back("fragment");

		optional.push_back("templateStem1");
		optional.push_back("templateStem2");
		optional.push_back("clashCheck");
		optional.push_back("numClashes");
		optional.push_back("includeTemplateStems");
		optional.push_back("checkCaCaDistances");
	}




	// Storage for the vales of each option
	string templatePDB;
        string templateStem1;
        string templateStem2;
        string fragmentPDB;
        bool clashCheck;
        bool includeTemplateStems;
        bool checkCaCaDistances;
	int numClashes;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
