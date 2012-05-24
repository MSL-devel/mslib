#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("template");
		optional.push_back("templateChain");
		required.push_back("fragment");
		required.push_back("fragmentChain");

		optional.push_back("templateStem1");
		optional.push_back("templateStem2");
		optional.push_back("clashCheck");
		optional.push_back("numClashes");
		optional.push_back("includeTemplateStems");
		optional.push_back("checkCaCaDistances");
		optional.push_back("outputRosettaFiles");
	}




	// Storage for the vales of each option
	string templatePDB;
	string templateChain;
        string templateStem1;
        string templateStem2;
        string fragmentPDB;
        bool clashCheck;
        bool includeTemplateStems;
        bool checkCaCaDistances;
	int numClashes;
        string fragmentChain;
        string outputRosettaFiles;
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
