#include <vector>
#include <string>
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// Renumber pdb by refpdb
		required.push_back("pdblist");
		required.push_back("distanceStem1");
		required.push_back("distanceStem2");

		required.push_back("loop1min");
		required.push_back("loop1max");
		required.push_back("loop2min");
		required.push_back("loop2max");

		required.push_back("stem1pdb");
		required.push_back("stem1res1");
		required.push_back("stem1res2");

		required.push_back("stem2pdb");
		required.push_back("stem2res1");
		required.push_back("stem2res2");

		optional.push_back("rmsdTol");
	}




	// Storage for the vales of each option
	string pdblist;
        double distanceStem1;
        double distanceStem2;
        int loop1min;
        int loop1max;
        int loop2min;
        int loop2max;
        string stem1pdb;
        string stem1res1;
        string stem1res2;

        string stem2pdb;
        string stem2res1;
        string stem2res2;

        double rmsdTol;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
