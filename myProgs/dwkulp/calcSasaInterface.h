#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input list
		required.push_back("list");
		optional.push_back("dist");
		optional.push_back("resfile");
		

	}

	// Storage for the vales of each optional
	string list;
        double dist;
        bool resfile;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
double mean(vector<double> &data);
double stddev(vector<double> &data);
