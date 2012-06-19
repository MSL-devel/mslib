#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB 
		required.push_back("pdb");

		required.push_back("position1");
		required.push_back("position2");

		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("hbondfile");
		optional.push_back("dielectric");
		optional.push_back("distanceDependentElectrostatics");
		optional.push_back("ctonnb");
		optional.push_back("ctofnb");
		optional.push_back("cutnb");
		optional.push_back("cuthb");


	}




	// Storage for the vales of each option
	string pdb;
        string position1;
        string position2;
        string topfile;
        string parfile;
        string hbondfile;
        double ctonnb;
        double ctofnb;
        double cutnb;
        double cuthb;
        double dielectric;
        double distanceDependentElectrostatics;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
