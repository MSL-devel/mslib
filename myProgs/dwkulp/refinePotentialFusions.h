#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){


		required.push_back("pdb");
		required.push_back("junction_positionABC_1");
		required.push_back("junction_positionABC_2");
		required.push_back("junction_positionXYZ_1");
		required.push_back("junction_positionXYZ_2");
		required.push_back("min_junction_distance");
		required.push_back("max_junction_distance");
		required.push_back("transZ_limit");
		required.push_back("rotZ_limit");
	}


	// Storage for the vales of each option
	string pdb;
	string junction_positionABC_1;
	string junction_positionABC_2;
	string junction_positionXYZ_1;
	string junction_positionXYZ_2;
	double min_junction_distance;
	double max_junction_distance;
        double transZ;
        double rotZ;
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

