#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Type of symmetry to apply
		required.push_back("symmetry");

		// super-helical Radius parameters (LOW, HIGH, STEP)
		required.push_back("superHelicalRadius");

		// alpha-helical phase angle parameters (LOW, HIGH, STEP)
		required.push_back("alphaHelicalPhaseAngle");

		// super-helical pitch angle (LOW, HIGH, STEP)
		optional.push_back("superHelicalPitchAngle");

		// Number of residues
		required.push_back("numberOfResidues");

		// Optional arguments
		optional.push_back("d2zTranslation");
		optional.push_back("superHelicalPhaseAngle");
		optional.push_back("name");


	}

	// Storage for the vales of each option
	string symmetry;
	vector<double> superHelicalRadius;
	vector<double> alphaHelicalPhaseAngle;
	int numberOfResidues;

	vector<double> d2zTranslation;
	vector<double> superHelicalPhaseAngle;
        vector<double> superHelicalPitchAngle;
	string name;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
