#include <string>
#include <vector>

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
	std::string symmetry;
	std::vector<double> superHelicalRadius;
	std::vector<double> alphaHelicalPhaseAngle;
	int numberOfResidues;

	std::vector<double> d2zTranslation;
	std::vector<double> superHelicalPhaseAngle;
        std::vector<double> superHelicalPitchAngle;
	std::string name;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
