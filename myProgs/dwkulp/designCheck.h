#include <vector>
#include "PyMolVisualization.h"
struct Options {

	// Set up options here...
	Options(){

		// PDB
		required.push_back("pdb");
		required.push_back("prosite");
		optional.push_back("ref");
	}




	// Storage for the vales of each option
	string pdb;
        string ref;
        string prosite;
        PyMolVisualization pymol;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
void reportPositionMetrics(System &_sys, Options &_opt);
void reportStructureMotifMetrics(System &_sys, Options &_opt);
void reportSequenceMotifMetrics(System &_sys, Options &_opt);

