#include <vector>
#include "PyMolVisualization.h"
struct Options {

	// Set up options here...
	Options(){

		// PDB
		required.push_back("pdb");
		optional.push_back("prosite");
		
		optional.push_back("topfile");
		optional.push_back("parfile");
		optional.push_back("dielectric");
		optional.push_back("distanceDependentElectrostatics");
		optional.push_back("vdwScale");
		optional.push_back("hbondfile");
		optional.push_back("baselinefile");
		optional.push_back("cuton");
		optional.push_back("cutoff");
		
		optional.push_back("ref");
	}




	// Storage for the vales of each option
	string pdb;
        string ref;
        string prosite;
        string topfile;
        string parfile;
        string hbondfile;
        string baselinefile;
        double dielectric;
        bool distanceDependentElectrostatics;
        double vdwScale;
        double cuton;
        double cutoff;
        PyMolVisualization pymol;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);
void reportPositionMetrics(System &_sys, Options &_opt);
void reportStructureMotifMetrics(System &_sys, Options &_opt);
void reportSequenceMotifMetrics(System &_sys, Options &_opt);
void reportEnergyMetrics(System &_sys, Options &_opt);
