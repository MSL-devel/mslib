#include <string>
#include <vector>
using namespace std;



// Define Input Options and Store them.
struct Options {

        enum vtype { freq=1,entropy=2,mostfreq=3};

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");

		// Input selection (part of pdb to score)
		required.push_back("sel");
		
		// Input fasta
		required.push_back("fasta");

		// Input reference counts
		required.push_back("refAACounts");


		// Number of rotamers
		optional.push_back("refSeqName");
		optional.push_back("refSeqOffset");
		optional.push_back("outPdb");
		optional.push_back("logodds");
		optional.push_back("applyToAllChains");
		optional.push_back("regexForFasta");
		optional.push_back("debug");
		optional.push_back("seq");
		optional.push_back("type");
		defaultArgs.push_back("configfile");


		valueMap["freq"] = freq;
		valueMap["entropy"] = entropy;
		valueMap["mostfreq"] = mostfreq;
	}

	// Storage for the vales of each optional
	string pdb;
	string sel;	
	string fasta;	
	string refAACounts;	
	string refSeqName;	
        int refSeqOffset;
	string outPdb;	
        string configfile;
        string seq;
        string regexForFasta;
	bool logodds;
	bool applyToAllChains;
	bool debug;
        int valueType;
        map<string,int> valueMap;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

