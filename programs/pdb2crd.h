#include <iostream>
#include <cstdlib>

#include "System.h"
#include "FormatConverter.h"
#include "OptionParser.h"
#include "CRDWriter.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  Conversion utility from PDB to CHARMM atom and residue names
 *******************************************************************/

string programName = "pdb2crd";
string programDescription = "This program converts a PDB file into a CHARMM crd file or a PDB file with CHARMM atom/residue names";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.0";
string programDate = "Apr 16 2013";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

/******************************************
 *      =======  STRUCTURE Options  =======
 * 
 *  This structure contains all the user
 *  defined options
 * 
 *  The function
 *     Options parseOptions(int theArgc, char * theArgv[], Options defaults)
 *  parses the argument (using the OptionParser object) and return
 *  a filled Options structure (opt) which is used during the program
 *
 *  Another Options structure (defaults) stores some default
 *  values that will be given to opt if the user did not input them as an
 *  argument
 * 
 ******************************************/
struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/
	// MOLECULE
	string input; // the input pdb or crd file
	string output; // the output pdb or crd file
	bool noNameTranslation; // change file format (crd to pdb or vice versa) but do not change names
	//bool crd2pdb; // convert from CHARMM to PDB (PDB to charmm is the default)
	bool pdb2; // use PDB v.2 names (default is v. 3)
	bool charmm19; // use CHARMM 19 names (default is CHARMM 22)

	/***** MANAGEMENT VARIABLES ******/
	bool version; // ask for the program version
	bool help; // ask for the program help


	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> required; //list of required options
	vector<string> allowed; //list of allowed options
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"

	string OPerrors; //the errors from the option parser

};
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	//cout << "   % " << programName << " --input <in.pdb|in.crd> [--output <out.pdb|out.crd>] [--noNameTranslation] [--crd2pdb] [--pdb2] [--charmm19]" << endl;
	cout << "   % " << programName << " --input <in.pdb|in.crd> [--output <out.pdb|out.crd>] [--noNameTranslation] [--pdb2] [--charmm19]" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
	usage();
}

Options parseOptions(int _argc, char * _argv[]) {

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	opt.required.push_back("input"); // the input pdb or crd file

	opt.allowed.push_back("output"); // the output pdb or crd filee 
	opt.allowed.push_back("noNameTranslation"); //  change file format (crd to pdb or vice versa) but do not change names
	//opt.allowed.push_back("crd2pdb"); // convert from CHARMM to PDB (PDB to charmm is the default) 
	opt.allowed.push_back("pdb2"); // use PDB v.2 names (default is v. 3) 
	opt.allowed.push_back("charmm19"); // use CHARMM 19 names (default is CHARMM 22)

	opt.defaultArgs.push_back("input");

	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("help"); // --help
	opt.allowed.push_back("v"); // -v is equivalent to --version
	opt.allowed.push_back("h"); // -h is equivalent to --help

	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.readArgv(_argc, _argv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version prints the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");

	if (opt.help) {
		help();
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;
	opt.errorMessages = "";
	opt.warningMessages = "";

	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}

	opt.input = OP.getString("input");

	opt.output = OP.getString("output");

	opt.noNameTranslation = OP.getBool("noNameTranslation");
	if(OP.fail()) {
		opt.noNameTranslation = false;
	}

//	opt.crd2pdb = OP.getBool("crd2pdb");
//	if(OP.fail()) {
//		opt.crd2pdb = false;
//	}

	opt.pdb2 = OP.getBool("pdb2");
	if(OP.fail()) {
		opt.pdb2 = false;
	}

	opt.charmm19 = OP.getBool("charmm19");
	if(OP.fail()) {
		opt.charmm19 = false;
	}

	
	// return the Options structure
	return opt;

}
