#include <iostream>

#include "System.h"
#include "Atom3DGrid.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "SasaCalculator.h"
#include "release.h"
#include "OptionParser.h"

using namespace std;

using namespace MSL;


string programName = "calculateSasa";
string programDescription = "This programs calculates the solvent exposed surface area (SASA) or a PDB";
string programAuthor = "Sabareesh Subramaniam, Alessandro Senes";
string programVersion = "1.0.1";
string programDate = "7 October 2009";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

 // variables to store time at start and end of application
time_t start_time, finish_time;

/******************************************
 *      =======  STRUCTURE Options  =======
 * 
 *  This structure contains all the user
 *  defined options
 * 
 *  The function
 *     Options parseOptions(int _argc, char * _argv[], Options defaults)
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
	string pdb; // first pdb
	double probeRadius; // normally 1.4 A
	string outputPdb; // the new pdb
	bool writePdb; // only calculate the rmsd between the selections
	string outputFile; // do not print the header output lines
	int sphereDensity; // number of points in the surface sphere
	bool reportByResidue; // if false it reports by atom
	//string outputdir;  // the directory with the output for the run
	string configfile;  // name of the configuration file

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	ofstream * cout_fs; // file stream for the redirected standard output
	ofstream * cerr_fs; // file stream for the redirected standard error
	streambuf * coutbuf; // pointer for storing the cout streambuffer
	streambuf * cerrbuf; // pointer for storing the cerr streambuffer
	streambuf * coutsbuf; // pointer to the filestr streambuffer
	streambuf * cerrsbuf; // pointer to the filestr streambuffer
	bool coutRedirected; // true if cout was redirected

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser



};




/******************************************
 *  
 *  =======  PREDECLARATION OF THE FUNCTIONS  =======
 *
 ******************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);


/******************************************
 *  
 *         =======  MAIN  =======
 *
 ******************************************/
int main(int argc, char* argv[]) {

	// store start time
	time(&start_time);

	cout << "Start: " << ctime(&start_time);
	cout << "Program: " << programName << endl;
	cout << "Program description: " << programDescription << endl;
	cout << "Program author: " << programAuthor << endl;
	cout << "Program version: " << programVersion << " " << programDate << endl;
	cout << "MSL version: " << mslVersion << " " << mslDate << endl;

	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are
	 *  not always required
	 ******************************************************************************/
	Options defaults;
	defaults.probeRadius = 1.4;
	defaults.sphereDensity = 2000;
	defaults.reportByResidue = false;

	/******************************************************************************
	 *                             === OPTION PARSING ===
	 *
	 *  Parse the command line (and possibly input file)
	 *  options. It will also create an input file based
	 *  on the current options in the output directory
	 *  that can be used to re-run the program with the same
	 *  configuration
	 ******************************************************************************/
	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		//printErrors(opt);
		usage();
		exit(1);
	}
	cout << endl;
	cout << "Pdb: " << opt.pdb << endl;
	cout << "Probe radius: " << opt.probeRadius << endl;
	cout << "Sphere density: " << opt.sphereDensity << endl;
	if (opt.reportByResidue) {
		cout << "Reporting by: residue" << endl;
	} else {
		cout << "Reporting by: atom" << endl;
	}
	if (opt.writePdb) {
		cout << "Output pdb: " << opt.outputPdb << endl;
	}
	cout << endl;

	System sys1;
	sys1.readPdb(opt.pdb);

	PDBReader pin;
	if(!pin.open(opt.pdb)) {
		cerr << "Unable to open pdb " << opt.pdb << endl;
		exit(1);
	} 
	pin.read();
	pin.close();

	System sys;
	sys.addAtoms(pin.getAtoms());
	AtomPointerVector atoms = sys.getAtoms();

	SasaCalculator b(atoms, opt.probeRadius, opt.sphereDensity); 
	if (opt.writePdb) {
		b.setTempFactorWithSasa(true);
	}
	b.calcSasa();
	string output;
	if (opt.reportByResidue) {
		output = b.getSasaTable(false);
		if (opt.writePdb) {
			for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end();k++) {
				// set the residue sasa in the b-factor
				(*k)->setTempFactor((*k)->getParentResidue()->getSasa());
			}
		}
	} else {
		output = b.getSasaTable();
	}

	if (opt.outputFile != "") {
		ofstream out_fs;
		out_fs.open(opt.outputFile.c_str());
		if (out_fs.fail()) {
			cerr << "Cannot write output file " << opt.outputFile << endl;
		} else {
			cout << "Output written to " << opt.outputFile << endl;
			out_fs << output << endl;
		}
	} else {
		cout << output << endl;
	}


	if (opt.writePdb) {
		PDBWriter writer;
		if (!writer.open(opt.outputPdb)) {
			cerr << "Cannot open " << opt.outputPdb << " for writing" << endl;
		}
		writer.write(atoms);
		writer.close();
	}

	
}

Options parseOptions(int _argc, char * _argv[], Options defaults) {

	/******************************************<< opt.outputFile << endl;
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
	vector<string> allowed;

	opt.required.push_back("pdb");
	opt.allowed.push_back("probeRadius");
	opt.allowed.push_back("writePdb");
	opt.allowed.push_back("sphereDensity");
	opt.allowed.push_back("outputPdb");
	opt.allowed.push_back("reportByResidue");
	opt.allowed.push_back("outputFile");
	//opt.allowed.push_back("outputdir");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("v"); // -v is equivalent to --version
	opt.allowed.push_back("help"); // --help
	opt.allowed.push_back("h"); // -h is equivalent to --help

	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	opt.defaultArgs.push_back("pdb");


	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
//	if (OP.fail()) {
//		opt.version = OP.getBool("v");
//	}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
//	if (OP.fail()) {
//		opt.help = OP.getBool("h");
//	}

	if (opt.help) {
		help(defaults);
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	//opt.commandName = OP.getCommandName();
	//if (!OP.checkOptions()) {
	//	opt.errorFlag = true;
	//	opt.disallowed = OP.getDisallowedOptions();
	//	opt.missing = OP.getMissingOptions();
	//	opt.ambiguous = OP.getAmbiguousOptions();
	//	return opt;
	//}
	opt.errorFlag = false;
	opt.warningFlag = false;


	opt.pdb = OP.getString("pdb");
	if (OP.fail()) {
		opt.errorMessages = "Option name of first pdb file \"pdb\" not specified";
		opt.errorFlag = true;
	}

	opt.probeRadius = OP.getDouble("probeRadius");
	if (OP.fail()) {
		opt.probeRadius = defaults.probeRadius;
	}

	opt.sphereDensity = OP.getInt("sphereDensity");
	if (OP.fail()) {
		opt.sphereDensity = defaults.sphereDensity;
	}

	opt.writePdb = OP.getBool("writePdb");

	opt.reportByResidue = OP.getBool("reportByResidue");
	if (OP.fail()) {
		opt.reportByResidue = defaults.reportByResidue;
	}


	opt.outputPdb = OP.getString("outputPdb");
	if (OP.fail()) {
		opt.warningMessages = "Option name of output aligned second pdb file \"outputPdb\" not specified";
		opt.warningFlag = true;
		string base = MslTools::pathTail(opt.pdb);
		base = MslTools::pathRoot(base);
		opt.outputPdb = base + (string)"-sasa.pdb";
	}

	opt.outputFile = OP.getString("outputFile");

	/*
	int index = 0;
	while (true) {
		string sele = OP.getString("sele1", index);
		if (OP.fail()) {
			break;
		}
		opt.sele1.push_back(sele);
		index++;
	}

	index = 0;
	while (true) {
		string sele = OP.getString("sele2", index);
		if (OP.fail()) {
			break;
		}
		opt.sele2.push_back(sele);
		index++;
	}

	opt.noAlign = OP.getBool("noAlign");
	*/

	// return the Options structure
	return opt;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help(Options defaults) {
       	cout << "Run  as:" << endl;
	cout << " % calculateSasa --pdb <pdbfile.pdb> [--probeRadius 1.4] [--writePdb] [--sphereDensity 2000] [--reportByResidue] [--outputFile <file.txt>]" << endl;
	cout << endl;
}

