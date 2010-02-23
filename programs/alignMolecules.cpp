#include "PDBReader.h"
#include "PDBWriter.h"
#include "AtomSelection.h"
#include "System.h"
#include "OptionParser.h"
#include "release.h"
#include "AtomSelection.h"
#include "Transforms.h"
#include "MslTools.h"

using namespace std;

using namespace MSL;


string programName = "alignMolecules";
string programDescription = "This programs aligns two PDB based on a subset of atoms";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.2";
string programDate = "23 September 2009";
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
	string pdb1; // first pdb
	string pdb2; // second pdb
	vector<string> sele1; // atom selection one
	vector<string> sele2; // atom selection two
	string outputPdb2; // the new pdb
	bool noAlign; // only calculate the rmsd between the selections
	string outputdir;  // the directory with the output for the run
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
int main(int argc, char *argv[]) {
	
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
	//defaults.charmmParam = "22";
	

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


	cout << "Read pdb 1: " << opt.pdb1 << endl;
	PDBReader pin;
	if(!pin.open(opt.pdb1)) {
		cerr << "Unable to open pdb " << opt.pdb1 << endl;
		exit(1);
	} 
	pin.read();
	pin.close();
	System sys1;
	sys1.addAtoms(pin.getAtoms());

	cout << "Read pdb 2: " << opt.pdb2 << endl;
	if(!pin.open(opt.pdb2)) {
		cerr << "Unable to open pdb " << opt.pdb2 << endl;
		exit(1);
	} 
	pin.read();
	pin.close();
	System sys2;
	sys2.addAtoms(pin.getAtoms());

	AtomPointerVector av1 = sys1.getAtoms();
	AtomPointerVector av2 = sys2.getAtoms();
	
	AtomPointerVector alignAtoms1;
	if (opt.sele1.size() == 0) {
		// select all atoms
		alignAtoms1.insert(alignAtoms1.end(), av1.begin(), av1.end());
	} else {
		AtomSelection sel1(av1);
		for (unsigned int i=0; i<opt.sele1.size(); i++) {
			char c [1000];
			sprintf(c, "keyatoms, %s", opt.sele1[i].c_str());
			AtomPointerVector selAtom = sel1.select(c);
			alignAtoms1.insert(alignAtoms1.end(), selAtom.begin(), selAtom.end());
		}
	}
	cout << "Selected " << alignAtoms1.size() << " reference atoms for pdb " << opt.pdb1 << endl;

	AtomPointerVector alignAtoms2;
	if (opt.sele2.size() == 0) {
		// select all atoms
		alignAtoms2.insert(alignAtoms2.end(), av2.begin(), av2.end());
	} else {
		AtomSelection sel2(av2);
		for (unsigned int i=0; i<opt.sele2.size(); i++) {
			char c [1000];
			sprintf(c, "keyatoms, %s", opt.sele2[i].c_str());
			AtomPointerVector selAtom = sel2.select(c);
			alignAtoms2.insert(alignAtoms2.end(), selAtom.begin(), selAtom.end());
		}
	}
	cout << "Selected " << alignAtoms2.size() << " reference atoms for pdb " << opt.pdb2 << endl;

	if (alignAtoms1.size() != alignAtoms2.size()) {
		cerr << "The number of atoms selected for pdb 1 (" << alignAtoms1.size() << ") does not match the number of atoms selected for pdb 2 (" << alignAtoms2.size() << ")" << endl;
		exit(1);
	}

	cout << endl;
	cout << "Set 1 ================================" << endl;
	cout << endl;
	cout << alignAtoms1;
	cout << endl;
	cout << "Set 2 ================================" << endl;
	cout << endl;
	cout << alignAtoms2;
	cout << endl;
	cout << "      ================================" << endl;

	if (opt.noAlign) {
		// only calc the rmsd
		double rmsd = alignAtoms2.rmsd(alignAtoms1);
		cout << endl;
		cout << opt.pdb2 << " was NOT aligned to " << opt.pdb1 << "." << endl;
		cout << "RMSD " << rmsd << endl;
	} else {
		Transforms tm;
		if (!tm.align(alignAtoms2,alignAtoms1,av2)) {
			cerr << "Alignment failed!" << endl;
			exit(1);
		}
		//double rmsd = tm.getRMSD();
		double rmsd = alignAtoms2.rmsd(alignAtoms1);
		Matrix rotMatrix = tm.getLastRotationMatrix();
		CartesianPoint translation = tm.getLastTranslation();
		cout << rotMatrix << endl;
		cout << translation << endl;
		PDBWriter writer;
		if (!writer.open(opt.outputPdb2)) {
			cerr << "Cannot open " << opt.outputPdb2 << " for writing" << endl;
		}
		writer.write(av2);
		writer.close();
		cout << endl;
		cout << "Aligned " << opt.pdb2 << " to " << opt.pdb1 << "." << endl;
		cout << "RMSD " << rmsd << endl;
		cout << "Written to " << opt.outputPdb2 << endl;
	}



	return 0;

}

Options parseOptions(int _argc, char * _argv[], Options defaults) {

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
	vector<string> allowed;

	opt.required.push_back("pdb1");
	opt.required.push_back("pdb2");
	opt.allowed.push_back("sele1");
	opt.allowed.push_back("sele2");
	opt.allowed.push_back("outputPdb2");
	opt.allowed.push_back("noAlign");
	opt.allowed.push_back("outputdir");
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

	opt.defaultArgs.push_back("configfile");


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
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

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


	opt.pdb1 = OP.getString("pdb1");
	if (OP.fail()) {
		opt.errorMessages = "Option name of first pdb file \"pdb1\" not specified";
		opt.errorFlag = true;
	}

	opt.pdb2 = OP.getString("pdb2");
	if (OP.fail()) {
		opt.errorMessages = "Option name of second pdb file \"pdb2\" not specified";
		opt.errorFlag = true;
	}

	opt.outputPdb2 = OP.getString("outputPdb2");
	if (OP.fail()) {
		opt.warningMessages = "Option name of output aligned second pdb file \"outputPdb2\" not specified";
		opt.warningFlag = true;
		string base = MslTools::pathTail(opt.pdb2);
		base = MslTools::pathRoot(base);
		opt.outputPdb2 = base + (string)"-aligned.pdb";
	}

	int index = 0;
	while (true) {
		string sele = OP.getString("sele1", index);
		if (OP.fail()) {
			break;
		}
		opt.sele1.push_back(sele);
		index++;
	}
	//if (opt.sele1.size() == 0) {
	//	opt.sele1.push_back("all");
	//}

	index = 0;
	while (true) {
		string sele = OP.getString("sele2", index);
		if (OP.fail()) {
			break;
		}
		opt.sele2.push_back(sele);
		index++;
	}
	//if (opt.sele2.size() == 0) {
	//	opt.sele2.push_back("all");
	//}

	opt.noAlign = OP.getBool("noAlign");

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
	cout << " % alignMolecules --pdb1 <pdbfile.pdb> --pdb2 <pdbToBeAligned.pdb> [--sele1 <atom selection> --sele2 <atom selection> [--sele1 <atom selection> --sele2 <atom selection>]] [--outputPdb2 <output.pdb>] [--noAlign]" << endl;
	cout << endl;
	cout << " **************************************************************************************" << endl;
	cout << " * NOTE: No order is assumed WITHIN a selection but if multiple --sele1/sele2 are     *" << endl;
	cout << " *       given order is preserved.                                                    *" << endl;
	cout << " *       If no selection is given, all atoms are used.                                *" << endl;
	cout << " *       if option --noAlign is given, the rmsd is calculated but no pdb is produced. *" << endl;
	cout << " **************************************************************************************" << endl;
	cout << endl;
}

