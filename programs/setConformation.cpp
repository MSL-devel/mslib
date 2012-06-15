/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include "AtomBondBuilder.h"
#include "System.h"
#include "Transforms.h"
#include "OptionParser.h"

/*******************************************************************************
 *  This program takes a PDB and edits the conformation of the molecule by changing
 *  one or more of its bond distances, angles or dihedrals
 *******************************************************************************/

using namespace std;
using namespace MSL;

string programName = "setConformation";
string programDescription = "This programs allows to edit degrees of freedom in a protein";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.1";
string programDate = "16 April 2010";
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
struct DoF {
	vector<string> atomNames;
	double value;
	string type; // bond, angle, dihedral, improper 
};
struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/
	// MOLECULE
	string pdb; // first pdb
	vector<DoF> edits; // atom selection one
	string outputPdb; // the new pdb
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


int main(int argc,char *argv[]) {

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
	Options defaults; // currently none
	

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
	System sys;
	if (!sys.readPdb(opt.pdb)) {
		cerr << "ERROR: Cannot read pdb file " << opt.pdb << endl;
		exit(1);
	}
	cout << endl;
	cout << "Read pdb file " << opt.pdb << endl;

	AtomBondBuilder abb;
	abb.buildConnections(sys.getAtomPointers());

	Transforms T;
	cout << endl;
	for (unsigned int i=0; i<opt.edits.size(); i++) {
		if (opt.edits[i].type == "bond") {
			if (opt.edits[i].atomNames.size() != 2) {
				cerr << "ERROR: Type bond requires 2 atom name" << endl;
				exit(1);
			}
			Atom * pAtom1 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[0])) {
				pAtom1 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[0] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom2 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[1])) {
				pAtom2 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[1] << " not found" << endl;
				exit(1);
			}
			T.setBondDistance(*pAtom1, *pAtom2, opt.edits[i].value);
			cout << i+1 << " Set bond [" << opt.edits[i].atomNames[0] << ", " << opt.edits[i].atomNames[1] << "] to " << opt.edits[i].value << " (" << pAtom1->distance(*pAtom2) << ")" << endl;
		} else if (opt.edits[i].type == "angle") {
			if (opt.edits[i].atomNames.size() != 3) {
				cerr << "ERROR: Type angle requires 3 atom name" << endl;
				exit(1);
			}
			Atom * pAtom1 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[0])) {
				pAtom1 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[0] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom2 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[1])) {
				pAtom2 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[1] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom3 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[2])) {
				pAtom3 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[2] << " not found" << endl;
				exit(1);
			}
			T.setBondAngle(*pAtom1, *pAtom2, *pAtom3, opt.edits[i].value);
			cout << i+1 << " Set angle [" << opt.edits[i].atomNames[0] << ", " << opt.edits[i].atomNames[1] << ", " << opt.edits[i].atomNames[2] << "] to " << opt.edits[i].value << " (" << pAtom1->angle(*pAtom2, *pAtom3) << ")" << endl;
		} else if (opt.edits[i].type == "dihedral" || opt.edits[i].type == "improper") {
			if (opt.edits[i].atomNames.size() != 4) {
				cerr << "ERROR: Type dihedral/improper requires 4 atom name" << endl;
				exit(1);
			}
			Atom * pAtom1 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[0])) {
				pAtom1 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[0] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom2 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[1])) {
				pAtom2 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[1] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom3 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[2])) {
				pAtom3 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[2] << " not found" << endl;
				exit(1);
			}
			Atom * pAtom4 = NULL;
			if (sys.atomExists(opt.edits[i].atomNames[3])) {
				pAtom4 = &(sys.getLastFoundAtom());
			} else {
				cout << "ERROR: Atom " << opt.edits[i].atomNames[3] << " not found" << endl;
				exit(1);
			}
			if (opt.edits[i].type == "dihedral") {
				T.setDihedral(*pAtom1, *pAtom2, *pAtom3, *pAtom4, opt.edits[i].value);
				cout << i+1 << " Set dihedral [" << opt.edits[i].atomNames[0] << ", " << opt.edits[i].atomNames[1] << ", " << opt.edits[i].atomNames[2] << ", " << opt.edits[i].atomNames[3] << "] to " << opt.edits[i].value << " (" << pAtom1->dihedral(*pAtom2, *pAtom3, *pAtom4) << ")" << endl;
			} else {
				T.setImproper(*pAtom1, *pAtom2, *pAtom3, *pAtom4, opt.edits[i].value);
				cout << i+1 << " Set improper [" << opt.edits[i].atomNames[0] << ", " << opt.edits[i].atomNames[1] << ", " << opt.edits[i].atomNames[2] << ", " << opt.edits[i].atomNames[3] << "] to " << opt.edits[i].value << opt.edits[i].value << " (" << pAtom1->dihedral(*pAtom2, *pAtom3, *pAtom4) << ")" << endl;
			}
		} else {
			cout << "ERROR: Degree of freedom type \"" << opt.edits[i].type << "\" not recognized" << endl;
			exit(1);
		}
	}
	cout << endl;

	if (!sys.writePdb(opt.outputPdb)) {
		cout << "ERROR: cannot write pdb file " << opt.outputPdb << endl;
		exit(1);
	}
	cout << "Output file written to " << opt.outputPdb << endl;

	/*
	AtomPointerVector av = sys.getAtomPointers();

	
	cout << endl;
	cout << "===========================================" << endl;
	cout << endl;
	Atom * pAtomN4 = NULL;
	if (sys.atomExists("A, 4, N")) {
		pAtomN4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomN4 << endl;
	} else {
		cout << "Atom A 4 N not found" << endl;
		exit(1);
	}
	Atom * pAtomCA4 = NULL;
	if (sys.atomExists("A, 4, CA")) {
		pAtomCA4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCA4 << endl;
	} else {
		cout << "Atom A 4 CA not found" << endl;
		exit(1);
	}
	Atom * pAtomC4 = NULL;
	if (sys.atomExists("A, 4, C")) {
		pAtomC4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomC4 << endl;
	} else {
		cout << "Atom A 4 C not found" << endl;
		exit(1);
	}

	// PHE 5
	Atom * pAtomN5 = NULL;
	if (sys.atomExists("A, 5, N")) {
		pAtomN5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomN5 << endl;
	} else {
		cout << "Atom A 5 N not found" << endl;
		exit(1);
	}
	Atom * pAtomCA5 = NULL;
	if (sys.atomExists("A, 5, CA")) {
		pAtomCA5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCA5 << endl;
	} else {
		cout << "Atom A 5 CA not found" << endl;
		exit(1);
	}
	Atom * pAtomCB5 = NULL;
	if (sys.atomExists("A, 5, CB")) {
		pAtomCB5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCB5 << endl;
	} else {
		cout << "Atom A 5 CB not found" << endl;
		exit(1);
	}
	Atom * pAtomCG5 = NULL;
	if (sys.atomExists("A, 5, CG")) {
		pAtomCG5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCG5 << endl;
	} else {
		cout << "Atom A 5 CG not found" << endl;
		exit(1);
	}

	// save the original coordinates
	cout << "Save the coordinates of the original pdb " << inputPdb << endl;
	av.saveCoor("orig");

	Transforms T;

	cout << " ==============================" << endl;
	cout << "Set the PHE 5 (CA, CB) bond to 1.65" << endl;
	cout << "  Bond before " << pAtomCA5->distance(*pAtomCB5) << endl;
	T.setBondDistance(*pAtomCA5, *pAtomCB5, 1.65);
	cout << "  Bond after " << pAtomCA5->distance(*pAtomCB5) << endl;
	string outPdb = "/tmp/testEdit-bond.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 (N, CA, CB) angle to 112.47" << endl;
	cout << "  Angle before " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	T.setBondAngle(*pAtomN5, *pAtomCA5, *pAtomCB5, 112.47);
	cout << "  Angle after " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	outPdb = "/tmp/testEdit-angle.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 chi1 (N, CA, CB, CG) to 178.34" << endl;
	cout << "  Chi1 before " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	T.setDihedral(*pAtomN5, *pAtomCA5, *pAtomCB5, *pAtomCG5, 178.34);
	cout << "  Chi1 after " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	outPdb = "/tmp/testEdit-chi1.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the THR 4 psi (N, CA, C, N+1) to 175.67" << endl;
	cout << "  Psi before " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	T.setDihedral(*pAtomN4, *pAtomCA4, *pAtomC4, *pAtomN5, 175.67);
	cout << "  Psi after " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	outPdb = "/tmp/testEdit-psi.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);
	*/

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

	opt.required.push_back("pdb");
	opt.required.push_back("type");
	opt.required.push_back("atomNames");
	opt.required.push_back("value");
	opt.allowed.push_back("outputPdb");
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


	opt.pdb = OP.getString("pdb");
	if (OP.fail()) {
		opt.errorMessages = "Option name of first pdb file \"pdb\" not specified";
		opt.errorFlag = true;
	}

	opt.outputPdb = OP.getString("outputPdb");
	if (OP.fail()) {
		opt.warningMessages = "Option name of output edited pdb file \"outputPdb\" not specified";
		opt.warningFlag = true;
		string base = MslTools::pathTail(opt.pdb);
		base = MslTools::pathRoot(base);
		opt.outputPdb = base + (string)"-edited.pdb";
	}

	int index = 0;
	while (true) {
		DoF tmp;
		tmp.type = OP.getString("type", index);
		if (OP.fail()) {
			break;
		}
		tmp.atomNames = OP.getStringVector("atomNames", index);
		if (OP.fail()) {
			break;
		}
		tmp.value = OP.getDouble("value", index);
		if (OP.fail()) {
			break;
		}
		opt.edits.push_back(tmp);
		index++;
	}

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
	cout << " % setConformation --pdb <pdbfile> --type <bond|angle|dihedral|improper> --atomNames <atomId atomId ...> [--type <bond|angle|dihedral|improper> --atomNames <atomId atomId ...>] [--outputPdb <output.pdb>]" << endl;
	cout << endl;
	cout << " **************************************************************************************" << endl;
	cout << " * NOTE: The atom IDs needs to be chain,number,name comma separated, no spaces and    *" << endl;
	cout << " *       spaces between the atoms. Example:                                           *" << endl;
	cout << " *       --type bond --atomNames A,23,CA A,23,CB --value 1.52                         *" << endl;
	cout << " *       --type angle --atomNames A,23,CA A,23,CB A,23,CG --value 120.0               *" << endl;
	cout << " *       --type dihedral --atomNames A,23,CA A,23,CB A,23,CG A,23,CD1 --value=-175.0  *" << endl;
	cout << " *       --type improper --atomNames A,23,N A,23,C A,23,CA A,23,CB --value=-120.0     *" << endl;
	cout << " *                                                                                    *" << endl;
	cout << " *       The last atom is the mobile atom.  All connected atoms will rotate           *" << endl;
	cout << " *       accordingly.                                                                 *" << endl;
	cout << " *                                                                                    *" << endl;
	cout << " *       NOTE: the = sign is necessary for negative values                            *" << endl;
	cout << " *       --value 125 OK   --value=125 OK   --value=-125 OK   --value -125 NOT OK      *" << endl;
	cout << " **************************************************************************************" << endl;
	cout << endl;
}

