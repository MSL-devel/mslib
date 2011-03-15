/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

#include <iostream>

#include "AtomContainer.h"
#include "AtomSelection.h"
#include "OptionParser.h"

using namespace std;

using namespace MSL;


string programName = "calculateDistanceOrAngle";
string programDescription = "This programs calculates a distance (2 atoms specified), an angle (3 atoms) or a dihedral (4 atoms)";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.0";
string programDate = "15 March 2011";
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
	vector<string> atoms; // the atoms by id (i.e. A,7,CA)
	vector<int> atomIndeces; // the atoms by index (i.e 3)

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
	vector<vector<string> > oneRequired; //list of required options
	vector<vector<string> > mutuallyExclusive; //list of required options

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
void help();


/******************************************
 *  
 *         =======  MAIN  =======
 *
 ******************************************/
int main(int argc, char* argv[]) {

	// store start time
	time(&start_time);

	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are
	 *  not always required
	 ******************************************************************************/
	Options defaults;

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

	AtomContainer mol;
	if (!mol.readPdb(opt.pdb)) {
		cerr << "Cannot read PDB file " << opt.pdb << endl;
		exit(1);
	}

	vector<Atom*> pAtoms;
	if (opt.atoms.size() > 0) {
		if (opt.atoms.size() < 2) {
			cerr << "Not enough atoms specified" << endl;
			exit(1);
		}

		if (mol.atomExists(opt.atoms[0])) {
			pAtoms.push_back(&mol.getLastFoundAtom());
		} else {
			cerr << "Atom 1 \"" << opt.atoms[0] << "\" not found" << endl;
			exit(1);
		}
		if (mol.atomExists(opt.atoms[1])) {
			pAtoms.push_back(&mol.getLastFoundAtom());
		} else {
			cerr << "Atom 2 \"" << opt.atoms[1] << "\" not found" << endl;
			exit(1);
		}
		if (opt.atoms.size() >= 3) {
			// an angle
			if (mol.atomExists(opt.atoms[2])) {
				pAtoms.push_back(&mol.getLastFoundAtom());
			} else {
				cerr << "Atom 3 \"" << opt.atoms[2] << "\" not found" << endl;
				exit(1);
			}
			if (opt.atoms.size() == 3) {
				cout << pAtoms[0]->angle(*pAtoms[1], *pAtoms[2]) << endl;
				exit(0);
			} else {
				// a dihedral
				if (mol.atomExists(opt.atoms[3])) {
					pAtoms.push_back(&mol.getLastFoundAtom());
				} else {
					cerr << "Atom 4 \"" << opt.atoms[3] << "\" not found" << endl;
					exit(1);
				}
				cout << pAtoms[0]->dihedral(*pAtoms[1], *pAtoms[2], *pAtoms[3]) << endl;;
				exit(0);
			}
		}
		// a distance
		cout << pAtoms[0]->distance(*pAtoms[1]) << endl;
	} else {
		if (opt.atomIndeces.size() < 2) {
			cerr << "Not enough atoms specified" << endl;
			exit(1);
		}

		if (mol.atomSize() >= opt.atomIndeces[0] && opt.atomIndeces[0] > 0) {
			pAtoms.push_back(&mol[opt.atomIndeces[0]-1]);
		} else {
			cerr << "Atom 1 \"" << opt.atomIndeces[0] << "\" out of index (" << mol.atomSize() << ")" << endl;
			exit(1);
		}
		if (mol.atomSize() >= opt.atomIndeces[1] && opt.atomIndeces[1] > 0) {
			pAtoms.push_back(&mol[opt.atomIndeces[1]-1]);
		} else {
			cerr << "Atom 2 \"" << opt.atomIndeces[1] << "\" not found" << endl;
			exit(1);
		}
		if (opt.atomIndeces.size() >= 3) {
			// an angle
			if (mol.atomSize() >= opt.atomIndeces[2] && opt.atomIndeces[2] > 0) {
				pAtoms.push_back(&mol[opt.atomIndeces[2]-1]);
			} else {
				cerr << "Atom 3 \"" << opt.atomIndeces[2] << "\" not found" << endl;
				exit(1);
			}
			if (opt.atomIndeces.size() == 3) {
				cout << pAtoms[0]->angle(*pAtoms[1], *pAtoms[2]) << endl;
				exit(0);
			} else {
				// a dihedral
				if (mol.atomSize() >= opt.atomIndeces[3] && opt.atomIndeces[3] > 0) {
					pAtoms.push_back(&mol[opt.atomIndeces[3]-1]);
				} else {
					cerr << "Atom 4 \"" << opt.atomIndeces[3] << "\" not found" << endl;
					exit(1);
				}
				cout << pAtoms[0]->dihedral(*pAtoms[1], *pAtoms[2], *pAtoms[3]) << endl;;
				exit(0);
			}
		}
		// a distance
		cout << pAtoms[0]->distance(*pAtoms[1]) << endl;
	}

	
	return 0;

	
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
	vector<string> oneRequired;

	opt.required.push_back("pdb");
	opt.allowed.push_back("atoms");
	opt.allowed.push_back("atomIndeces");
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
	opt.oneRequired.push_back(vector<string>());
	opt.oneRequired.back().push_back("atoms");
	opt.oneRequired.back().push_back("atomIndeces");

	opt.defaultArgs.push_back("pdb");


	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.setOneRequired(opt.oneRequired);
	OP.setMutualExclusive(opt.oneRequired);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
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
		help();
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
	//if (OP.fail()) {
	//	opt.errorMessages = "Option name of first pdb file \"pdb\" not specified";
	//	opt.errorFlag = true;
	//}

	opt.atoms = OP.getStringVectorJoinAll("atoms");
	//if (OP.fail()) {
	//	opt.errorMessages = "Option name of first pdb file \"pdb\" not specified";
	//	opt.errorFlag = true;
	//}

	opt.atomIndeces = OP.getIntVectorJoinAll("atomIndeces");
	//if (OP.fail()) {
	//	opt.errorMessages = "Option name of first pdb file \"pdb\" not specified";
	//	opt.errorFlag = true;
	//}

	/*
	int index = 0;
	while (true) {
		string atom = OP.getString("atoms", index);
		if (OP.fail()) {
			break;
		}
		opt.atoms.push_back(atom);
		index++;
	}
	*/

	// return the Options structure
	return opt;

}

void usage() {
	help();
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
	cout << "Program: " << programName << endl;
	cout << "Program description: " << programDescription << endl;
	cout << "Program author: " << programAuthor << endl;
	cout << "Program version: " << programVersion << " " << programDate << endl;
	cout << "MSL version: " << mslVersion << " " << mslDate << endl;
	cout << endl;

   	cout << "You can run either with atom identifies (chain,resnum,atomName comma separated without spaces, i.e. \"N,37,CA\"):" << endl;
	cout << endl;
	cout << " % calculateDistanceOrAngle --pdb <pdbfile.pdb> --atoms <atomId> <atomId> [<atomId> [<atomId>]]" << endl;
	cout << endl;
	cout << "... or with the atom numbers in the PDB" << endl;
	cout << endl;
	cout << " % calculateDistanceOrAngle --pdb <pdbfile.pdb> --atoms <N> <N> [<N> [<N>]]" << endl;
	cout << endl;
	cout << "2 atoms: returns a distance" << endl;
	cout << "3 atoms: returns a angle" << endl;
	cout << "4 atoms: returns a dihedral" << endl;
	cout << endl;
}

