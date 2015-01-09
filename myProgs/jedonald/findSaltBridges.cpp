/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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

#include <iostream>

#include "System.h"
#include "DSSPReader.h"
#include "StrideReader.h"
#include "CartesianPoint.h"
#include "PDBTopology.h"
#include "ResidueSelection.h"
#include "release.h"
#include "OptionParser.h"

using namespace std;

using namespace MSL;


string programName = "findSaltBridges";
string programDescription = "This program finds and aligns residues around a residue/substrate in a pdb";
string programAuthor = "Jason Donald";
string programVersion = "1.0.0";
string programDate = "18 April 2014";
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
	// DSSP/Stride and PDB information
	string pdb; // pdb file

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

		usage();
		exit(1);
	}

	System sys;
	sys.readPdb(opt.pdb);

	char resSel[100];
	sprintf(resSel,"resn %s", "ARG");

	ResidueSelection sel(sys);
	vector<Residue *> aspResidues = sel.select("resn ASP");
	vector<Residue *> gluResidues = sel.select("resn GLU");
	vector<Residue *> argResidues = sel.select("resn ARG");
	vector<Residue *> lysResidues = sel.select("resn LYS");
	vector<Residue *> hisResidues = sel.select("resn HIS or resn HID or resn HIE or resn HSD or resn HSE");

	// Define atoms of interest
	vector<string> argNitrogens;
	argNitrogens.push_back("NE");
	argNitrogens.push_back("NH1");
	argNitrogens.push_back("NH2");

	vector<string> lysNitrogens;
	lysNitrogens.push_back("NZ");

	vector<string> hisNitrogens;
	hisNitrogens.push_back("ND1");
	hisNitrogens.push_back("NE2");

	vector<string> aspOxygens;
	aspOxygens.push_back("OD1");
	aspOxygens.push_back("OD2");

	vector<string> gluOxygens;
	gluOxygens.push_back("OE1");
	gluOxygens.push_back("OE2");

	// For each residue
	for (uint r = 0; r < argResidues.size();r++){
		Residue *basicRes = argResidues[r];

		for (uint d = 0; d < aspResidues.size();d++){
			Residue *acidicRes = aspResidues[d];

			bool selectRes = false;
			for (uint i = 0; i < argNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(argNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(argNitrogens[i]);
					for (uint j = 0; j < aspOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(aspOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(aspOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "ASP " << acidicRes->getResidueNumber() << " ARG " << basicRes->getResidueNumber() << endl;
			}
		}
		for (uint e = 0; e < gluResidues.size();e++){
			Residue *acidicRes = gluResidues[e];

			bool selectRes = false;
			for (uint i = 0; i < argNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(argNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(argNitrogens[i]);
					for (uint j = 0; j < gluOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(gluOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(gluOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "GLU " << acidicRes->getResidueNumber() << " ARG " << basicRes->getResidueNumber() << endl;
			}
		}
	}

	for (uint k = 0; k < lysResidues.size();k++){
		Residue *basicRes = lysResidues[k];

		for (uint d = 0; d < aspResidues.size();d++){
			Residue *acidicRes = aspResidues[d];

			bool selectRes = false;
			for (uint i = 0; i < lysNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(lysNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(lysNitrogens[i]);
					for (uint j = 0; j < aspOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(aspOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(aspOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "ASP " << acidicRes->getResidueNumber() << " LYS " << basicRes->getResidueNumber() << endl;
			}
		}
		for (uint e = 0; e < gluResidues.size();e++){
			Residue *acidicRes = gluResidues[e];

			bool selectRes = false;
			for (uint i = 0; i < lysNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(lysNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(lysNitrogens[i]);
					for (uint j = 0; j < gluOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(gluOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(gluOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "GLU " << acidicRes->getResidueNumber() << " LYS " << basicRes->getResidueNumber() << endl;
			}
		}
	}

	for (uint h = 0; h < hisResidues.size();h++){
		Residue *basicRes = hisResidues[h];

		for (uint d = 0; d < aspResidues.size();d++){
			Residue *acidicRes = aspResidues[d];

			bool selectRes = false;
			for (uint i = 0; i < hisNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(hisNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(hisNitrogens[i]);
					for (uint j = 0; j < aspOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(aspOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(aspOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "ASP " << acidicRes->getResidueNumber() << " HIS " << basicRes->getResidueNumber() << endl;
			}
		}
		for (uint e = 0; e < gluResidues.size();e++){
			Residue *acidicRes = gluResidues[e];

			bool selectRes = false;
			for (uint i = 0; i < hisNitrogens.size() and !selectRes; i++) {
				if (basicRes->atomExists(hisNitrogens[i])) {
					Atom basicAtom = basicRes->getAtom(hisNitrogens[i]);
					for (uint j = 0; j < gluOxygens.size() and !selectRes; j++) {
						if (acidicRes->atomExists(gluOxygens[j])) {
							Atom acidicAtom = acidicRes->getAtom(gluOxygens[j]);
							if (basicAtom.distance(acidicAtom) < 4) {
								selectRes = true;
							}
						}
					}
				}
			}
			if (selectRes) {
				cout << "GLU " << acidicRes->getResidueNumber() << " HIS " << basicRes->getResidueNumber() << endl;
			}
		}
	}
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
	OP.setDefaultArguments(opt.defaultArgs);
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
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");

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
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.pdb = OP.getString("pdb");
	if (OP.fail()) {
		opt.errorMessages = "Option name of pdb file \"pdb\" not specified";
		opt.errorFlag = true;
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
	cout << " % findSaltBridges --pdb <pdb file>" << endl;
	cout << endl;
}

