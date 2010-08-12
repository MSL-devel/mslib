/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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

#include "CoiledCoils.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "Transforms.h"
#include "OptionParser.h"
#include "MslOut.h"
#include <cstdlib>
#include <string>

using namespace std;
using namespace MSL;

string programName = "testCoiledCoils";
string programDescription = "This program generates Coiled Coils based on a given set of parameters";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.2";
string programDate = "28 April 2010";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/
	// MOLECULE
	double r0;
	double w0;
	double alpha;
	double r1;
	double w1;
	double phi1;
	//double dZ;
	double rpr;
	double pitch;
	int nRes;
	string fileName;
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

	CoiledCoils cc;
	PDBWriter pout("/tmp/test.pdb");


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



	/******************************************************************************
	 *
	 *                   === GEVORG/CRICKS CC GENERATION ===
	 *
	 ******************************************************************************/

	//cc.degGevorgCoiledCoil(opt.r0, opt.w0, opt.alpha, opt.r1, opt.w1, opt.phi1, opt.dZ, opt.nRes);
	//AtomPointerVector gevorg = cc.getAtomPointers();

	//pout.open(opt.fileName);
	//pout.write(gevorg);
	//pout.close();


	//cc.degNorthCoiledCoil(opt.r0, opt.rpr, opt.pitch, opt.r1, opt.nRes, opt.w1, opt.phi1);
	//AtomPointerVector norths = cc.getAtomPointers();

	//pout.open(opt.fileName);
	//pout.write(norths);
	//pout.close();

	//cc.getCoiledCoilBundle(opt.r0, opt.rpr, opt.pitch, opt.r1, opt.nRes, opt.w1, opt.phi1, 'c', 2);
	//AtomPointerVector norths = cc.getAtomPointers();

	//pout.open(opt.fileName);
	//pout.write(norths);
	//pout.close();

	//norths.subdivideByChainAndPosition(norths);
	//norths.subdivideByChainPositionAndIndentity(norths);

	//System sys;
	//sys.readPdb("/data00/bkmueller/gevorgImplementation/1a93temp.pdb");
	//std::vector<string> paramList;
	//paramList.push_back("A,1");
	//paramList.push_back("B,1");
	//paramList.push_back("C,1");

	//cc.primarySequenceToCoiledCoil(opt.r0, opt.rpr, opt.pitch, opt.r1, opt.nRes, opt.w1, opt.phi1, 'c', 3, sys, paramList);

	//sys.writePdb(opt.fileName);
	cc.getCoiledCoilCricks(4.910, -4.027, -13.530, 2.267, 102.806, 149.984, 0.0, 26);
	AtomPointerVector cricks = cc.getAtomPointers();
	pout.open("/tmp/cricks.pdb");
	pout.write(cricks);
	pout.close();

	cc.getCoiledCoil(4.910, 1.475, 128.206, 2.267, 102.806, 149.984, 0.0, 26);
	AtomPointerVector norths = cc.getAtomPointers();
	pout.open("/tmp/norths.pdb");
	pout.write(norths);
	pout.close();


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

	opt.required.push_back("r0");
	opt.required.push_back("w0");
	opt.required.push_back("alpha");
	opt.required.push_back("r1");
	opt.required.push_back("w1");
	opt.required.push_back("phi1");
	//opt.required.push_back("dZ");
	opt.required.push_back("rpr");
	opt.required.push_back("pitch");
	opt.required.push_back("nRes");
	opt.required.push_back("fileName");
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

	opt.r0 = OP.getDouble("r0");
	if (OP.fail()) {
		opt.errorMessages = "r0 not specified";
		opt.errorFlag = true;
	}
	opt.w0 = OP.getDouble("w0");
	if (OP.fail()) {
		opt.errorMessages = "w0 not specified";
		opt.errorFlag = true;
	}
	opt.alpha = OP.getDouble("alpha");
	if (OP.fail()) {
		opt.errorMessages = "alpha not specified";
		opt.errorFlag = true;
	}
	opt.r1 = OP.getDouble("r1");
	if (OP.fail()) {
		opt.errorMessages = "r1 not specified";
		opt.errorFlag = true;
	}
	opt.w1 = OP.getDouble("w1");
	if (OP.fail()) {
		opt.errorMessages = "w1 not specified";
		opt.errorFlag = true;
	}
	opt.phi1 = OP.getDouble("phi1");
	if (OP.fail()) {
		opt.errorMessages = "phi1 not specified";
		opt.errorFlag = true;
	}
	//opt.dZ = OP.getDouble("dZ");
	//if (OP.fail()) {
	//	opt.errorMessages = "dZ not specified";
	//	opt.errorFlag = true;
	//}
	opt.rpr = OP.getDouble("rpr");
	if (OP.fail()) {
		opt.errorMessages = "rpr not specified";
		opt.errorFlag = true;
	}
	opt.pitch = OP.getDouble("pitch");
	if (OP.fail()) {
		opt.errorMessages = "pitch not specified";
		opt.errorFlag = true;
	}
	opt.nRes = OP.getInt("nRes");
	if (OP.fail()) {
		opt.errorMessages = "dZ not specified";
		opt.errorFlag = true;
	}
	opt.fileName = OP.getString("fileName");
	if (OP.fail()) {
		opt.errorMessages = "fileName not specified";
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
	cout << " % testCoiledCoils --r0 <superhelical radius> --w0 <superhelical frequency> --alpha <pitch angle> --r1 <alpha helical radius> --w1 <alpha helical frequency> --phi1 <alpha helical phase> --rpr <rise per residue> --pitch <superhelical pitch> --nRes <number of residues> --fileName <nameOfOutPutFile.pdb>" << endl;
	cout << endl;
}


