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

#include <string>
#include <vector>
#include <sstream>
#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "Quench.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"
#include "SystemRotamerLoader.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "SysEnv.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

string programName = "runQuench";
string programDescription = "This program further optimizes a protein structure by running a side chain quenching algorithm";
string programAuthor = "Jason Donald";
string programVersion = "1.0.2";
string programDate = "Sep 18 2009";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

// Define Input Options and Store them.
struct Options {

	// Storage for the vales of each option
	std::string pdb;
	std::string topfile;
	std::string parfile;
	std::string configfile;
	std::string rotlib;
	std::string outfile;
	std::string rotLevel;
	std::vector<std::string> positions;
	int largeRotNum;
	int smallRotNum;
	bool autoFindPositions;

	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages
	string OPerrors; //the errors from the option parser

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> optional;
};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " --help" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
	cout << "Run as: " << endl;
	cout << " % runQuench \n --pdb <PDB> " << endl;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << " --topfile <TOPFILE> " << endl;
	cout << " --parfile <PARFILE> " << endl;
	cout << " --rotlib <ROTLIB> " << endl;
	cout << " --outfile <OUTPDB> " << endl;
	cout << " --positions <CHAIN_POS> " << endl;
	cout << " --rotlevel <SL85.00> " << endl;
	cout << " --smallRotNum <NUM> " << endl;
	cout << " --largeRotNum <NUM> " << endl;
	cout << " --autoFinPositions" << endl;
	cout << endl;
}

Options setupOptions(int theArgc, char * theArgv[]){
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
	 ******************************************/
	// Required arguments
	opt.required.push_back("pdb");

	// Optional arguments
	opt.optional.push_back("topfile");
	opt.optional.push_back("parfile");
	opt.optional.push_back("rotlib");
	opt.optional.push_back("outfile");
	opt.optional.push_back("positions");
	opt.optional.push_back("largeRotNum");
	opt.optional.push_back("smallRotNum");
	opt.optional.push_back("autoFindPositions");
	opt.optional.push_back("rotlevel"); 

	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);

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

	opt.configfile = OP.getString("configfile");
	
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			return opt;
		}
	}
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		opt.errorMessages = "pdb not specified\n";
		opt.errorFlag = true;
		return opt;
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		opt.topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
		opt.warningMessages += "charmmtopfile not specified, using " + opt.topfile + "\n";
		opt.warningFlag = true;
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		opt.parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
		opt.warningMessages += "charmmparfile not specified, using " + opt.parfile + "\n";
		opt.warningFlag = true;
	}

	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		opt.rotlib = SYSENV.getEnv("MSL_ROTLIB");
		opt.warningMessages += "rotlib not specified, using " + opt.rotlib + "\n";
		opt.warningFlag = true;
	}

	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		opt.outfile = "currentConformation.pdb";
		opt.warningMessages += "outfile not specified, using " + opt.outfile + "\n";
		opt.warningFlag = true;
	}

	opt.positions = OP.getStringVectorJoinAll("positions");

	opt.rotLevel = OP.getString("rotlevel");
	if(OP.fail()) {
		opt.rotLevel = "";
	}

	opt.largeRotNum = OP.getInt("largeRotNum");
	if ((opt.rotLevel == "") && OP.fail()){
		opt.warningMessages += "largeRotNum not specified, using 50\n";
		opt.largeRotNum = 50;
		opt.warningFlag = true;
	}

	opt.smallRotNum = OP.getInt("smallRotNum");
	if ((opt.rotLevel == "") && OP.fail()){
		opt.warningMessages += "smallRotNum not specified, using 5\n";
		opt.smallRotNum = 5;
		opt.warningFlag = true;
	}
	if (opt.largeRotNum < 0 || opt.smallRotNum < 0) {
		opt.errorFlag = true;
		opt.errorMessages += "Need a positive rotamer number\n";
		return opt;
	}

	opt.autoFindPositions = OP.getBool("autoFindPositions");
	if (OP.fail()){
		opt.autoFindPositions = false;
	}

	return opt;
}
