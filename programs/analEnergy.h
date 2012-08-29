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
#include "EnergeticAnalysis.h"
#include "CharmmSystemBuilder.h"
#include "PolymerSequence.h"
#include "OptionParser.h"
#include "ResidueSelection.h"
#include "release.h"
#include "CharmmEnergyCalculator.h"
#include "AtomSelection.h"
#include "SysEnv.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

string programName = "analEnergy";
string programDescription = "This program analyzes the energy of the total protein, or if selections are given, the pairwise energy of residues";
string programAuthor = "Dan Kulp";
string programVersion = "1.0.2";
string programDate = "Sep 18 2009";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

struct Options {

	// Storage for the vales of each option
	std::string pdb;
	std::string topfile;
	std::string parfile;
	std::string configfile;
        std::string selection;
        std::string selection2;
	bool pymolOutput;
	std::vector<std::string> positions;
	double cuton;
	double cutoff;
	double cutnb;

	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

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
        cout << " % analEnergy \n --pdb <PDB> " << endl;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << " --positions <CHAIN_RESNUM>" << endl;
	cout << " --select <SEL>" << endl;
	cout << " --select2 <SEL>" << endl;
	cout << " --topfile <TOPFILE>" << endl;
	cout << " --parfile <PARFILE>" << endl;
	cout << " --cuton <nbcuton>" << endl;
	cout << " --cutoff <nbcutoff>" << endl;
	cout << " --cutnb <nbcutnb>" << endl;
	cout << " --pymolOutput" << endl;
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

	// Options
	opt.optional.push_back("positions");
	opt.optional.push_back("select");
	opt.optional.push_back("select2");
	opt.optional.push_back("parfile");
	opt.optional.push_back("topfile");
	opt.optional.push_back("pymolOutput");
	opt.optional.push_back("cuton");
	opt.optional.push_back("cutoff");
	opt.optional.push_back("cutnb");
	opt.optional.push_back("help");
	opt.optional.push_back("version"); // --version

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

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		opt.errorMessages = "pdb not specified\n";
		opt.errorFlag = true;
		return opt;
	}

	opt.positions = OP.getStringVector("positions");
	if (OP.fail()){
		opt.warningMessages += "no positions chosen to analyze\n";
		opt.warningFlag = true;
	}	
	
	opt.selection = OP.getString("select");
	if (OP.fail()){
		opt.warningMessages += "no selections to analyze\n";
		opt.warningFlag = true;
	}
	opt.selection2 = OP.getString("select2");
	
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
	
	opt.pymolOutput = OP.getBool("pymolOutput");
	if (OP.fail()){
		opt.pymolOutput = false;
	}
	
	opt.cuton = OP.getDouble("cuton");
	if(OP.fail()) {
		opt.warningMessages += "cuton not specified, using 9.0\n";
		opt.warningFlag = true;
		opt.cuton = 9.0;
	}
	opt.cutoff = OP.getDouble("cutoff");
	if(OP.fail()) {
		opt.warningMessages += "cutoff not specified, using 10.0\n";
		opt.warningFlag = true;
		opt.cutoff = 10.0;
	}
	opt.cutnb = OP.getDouble("cutnb");
	if(OP.fail()) {
		opt.warningMessages += "cutnb not specified, using 11.0\n";
		opt.warningFlag = true;
		opt.cutnb = 11.0;
	}

	return opt;
}

