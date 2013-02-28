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
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "HydrogenBondBuilder.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SidechainOptimizationManager.h"
#include "SasaCalculator.h"
#include "MslTools.h"
#include "SysEnv.h"
#include "release.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

string programName = "repackSideChains";
string programDescription = "This program repacks all positions in a given protein using a given rotamer library and prints out the repacked structure in PDB format";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.2";
string programDate = "Feb 11 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct Options {

        Options(){
		annealShapeMap["CONSTANT"]    = MonteCarloManager::CONSTANT;
		annealShapeMap["LINEAR"]      = MonteCarloManager::LINEAR;
		annealShapeMap["SIGMOIDAL"]   = MonteCarloManager::SIGMOIDAL;
		annealShapeMap["EXPONENTIAL"] = MonteCarloManager::EXPONENTIAL;
		annealShapeMap["SOFT"]        = MonteCarloManager::SOFT;
		shape = annealShapeMap["LINEAR"];
	}

	string commandName; // name of this program
	string pdbFile; // file containing the atom coordinates in PDB format
	string rotlibFile; // file containing the atom coordinates in PDB format
	string charmmTopFile;
	string charmmParFile;
	string hbondParFile;
	string solvFile;
	string solvent;

	double cuton;
	double cutoff;
	double cutnb;
	double cuthb;

	string outputPDBFile; // outputPDBFile
	string logfile; // logFile
	string configfile; // configFile

	// Algorithms for repack
	bool runGoldsteinSingles;
	bool runGoldsteinPairs;
	bool runSCMF;
	bool runSCMFBiasedMC;
	bool runUnbiasedMC;
	
	bool runGreedy; // run the greedyOptimizer
	int greedyCycles;

	bool onTheFly;

	unsigned int seed;
	bool useTimeToSeed; // true if no seed is specified

	bool includeCR; // include the crystal rotamer
	bool verbose;

	// MC Parameters
	double startT; 
	double endT; 
	int nCycles; 
        //int shape; 
        MonteCarloManager::ANNEALTYPES shape;
        map<string,MonteCarloManager::ANNEALTYPES> annealShapeMap;
	int maxReject;
	int deltaSteps; 
	double minDeltaE;

	map<string,int> numRots; // number of rotamers for each residue type 
	string rotLevel;
	vector<string> excludeTerms;

	vector<string> fixedPos;
	map<string,bool> fixedPosMap;

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

	string OPerrors; //the errors from the option parser
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
       	cout << "Run  as:" << endl;
	cout << " % repackSideChains \n --pdbfile <pdbfile> " << endl;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << " --rotlibfile <rotlibfile> \n --charmmtopfile <charmmTopFile> \n --charmmparfile <charmmParFile> \n --hbondparfile <hbondParFile> \n --solvfile <solvationFile> --solvent <string> \n --outputpdbfile <outputpdbfile> \n --logfile <logfile> \n --verbose <true/false> \n --cuton <nbcuton> \n --cutoff <nbcutoff> \n --cutnb <nbcutnb> \n --includecrystalrotamer <true/false> (include crystal rotamer)" << endl;
	cout << " --configfile <configfile> \n --rungoldsteinsingles <true/false> \n --rungoldsteinpairs <true/false> \n --runscmf <true/false> \n --runscmfbiasedmc <true/false> \n --rununbiasedmc <true/false> --rungreedy <true/false> --greedyCycles <int>" << endl;
	cout << "--excludeenergyterm <term1> --excludeenergyterm <term2> \n   [Terms can be CHARMM_ANGL,CHARMM_BOND,CHARMM_DIHE,CHARMM_ELEC,CHARMM_IMPR,CHARMM_U-BR,CHARMM_VDW,SCWRL4_HBOND] All terms are implemented by default " << endl;
	cout << endl;
	cout << "Optional MC Parameters " << endl;
	cout << " --mcstarttemp <startT> \n --mcendtemp <endT> \n --mccycles <numCycles> \n --mcshape <CONSTANT/LINEAR/EXPONENTIAL/SIGMOIDAL/SOFT> \n --mcmaxreject <numReject> \n --mcdeltasteps <numsteps> \n --mcmindeltaenergy <minEnergy>" << endl;
	cout << endl;
	cout << "Optional - Num of rotamers per residue type (-ve or 0 means the wt rotamer alone will be used)" << endl;
	cout << " --rotlevel <SL85.00>\n --ALA <nALA>\n --ARG <nARG> \n .........\n --HSD <nHSD>\n --HSE <nHSD>\n --HSP <nHSP> \n ....... \n --VAL <nVAL> " << endl;
	cout << endl;
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
	opt.required.push_back("pdbfile"); // PDB

	opt.allowed.push_back("rotlibfile"); // rotamerLibrary 
	opt.allowed.push_back("charmmtopfile"); // rotamerLibrary 
	opt.allowed.push_back("charmmparfile"); // rotamerLibrary 
	opt.allowed.push_back("hbondparfile"); // rotamerLibrary 

	opt.allowed.push_back("outputpdbfile"); // repacked structure will be written to this file
	opt.allowed.push_back("logfile"); // all output will be redirected to this logFile
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("excludeenergyterm");
	opt.allowed.push_back("solvfile");
	opt.allowed.push_back("solvent");
	opt.allowed.push_back("excludeenergyterm");
	opt.allowed.push_back("verbose");
	opt.allowed.push_back("onthefly"); // 
	opt.allowed.push_back("seed"); // 
	opt.allowed.push_back("cuton");
	opt.allowed.push_back("cutoff");
	opt.allowed.push_back("cutnb");
	opt.allowed.push_back("cuthb");

	opt.allowed.push_back("rungoldsteinsingles"); // 
	opt.allowed.push_back("rungoldsteinpairs"); // 
	opt.allowed.push_back("runscmf"); // 
	opt.allowed.push_back("runscmfbiasedmc"); // 
	opt.allowed.push_back("rununbiasedmc"); // 
	opt.allowed.push_back("rungreedy"); // 
	opt.allowed.push_back("greedycycles"); // 
	opt.allowed.push_back("includecrystalrotamer");

	opt.allowed.push_back("mcstarttemp"); // 
	opt.allowed.push_back("mcendtemp"); // 
	opt.allowed.push_back("mccycles"); // 
	opt.allowed.push_back("mcshape"); // 
	opt.allowed.push_back("mcmaxreject"); // 
	opt.allowed.push_back("mcdeltasteps"); // 
	opt.allowed.push_back("mcmindeltaenergy"); // 

	opt.allowed.push_back("fixedpos");

	// To specify the number of rotamers for each amino acid type
	opt.allowed.push_back("rotlevel"); 
	opt.allowed.push_back("ALA"); 
	opt.allowed.push_back("ARG");
	opt.allowed.push_back("ASN");
	opt.allowed.push_back("ASP");
	opt.allowed.push_back("CYS");
	opt.allowed.push_back("GLN");
	opt.allowed.push_back("GLU");
	opt.allowed.push_back("GLY");
	opt.allowed.push_back("HSD");
	opt.allowed.push_back("HSE");
	opt.allowed.push_back("HSP");
	opt.allowed.push_back("ILE");
	opt.allowed.push_back("LEU");
	opt.allowed.push_back("LYS");
	opt.allowed.push_back("MET");
	opt.allowed.push_back("PHE");
	opt.allowed.push_back("PRO");
	opt.allowed.push_back("SER");
	opt.allowed.push_back("THR");
	opt.allowed.push_back("TRP");
	opt.allowed.push_back("TYR");
	opt.allowed.push_back("VAL");

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

	opt.logfile = OP.getString("logfile");
	if (!OP.fail()) {
		if(!freopen(opt.logfile.c_str(),"w",stdout)) {
			opt.errorFlag = true;
			opt.errorMessages +=  "Unable to redirect output to " + opt.logfile + "\n";
			return opt;
		}
	}

	opt.verbose = OP.getBool("verbose");
	if(OP.fail()) {
		opt.verbose = false;
		opt.warningMessages += "verbose not specified, using false\n";
		opt.warningFlag = true;
	}
	opt.includeCR = OP.getBool("includecrystalrotamer");
	if(OP.fail()) {
		opt.includeCR = false;
		opt.warningMessages += "includecrystalrotamer not specified, using false\n";
		opt.warningFlag = true;
	}
	opt.onTheFly = OP.getBool("onthefly");
	if(OP.fail()) {
		opt.onTheFly = false;
		opt.warningMessages += "onthefly not specified, using false\n";
		opt.warningFlag = true;
	}

	opt.seed = OP.getBool("seed");
	if(OP.fail()) {
		opt.useTimeToSeed = true;
		opt.warningMessages += "seed not specified, using time based seed\n";
		opt.warningFlag = true;
	} else {
		opt.useTimeToSeed = false;
	}
	opt.outputPDBFile = OP.getString("outputpdbfile");
	if(OP.fail()) {
		opt.outputPDBFile = "";
	}
	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	
	opt.pdbFile= OP.getString("pdbfile");
	if (OP.fail()) {
		opt.errorMessages = "pdbfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.rotlibFile = OP.getString("rotlibfile");
	if (OP.fail()) {
		opt.rotlibFile = SYSENV.getEnv("MSL_ROTLIB");
		opt.warningMessages += "rotlibfile not specified, using " + opt.rotlibFile + "\n";
		opt.warningFlag = true;
	}
	opt.charmmTopFile = OP.getString("charmmtopfile");
	if (OP.fail()) {
		opt.charmmTopFile = SYSENV.getEnv("MSL_CHARMM_TOP");
		opt.warningMessages += "charmmtopfile not specified, using " + opt.charmmTopFile + "\n";
		opt.warningFlag = true;
	}
	opt.charmmParFile = OP.getString("charmmparfile");
	if (OP.fail()) {
		opt.charmmParFile = SYSENV.getEnv("MSL_CHARMM_PAR");
		opt.warningMessages += "charmmparfile not specified, using " + opt.charmmParFile + "\n";
		opt.warningFlag = true;
	}
	opt.hbondParFile = OP.getString("hbondparfile");
	if (OP.fail()) {
		opt.hbondParFile = SYSENV.getEnv("MSL_HBOND_PAR");
		opt.warningMessages += "hbondparfile not specified, using " + opt.hbondParFile + "\n";
		opt.warningFlag = true;
	}
	opt.solvFile = OP.getString("solvfile");
	if (OP.fail()) {
		opt.solvFile = "";
		opt.warningMessages += "solvfile not specified, not using solvation\n";
		opt.warningFlag = true;
	}
	opt.solvent = OP.getString("solvent");
	if (OP.fail()) {
		opt.solvent = "WATER";
		opt.warningMessages += "hbondparfile not specified, using " + opt.solvent + "\n";
		opt.warningFlag = true;
	}
	opt.excludeTerms = OP.getMultiString("excludeenergyterm");
	if(OP.fail()) {
		opt.warningMessages += "excludeenergyterm not specified, using all terms\n";
		opt.warningFlag = true;
		opt.excludeTerms.clear();
	}

	// Algorithms for repack
	opt.runGoldsteinSingles = OP.getBool("rungoldsteinsingles");
	if (OP.fail()) {
		opt.warningMessages += "rungoldsteinsingles not specified, using true\n";
		opt.warningFlag = true;
		opt.runGoldsteinSingles = true;
	}
			
	opt.runGoldsteinPairs = OP.getBool("rungoldsteinpairs");
	if (OP.fail()) {
		opt.warningMessages += "rungoldsteinpairs not specified, using false\n";
		opt.warningFlag = true;
		opt.runGoldsteinPairs = false;
	}

	opt.runSCMF = OP.getBool("runscmf");
	if (OP.fail()) {
		opt.warningMessages += "runscmf not specified, using true\n";
		opt.warningFlag = true;
		opt.runSCMF = true;
	}
			
	opt.runSCMFBiasedMC = OP.getBool("runscmfbiasedmc");
	if (OP.fail()) {
		opt.warningMessages += "runscmfbiasedmc not specified, using true\n";
		opt.warningFlag = true;
		opt.runSCMFBiasedMC = true;
	}
			
	opt.runUnbiasedMC = OP.getBool("rununbiasedmc");
	if (OP.fail()) {
		opt.warningMessages += "rununbiasedmc not specified, using true\n";
		opt.warningFlag = true;
		opt.runUnbiasedMC = true;
	}

	opt.runGreedy = OP.getBool("rungreedy");
	if (OP.fail()) {
		opt.warningMessages += "rungreedy not specified, using false\n";
		opt.warningFlag = true;
		opt.runGreedy = false;
	}

	opt.greedyCycles = OP.getInt("greedycycles"); 
	if(OP.fail()) {
		opt.warningMessages += "greedycycles not specified, using 20\n";
		opt.warningFlag = true;
		opt.greedyCycles = 20;
	}

	// MC Parameters
	opt.startT = OP.getDouble("mcstartt"); 
	if(OP.fail()) {
		opt.warningMessages += "mcstartt not specified, using 2000\n";
		opt.warningFlag = true;
		opt.startT = 2000.0;
	}
	opt.endT = OP.getDouble("mcendt"); 
	if(OP.fail()) {
		opt.warningMessages += "mcendt not specified, using 0.5\n";
		opt.warningFlag = true;
		opt.endT = 0.5;
	}
	opt.nCycles = OP.getInt("mcnumcycles"); 
	if(OP.fail()) {
		opt.warningMessages += "mccycles not specified, using 10000\n";
		opt.warningFlag = true;
		opt.nCycles = 10000;
	}
	string shape = OP.getString("mcshape");
	if(OP.fail()) {
		opt.warningMessages += "mcshape not specified, using EXPONENTIAL\n";
		opt.warningFlag = true;
		shape = "EXPONENTIAL";
	}

	map<string,MonteCarloManager::ANNEALTYPES>::iterator it;
	it = opt.annealShapeMap.find(shape);
	if (it == opt.annealShapeMap.end()){
	  cerr << "ERROR 1111 anneal shape "<<shape<<" is not known.\n";
	  exit(1111);
	} else {
	  opt.shape = it->second;
	}

	/*
	if(shape == "EXPONENTIAL") {
		opt.shape = EXPONENTIAL;
	} else if (shape == "CONSTANT") {
		opt.shape = CONSTANT; 
	} else if (shape == "LINEAR") {
		opt.shape = LINEAR;
	} else if (shape == "SIGMOIDAL") {
		opt.shape = SIGMOIDAL;
	} else if (shape == "SOFT") {
		opt.shape = SOFT;
	} else {
		opt.warningMessages += "mcshape not specified, using EXPONENTIAL\n";
		opt.warningFlag = true;
		opt.shape = EXPONENTIAL;
	}
	*/
	opt.maxReject = OP.getInt("mcmaxreject"); 
	if(OP.fail()) {
		opt.warningMessages += "mcmaxreject not specified, using 2000\n";
		opt.warningFlag = true;
		opt.maxReject = 2000;
	}
	opt.deltaSteps = OP.getInt("mcdeltasteps"); 
	if(OP.fail()) {
		opt.warningMessages += "mcdeltasteps not specified, using 100\n";
		opt.warningFlag = true;
		opt.deltaSteps = 100;
	}
	opt.minDeltaE = OP.getDouble("mcmindeltaenergy");
	if(OP.fail()) {
		opt.warningMessages += "mindeltaE not specified, using 0.01\n";
		opt.warningFlag = true;
		opt.minDeltaE = 0.01;
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

	opt.cuthb = OP.getDouble("cuthb");
	if(OP.fail()) {
		opt.warningMessages += "cuthb not specified, using 10.0\n";
		opt.warningFlag = true;
		opt.cuthb = 10.0;
	}
	
	opt.rotLevel = OP.getString("rotlevel");
	if(OP.fail()) {
		opt.rotLevel = "";
	}

	int num = OP.getInt("ALA"); 
	if(!OP.fail()) {
		opt.numRots["ALA"] = num;
	}
	num = OP.getInt("ARG");
	if(!OP.fail()) {
		opt.numRots["ARG"] = num;
	}
	num = OP.getInt("ASN");
	if(!OP.fail()) {
		opt.numRots["ASN"] = num;
	}
	num = OP.getInt("ASP");
	if(!OP.fail()) {
		opt.numRots["ASP"] = num;
	}
	num = OP.getInt("CYS");
	if(!OP.fail()) {
		opt.numRots["CYS"] = num;
	}
	num = OP.getInt("GLN");
	if(!OP.fail()) {
		opt.numRots["GLN"] = num;
	}
	num = OP.getInt("GLU");
	if(!OP.fail()) {
		opt.numRots["GLU"] = num;
	}
	num = OP.getInt("GLY");
	if(!OP.fail()) {
		opt.numRots["GLY"] = num;
	}
	num = OP.getInt("HSD");
	if(!OP.fail()) {
		opt.numRots["HSD"] = num;
	}
	num = OP.getInt("HSE");
	if(!OP.fail()) {
		opt.numRots["HSE"] = num;
	}
	num = OP.getInt("HSP");
	if(!OP.fail()) {
		opt.numRots["HSP"] = num;
	}
	num = OP.getInt("ILE");
	if(!OP.fail()) {
		opt.numRots["ILE"] = num;
	}
	num = OP.getInt("LEU");
	if(!OP.fail()) {
		opt.numRots["LEU"] = num;
	}
	num = OP.getInt("LYS");
	if(!OP.fail()) {
		opt.numRots["LYS"] = num;
	}
	num = OP.getInt("MET");
	if(!OP.fail()) {
		opt.numRots["MET"] = num;
	}
	num = OP.getInt("PHE");
	if(!OP.fail()) {
		opt.numRots["PHE"] = num;
	}
	num = OP.getInt("PRO");
	if(!OP.fail()) {
		opt.numRots["PRO"] = num;
	}
	num = OP.getInt("SER");
	if(!OP.fail()) {
		opt.numRots["SER"] = num;
	}
	num = OP.getInt("THR");
	if(!OP.fail()) {
		opt.numRots["THR"] = num;
	}
	num = OP.getInt("TRP");
	if(!OP.fail()) {
		opt.numRots["TRP"] = num;
	}
	num = OP.getInt("TYR");
	if(!OP.fail()) {
		opt.numRots["TYR"] = num;
	}
	num = OP.getInt("VAL");
	if(!OP.fail()) {
		opt.numRots["VAL"] = num;
	}

	opt.fixedPos = OP.getStringVector("fixedpos");
	if(!OP.fail()) {
		for(int i = 0; i < opt.fixedPos.size(); i++) {
			opt.fixedPosMap[opt.fixedPos[i]] = true;
		}
	}
	
	// Print Configuration File / Commmand Line Options
	OP.printConfFile();

	
	// return the Options structure
	return opt;

}


