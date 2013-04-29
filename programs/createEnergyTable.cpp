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
#include "SelfPairManager.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"

using namespace std;
using namespace MSL;


// Read a rotamer library and list of environments
// and create a table of energies with conformers along columns and environments along row

string programName = "createEnergyTable";
string programDescription = "This program reads a rotamer library and  list of environments.\n It remodels each environement using conformers from the library and prints out the energies for each combination of environment and conformer as a table with rows = environments and cols = crystalConformation + library conformers.";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "April 02 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct Options {
	string pdbDir; // dir containing the pdbs in the format pdbId.pdb coordinates in PDB format
	string rotlibFile; // file containing the atom coordinates in PDB format
	string envListFile; // file containing the atom coordinates in PDB format
	string charmmTopFile;
	string charmmParFile;

	string hbondParFile;
	string solvParFile;
	string solvent;


	map<string,double> termWeights; // default is use all terms with weight 1.0

	double cuton;
	double cutoff;
	double cutnb;
	double cuthb;
	double vdwRescalingFactor;

	string rawTableFile; 
	string sortedTableFile;
	string configfile; // configFile


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
	cout << " % createEnergyTable \n --pdbdir <pdbDir> \n --envlistfile <envlistfile> \n --rotlibfile <rotlibfile> \n --charmmtopfile <charmmTopFile> \n --charmmparfile <charmmParFile>" << endl ;
	cout << " --rawtablefile <rawtablefile> \n --sortedtablefile <sortedtablefile> \n " << endl;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << "--hbondparfile <hbondParFile=""> --solvparfile <solvParFile=""> --solvent <=WATER>--cuton <cuton=9.0> --cutoff <cutoff=10.0> --cutnb <cutnb=11.0> --cuthb <cuthb=10.0> --vdwrescalingfactor <radius = 1.0x>" << endl;
	cout << " --configfile <configfile> " << endl;
	cout << "--weight <term1,weight> --weight <term2,weight>" << endl;
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
	opt.required.push_back("pdbdir"); // PDB
	opt.required.push_back("envlistfile"); // envlistfile 
	opt.required.push_back("rotlibfile"); // rotamerLibrary 
	opt.required.push_back("charmmtopfile"); // charmmtopology 
	opt.required.push_back("charmmparfile"); // charmmparameter 

	opt.allowed.push_back("hbondparfile"); // hydrogenbond
	opt.allowed.push_back("solvparfile"); // solvation
	opt.allowed.push_back("solvent"); // solvation

	opt.allowed.push_back("rawtablefile"); // repacked structure will be written to this file
	opt.allowed.push_back("sortedtablefile"); // all output will be redirected to this logFile
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("cuton");
	opt.allowed.push_back("cutoff");
	opt.allowed.push_back("cutnb");
	opt.allowed.push_back("cuthb");
	opt.allowed.push_back("vdwrescalingfactor");
	opt.allowed.push_back("weight");

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

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	
	opt.pdbDir= OP.getString("pdbdir");
	if (OP.fail()) {
		opt.errorMessages = "pdbdir not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.envListFile= OP.getString("envlistfile");
	if (OP.fail()) {
		opt.errorMessages = "envlistfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.rotlibFile = OP.getString("rotlibfile");
	if (OP.fail()) {
		opt.errorMessages = "rotlibfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.charmmTopFile = OP.getString("charmmtopfile");
	if (OP.fail()) {
		opt.errorMessages = "charmmtopfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.charmmParFile = OP.getString("charmmparfile");
	if (OP.fail()) {
		opt.errorMessages = "charmmparfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.rawTableFile= OP.getString("rawtablefile");
	if (OP.fail()) {
		opt.errorMessages = "rawtablefile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.sortedTableFile= OP.getString("sortedtablefile");
	if (OP.fail()) {
		opt.errorMessages = "sortedtablefile not specified\n";
		opt.errorFlag = true;
		return opt;
	}

	opt.hbondParFile = OP.getString("hbondparfile");
	if (OP.fail()) {
		opt.warningMessages += "hbondparfile not specified\n";
		opt.warningFlag = true;
		opt.hbondParFile = "";
	}
	opt.solvParFile = OP.getString("solvparfile");
	if (OP.fail()) {
		opt.warningMessages += "solvparfile not specified\n";
		opt.warningFlag = true;
		opt.solvParFile ="";
	}

	opt.solvent = OP.getString("solvent");
	if (OP.fail()) {
		opt.warningMessages += "solvent not specified using WATER\n";
		opt.warningFlag = true;
		opt.solvent ="WATER";
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
	opt.vdwRescalingFactor = OP.getDouble("vdwrescalingfactor");
	if(OP.fail()) {
		opt.warningMessages += "vdwrescalingfactor not specified, using 1.0\n";
		opt.warningFlag = true;
		opt.vdwRescalingFactor = 1.0;
	}

	vector<string> weights = OP.getMultiString("weight");
	if(!OP.fail()) {
		for(int i = 0; i < weights.size(); i++) {
			vector<string> tmp = MslTools::tokenizeAndTrim(weights[i],",");
			if(tmp.size() == 2) {
				opt.termWeights[tmp[0]] = MslTools::toDouble(tmp[1]); 
			} else {
				opt.errorMessages = "syntax error in weight specification [Use Format: --weight CHARMM_ELEC,0.5 --weight CHARMM_BOND,0.5]\n";
				opt.errorFlag = true;
				return opt;
			}
		}
	}
		
	// Print Configuration File / Commmand Line Options
	OP.printConfFile();

	
	// return the Options structure
	return opt;

}

void readEnvListFile (string _fileName, map<string,vector<string> >& _envInPdb) {
	ifstream file;
	file.open(_fileName.c_str());
	string tmpLine;

	// Format : 
	//	1A2J A,44,LEU A,56,ARG
	//      1A3A B,45E,SER
	//	1A2J C,47,LEU

	while(file) {
		getline(file, tmpLine);
		tmpLine = MslTools::uncomment(tmpLine);
		if (tmpLine.length() > 1) {
			vector<string> token = MslTools::tokenizeAndTrim(tmpLine);
			for(int i = 1; i < token.size();i++) {
				_envInPdb[token[0]].push_back(token[i]);
			}
		}
	}
	file.close();

}



void getEnergies(System& _sys, Position* _pPos, RotamerLibrary& _rotlib, vector<double>& _energies) {
	// Build the rotamers onto the given position
	string resName = _pPos->getResidueName();
	Residue& res = _pPos->getCurrentIdentity();
	SystemRotamerLoader sysRot;
	sysRot.setSystem(_sys);
	sysRot.setRotamerLibrary(&_rotlib);

	if(!sysRot.loadRotamers(_pPos,resName,_rotlib.size("",resName),"",true)) {
		cerr << "Cannot load rotamers " << _pPos->getPositionId() << "," << resName << endl;
		return;
	}
	// only the environment position will have multiple conformations/rotamers

	SelfPairManager spm;
	spm.setSystem(&_sys);
	spm.setOnTheFly(true);
	spm.calculateEnergies();

	// set current state to 0th rotamer for all variable positions - here 1
	vector<unsigned int> currState(1,0); 

	_energies.resize(_pPos->getTotalNumberOfRotamers());
	for(unsigned int i = 0; i < _pPos->getTotalNumberOfRotamers();i++) {
		currState[0] = i;
		_energies[i] = spm.getStateEnergy(currState) - spm.getFixedEnergy();
	}

	// Make sure to remove all the loaded conformations
	res.setActiveConformation(0);
	res.removeAllAltConformations();

}

void createTable(Options& _opt,string _pdbName, vector<string> _residueIds, RotamerLibrary& _rotlib, ofstream& _rawOut, ofstream& _sortedOut) {

	System sys;

	// Read all the charmm parameters
	CharmmSystemBuilder csb(sys,_opt.charmmTopFile,_opt.charmmParFile);
	if(_opt.solvParFile != "") {
		// if solvation is required read solvation parameters
		if(!csb.readSolvation(_opt.solvParFile)) {
			cerr << "Unable to read " << _opt.solvParFile << endl;
			exit(0);
		}
		csb.setSolvent(_opt.solvent);
	}

	// Build all the charmm interactions
	csb.setVdwRescalingFactor(_opt.vdwRescalingFactor);
	csb.setBuildNonBondedInteractions(false);
	if(!csb.buildSystemFromPDB(_opt.pdbDir + "/" + _pdbName + ".pdb")) {
		cerr << "Unable to build charmm interactions for " << _opt.pdbDir << "/" << _pdbName << ".pdb" << endl;
		exit(0);
	}
	sys.buildAllAtoms();
	csb.updateNonBonded(_opt.cuton,_opt.cutoff,_opt.cutnb);

	HydrogenBondBuilder hb; 
	// build hbond interactions only if hbondparfile is specified
	if(_opt.hbondParFile != "") {
		hb.setSystem(sys);
		if(!hb.readParameters(_opt.hbondParFile)) {
			cerr << "Unable to read " << _opt.hbondParFile << endl;
			exit(0); 
		}
		if(!hb.buildInteractions(_opt.cuthb)) {
			cerr << "Unable to build hydrogen bond interactions" << endl;
			exit(0);
		}
	}

	EnergySet* eSet = sys.getEnergySet();
	// Reweight each term as specified
	for(map<string,double>::iterator it = (_opt.termWeights).begin(); it != (_opt.termWeights).end(); it++) {
		//cout << it->first << " " << it ->second << endl;
		eSet->setWeight(it->first,it->second);
	}

	for(int i = 0; i < _residueIds.size(); i++) {
		if(!sys.identityExists(_residueIds[i])) {
			cerr << "Missing environment " << _pdbName <<  " " << _residueIds[i] << endl;
			continue;
		}

		Position* pPos  = (sys.getLastFoundIdentity()).getParentPosition();

		vector<double> energies;
		getEnergies(sys,pPos,_rotlib,energies);
		
		// print the pdbname residue Id and crystal energies to both output files 
		_rawOut << _pdbName << " " << _residueIds[i] << " " << energies[0]; 
		_sortedOut << _pdbName << " " << _residueIds[i] << " " << energies[0]; 

		double bestSoFar = energies[1];
		// print raw and sorted energies
		for(int j = 1; j < energies.size(); j++) {
			_rawOut << " " << energies[j];
			if(energies[j] < bestSoFar) {
				bestSoFar = energies[j];
			}
			_sortedOut << " " << bestSoFar;
		}
		_rawOut << endl;
		_sortedOut << endl;

	}
}

int main(int argc, char* argv[]) {

	Options opt = parseOptions(argc, argv);
	if(opt.errorFlag) {
		cout << opt.OPerrors ;
		cout << opt.errorMessages ;
		exit(0);
	}

	if(opt.warningFlag) {
		cout << opt.warningMessages;
	}


	map<string,vector<string> > envInPdb; // map from pdb id to residueIds of environments
	readEnvListFile(opt.envListFile,envInPdb);

	ofstream rawTable;
	rawTable.open(opt.rawTableFile.c_str());

	if(!rawTable.is_open()) {
		cerr << "Unable to open " << opt.rawTableFile << endl;
		exit(0);
	}
	ofstream sortedTable;
	sortedTable.open(opt.sortedTableFile.c_str());

	if(!sortedTable.is_open()) {
		cerr << "Unable to open " << opt.rawTableFile << endl;
		exit(0);
	}

	// Read rotamer library once
	RotamerLibrary rotlib;
	if(!rotlib.readFile(opt.rotlibFile)) {
		cerr << "Unable to read " << opt.rotlibFile << endl;
		exit(0);
	}

	for (map<string,vector<string> >::iterator it = envInPdb.begin(); it != envInPdb.end(); it++) {
		createTable(opt, it->first,it->second,rotlib,rawTable,sortedTable);
	}
	
}
