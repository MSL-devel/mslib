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
#include <iostream>

#include "CoiledCoils.h"
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "Enumerator.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"


using namespace MSL;
using namespace std;

string programName = "coiledCoilBuilder";
string programDescription = "This program generates Coiled Coils based on a given set of parameters";
string programAuthor = "Benjamin Keymar Mueller";
string programVersion = "1.0.0";
string programDate = "31 August 2010";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

 // variables to store time at start and end of application
time_t start_time, finish_time;


struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/

	// CC GEOMETRY AND SEQUENCE VARIABLES
	string sequence;
	vector<string> startPos;
	string startRegPos;
	double r0_start;
	double r0_end;
	double r0_step;
	double r1_start;
	double r1_end;
	double r1_step;
	double w1_start;
	double w1_end;
	double w1_step;
	double phi1_start;
	double phi1_end;
	double phi1_step;
	double rpr_start;
	double rpr_end;
	double rpr_step;
	double pitch_start;
	double pitch_end;
	double pitch_step;
	int nRes;
	string symmetry;
	int N;

	//int testNo;

	// RUN SETTINGS
	string rotamerSamplingSize;
	int seed;
	bool setState;
	vector<unsigned int> rotamerStates;
	string rotamerLibrary; //filename of rotamer library file

	// DEAD END ELIMINATION SETTINGS
	bool runDEE;

	// ENUMERATION SETTINGS
	bool runEnum;

	// SELF CONSISTENT MEAN FIELD SETTINGS
	bool runSCMF;

	// MONTECARLO SETTINGS
	bool runMCO;


	// FILE VARIABLES
	string outputdir;  // the directory with the output for the run
	string outputfile; //filename of output file
	string configfile;
	string rerunConfFile;
	string resultsFile; // sorted energies and states
	string resultsFileKey; // sorted energies and states


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
	string rerunConf; // data for a configuration file that would rerun the job as the current run
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
void saveMin(double _boundE, vector<unsigned int> _stateVec, vector<double> & _minBound, vector<vector<unsigned int> > & _minStates, int _maxSaved, Options & _opt);


/******************************************
 *  
 *  =======  MAIN METHOD  =======
 *
 ******************************************/

int main(int argc, char *argv[]) {
	// store the start time
	time(&start_time);

	System sys;
	CoiledCoils cc;
	OptionParser op;
	//vector<string> startingPositions;
	vector<unsigned int> rots;
	vector<string> stringVect;

	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are	
	 *  not always required                            
	 ******************************************************************************/
	Options defaults;                                  
   	defaults.outputdir = "coiledCoilBuilder";
                                                        
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
	} else {
		/*****************************************
		 *  IF THE OUTPUT DIR HAS A RANDOMLY GENERATED BIT (%r)
		 *  THE FOLLOWING ENSURE THAT IT WILL NOT PICK A DIRECTORY
		 *  THAT ALREADY EXISTS
		 *****************************************/
		string prevOutputdir = opt.outputdir;
		opt.outputdir = MslTools::outputFileNameParser(opt.outputdir);
		while (true) {
			if (prevOutputdir == opt.outputdir) {
				// the outputdir doesn't have randomly generated parts
				break;
			}
			ifstream test_fs;
			test_fs.open(opt.outputdir.c_str());
			if (test_fs.is_open()) {
				// directory exists, try again
				prevOutputdir = opt.outputdir;
				opt.outputdir = MslTools::outputFileNameParser(opt.outputdir);
				test_fs.close();
			} else {
				test_fs.close();
				break;
			}
		}

		opt.outputfile = opt.outputdir + "/" + opt.outputfile;
		opt.resultsFile = opt.outputdir + "/sortedEnergies.txt";
		opt.resultsFileKey = opt.outputdir + "/sortedEnergiesKey.txt";
		opt.rerunConfFile = opt.outputdir + "/rerun_conf.txt";

		cout << "Output in " << opt.outputdir << endl;
		MslTools::mkNestedDir(opt.outputdir, 0755); 
		//chdir(opt.outputdir.c_str());
		//log_fs.open(opt.logfile.c_str());
		//if (log_fs.fail()) {
		//	cerr << "ERROR: Cannot write " << opt.logfile << endl;
		//	exit(1);
		//}
	}

	/******************************************************************************
	 *                            === OUTPUT REDIRECTION ===
	 *
	 *   The standard output can be redirected to a file 
	 *   (option "--coutRedirected <filename>"
	 *   
	 *   if we redirectint, open the output file
	 ******************************************************************************/
	if (opt.coutRedirected) {
		string ofile = opt.outputfile + ".out";
		string efile = opt.outputfile + ".err";
		opt.cout_fs->open(ofile.c_str());
		if (opt.cout_fs->fail()) {
			cerr << "ERROR: cannot open for writing standard output file " << ofile << endl;
			exit(1);
		}
		opt.cerr_fs->open(efile.c_str());
		if (opt.cerr_fs->fail()) {
			cerr << "ERROR: cannot open for writing standard error file " << efile << endl;
			exit(1);
		}

		/****************************************************
		 * NOTE: REMEMBER TO CALL THE endjob() FUNCTION AT THE END
		 ****************************************************/
		// save the streambuf of cout on a pointer
		opt.coutbuf = cout.rdbuf();
		opt.cerrbuf = cerr.rdbuf();

		// JOIN THE TWO A NOT USE coutsbuf?
		// get the streambuf of the file
		opt.coutsbuf = opt.cout_fs->rdbuf();
		opt.cerrsbuf = opt.cerr_fs->rdbuf();

		// assign the streambuf to cout
		cout.rdbuf(opt.coutsbuf);
		cerr.rdbuf(opt.cerrsbuf);
		
		//opt.coutRedirected = true;
	}
	
	cout << "Date:" << (string)ctime(&start_time) << endl;
	cout << "Program: " << programName << endl;
	cout << "Program description: " << programDescription << endl;
	cout << "Program author: " << programAuthor << endl;
	cout << "Program version: " << programVersion << " " << programDate << endl;
	cout << "MSL version: " << mslVersion << " " << mslDate << endl;

	/************************************************************
	 *   Output a configuration file that would rerun the 
	 *   program with the same options of the current run (this is
	 *   useful when mixing options between the configuration file and
	 *   the command line
	 ************************************************************/
	ofstream rerunConf_fs;
	rerunConf_fs.open(opt.rerunConfFile.c_str());
	if (rerunConf_fs.fail()) {
		cerr << "ERROR: Cannot write " << opt.rerunConfFile << endl;
		exit(1);
	}
	rerunConf_fs << opt.rerunConf;
	rerunConf_fs.close();

	// Total Time
	time_t startTotalTime, endTotalTime;
	double totalTime;

	time (&startTotalTime);

	// Map for register positions
	vector<string> registerPos;
		registerPos.push_back("H");
		registerPos.push_back("N");
		registerPos.push_back("N");
		registerPos.push_back("H");
		registerPos.push_back("M");
		registerPos.push_back("N");
		registerPos.push_back("M");

	// Map containing list of amino acids and their None, Low, Medium and High rotamer counts
	map<string, map<string, int> > rotLevel;
		rotLevel["ALA"]["N"] = 0;
		rotLevel["ALA"]["L"] = 0;
		rotLevel["ALA"]["M"] = 0;
		rotLevel["ALA"]["H"] = 0;

		rotLevel["ARG"]["N"] = 0;
		rotLevel["ARG"]["L"] = 6;
		rotLevel["ARG"]["M"] = 12;
		rotLevel["ARG"]["H"] = 24;

		rotLevel["ASN"]["N"] = 0;
		rotLevel["ASN"]["L"] = 4;
		rotLevel["ASN"]["M"] = 7;
		rotLevel["ASN"]["H"] = 14;
		
		rotLevel["ASP"]["N"] = 0;
		rotLevel["ASP"]["L"] = 3;
		rotLevel["ASP"]["M"] = 5;
		rotLevel["ASP"]["H"] = 10;

		rotLevel["CYS"]["N"] = 0;
		rotLevel["CYS"]["L"] = 3;
		rotLevel["CYS"]["M"] = 5;
		rotLevel["CYS"]["H"] = 14;

		rotLevel["GLU"]["N"] = 0;
		rotLevel["GLU"]["L"] = 5;
		rotLevel["GLU"]["M"] = 10;
		rotLevel["GLU"]["H"] = 20;
                                          
		rotLevel["GLN"]["N"] = 0;
		rotLevel["GLN"]["L"] = 6;
		rotLevel["GLN"]["M"] = 12;
		rotLevel["GLN"]["H"] = 14;
                                          
		rotLevel["GLY"]["N"] = 0;
		rotLevel["GLY"]["L"] = 0;
		rotLevel["GLY"]["M"] = 0;
		rotLevel["GLY"]["H"] = 0;
                                          
		rotLevel["HSD"]["N"] = 0;
		rotLevel["HSD"]["L"] = 7;
		rotLevel["HSD"]["M"] = 15;
		rotLevel["HSD"]["H"] = 30;
                                          
		rotLevel["ILE"]["N"] = 0;
		rotLevel["ILE"]["L"] = 5;
		rotLevel["ILE"]["M"] = 10;
		rotLevel["ILE"]["H"] = 20;
                                       
		rotLevel["LEU"]["N"] = 0;
		rotLevel["LEU"]["L"] = 5;
		rotLevel["LEU"]["M"] = 10;
		rotLevel["LEU"]["H"] = 20;
                                          
		rotLevel["LYS"]["N"] = 0;
		rotLevel["LYS"]["L"] = 6;
		rotLevel["LYS"]["M"] = 12;
		rotLevel["LYS"]["H"] = 24;
                                          
		rotLevel["MET"]["N"] = 0;
		rotLevel["MET"]["L"] = 4;
		rotLevel["MET"]["M"] = 7;
		rotLevel["MET"]["H"] = 14;
                                          
		rotLevel["PHE"]["N"] = 0;
		rotLevel["PHE"]["L"] = 7;
		rotLevel["PHE"]["M"] = 15;
		rotLevel["PHE"]["H"] = 30;
                                          
		rotLevel["PRO"]["N"] = 0;
		rotLevel["PRO"]["L"] = 0;
		rotLevel["PRO"]["M"] = 0;
		rotLevel["PRO"]["H"] = 0;
                                       
		rotLevel["SER"]["N"] = 0;
		rotLevel["SER"]["L"] = 3;
		rotLevel["SER"]["M"] = 5;
		rotLevel["SER"]["H"] = 10;
                                          
		rotLevel["THR"]["N"] = 0;
		rotLevel["THR"]["L"] = 3;
		rotLevel["THR"]["M"] = 5;
		rotLevel["THR"]["H"] = 10;
                                          
		rotLevel["TRP"]["N"] = 0;
		rotLevel["TRP"]["L"] = 7;
		rotLevel["TRP"]["M"] = 15;
		rotLevel["TRP"]["H"] = 30;
                                          
		rotLevel["TYR"]["N"] = 0;
		rotLevel["TYR"]["L"] = 7;
		rotLevel["TYR"]["M"] = 15;
		rotLevel["TYR"]["H"] = 30;
                                          
		rotLevel["VAL"]["N"] = 0;
		rotLevel["VAL"]["L"] = 3;
		rotLevel["VAL"]["M"] = 5;
		rotLevel["VAL"]["H"] = 10;


	// Read in Sequence
	cout << "Reading Polymer sequence" << endl;
	PolymerSequence seq(opt.sequence);

	// Build Sequence
	CharmmSystemBuilder CSB(sys, "/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp");
	CSB.setDielectricConstant(80.0);
	CSB.setVdwRescalingFactor(0.95);
	CSB.buildSystem(seq);

	// Add Hydrogen Bonds
	HydrogenBondBuilder HBB(sys, "/data00/bkmueller/dataFiles/hbondlist.txt");
	HBB.buildInteractions(10.0);


	vector<unsigned int> enumCounter;

	// Create Enumerator to Loop over states
	// Adding 0.00001 is to avoid rounding errors and have the loop go to the end step, this will be fixed later
	vector<double> r0_vec;
	for (double r0 = opt.r0_start; r0 <= opt.r0_end + 0.00001; r0+=opt.r0_step) {
		r0_vec.push_back(r0);
	}
	enumCounter.push_back(r0_vec.size());
	vector<double> r1_vec;
	for (double r1 = opt.r1_start; r1 <= opt.r1_end + 0.00001; r1+=opt.r1_step) {
		r1_vec.push_back(r1);
	}
	enumCounter.push_back(r1_vec.size());
	vector<double> w1_vec;
	for (double w1 = opt.w1_start; w1 <= opt.w1_end + 0.00001; w1+=opt.w1_step) {
		w1_vec.push_back(w1);
	}
	enumCounter.push_back(w1_vec.size());
	vector<double> phi1_vec;
	for (double phi1 = opt.phi1_start; phi1 <= opt.phi1_end + 0.00001; phi1+=opt.phi1_step) {
		phi1_vec.push_back(phi1);
	}
	enumCounter.push_back(phi1_vec.size());
	vector<double> rpr_vec;
	for (double rpr = opt.rpr_start; rpr <= opt.rpr_end + 0.00001; rpr+=opt.rpr_step) {
		rpr_vec.push_back(rpr);
	}
	enumCounter.push_back(rpr_vec.size());
	vector<double> pitch_vec;
	for (double pitch = opt.pitch_start; pitch <= opt.pitch_end + 0.00001; pitch+=opt.pitch_step) {
		pitch_vec.push_back(pitch);
	}
	enumCounter.push_back(pitch_vec.size());

	Enumerator coilStates(enumCounter);
	vector<double> savedEnergyVector;
	vector<vector<unsigned int> > savedMCOState;
	vector<unsigned int> energyVectorIndex;

	if (opt.setState) {
		cout << "Setting specific rotameric state..." << endl;
		SelfPairManager spm;

		// Build Coiled Coil
		sys.wipeAllCoordinates();
		vector<string> startingPositions;

		cout << endl;
		cout << "Building Coiled Coil" << endl;
		for (int i = 0; i < opt.startPos.size(); i++) {
			startingPositions.push_back(opt.startPos[i]);
		}

		if(!cc.setSystemToCoiledCoil(opt.r0_start, opt.rpr_start, opt.pitch_start, opt.r1_start, opt.w1_start, opt.phi1_start, 0.0, opt.nRes, opt.symmetry, opt.N, sys, startingPositions)){
			cerr << "Error.256" << endl;
			exit(0);
		}

		// Add Side Chains
		sys.buildAtoms();

		//Build rotamers
		cout << "Adding Rotamers" << endl;
		vector<unsigned int> specificState = opt.rotamerStates;

		SystemRotamerLoader sysRot(sys, opt.rotamerLibrary);

		for (uint i = 0; i < sys.positionSize(); i++){
			Position &pos = sys.getPosition(i);

			//if (!sysRot.loadRotamers(&pos, "BALANCED-200", sys.getPosition(i).getResidueName(), 0, specificState[i])) {
			if (!sysRot.loadRotamers(&pos, sys.getPosition(i).getResidueName(), 0, specificState[i])) {
				cerr << "ERROR 1: Cannot load rotamers " << sys.getPosition(i).getResidueName() << endl;
				exit(1);
			}
		}

		// Calculating Energies
		vector<unsigned int> testVector;
		for (int k = 0; k < sys.positionSize(); k++) {
			testVector.push_back(k);
		}
		sys.setVariablePositions(testVector);

		spm.setSystem(&sys);
		spm.calculateEnergies(); 

		// Set State
		sys.setActiveRotamers(specificState);

		cout << "Energy: " << setiosflags(ios::fixed) << setprecision(10) << spm.getStateEnergy(specificState) << endl;
		//sys.writePdb(opt.outputFile);
		//cout << "Output file written to: " << opt.outputFile << endl;

		exit(0);
	}

	// START LOOP!
	for (int k=0; k < coilStates.size(); k++) {
		SelfPairManager spm;
		time_t startTime, endTime;
		double diffTime;

		time (&startTime);

		// Build Coiled Coil
		sys.wipeAllCoordinates();
		vector<string> startingPositions;
		//startingPositions.clear();

		cout << endl;
		cout << "Building Coiled Coil" << endl;
		for (int i = 0; i < opt.startPos.size(); i++) {
			startingPositions.push_back(opt.startPos[i]);
		}

		if(!cc.setSystemToCoiledCoil(r0_vec[coilStates[k][0]], rpr_vec[coilStates[k][4]], pitch_vec[coilStates[k][5]], r1_vec[coilStates[k][1]], w1_vec[coilStates[k][2]], phi1_vec[coilStates[k][3]], 0.0, opt.nRes, opt.symmetry, opt.N, sys, startingPositions)){
			cerr << "Error.256" << endl;
			exit(0);
		}

		cout << "#################################" << endl;
		cout << "CYCLE: " << k << " / " << coilStates.size()-1 << endl;
		cout << "r0: " << setiosflags(ios::fixed) << setprecision(2) << r0_vec[coilStates[k][0]] << endl;
		cout << "r1: " << r1_vec[coilStates[k][1]] << endl;
		cout << "w1: " << w1_vec[coilStates[k][2]] << endl;
		cout << "phi1: " << phi1_vec[coilStates[k][3]] << endl;
		cout << "rpr: " << rpr_vec[coilStates[k][4]] << endl;
		cout << "pitch: " << pitch_vec[coilStates[k][5]] << endl;
		cout << "#################################" << endl;


		// Add Side Chains
		sys.buildAtoms();
		//Build rotamers
		time_t startRotTime, endRotTime;
		double rotTime;

		time (&startRotTime);

		cout << "Adding Rotamers" << endl;
		SystemRotamerLoader sysRot(sys, opt.rotamerLibrary);

		string rotSamplingLevel;
		int rotamerAdderPosition = 0;
		if (!opt.startRegPos.compare("A")) { 
			rotSamplingLevel = registerPos[0]; 
			rotamerAdderPosition = 0;
		}
		else if (!opt.startRegPos.compare("B")) { 
			rotSamplingLevel = registerPos[1];
			rotamerAdderPosition = 1;
		}
		else if (!opt.startRegPos.compare("C")) { 
			rotSamplingLevel = registerPos[2];
			rotamerAdderPosition = 2;
		}
		else if (!opt.startRegPos.compare("D")) { 
			rotSamplingLevel = registerPos[3];
			rotamerAdderPosition = 3;
		}
		else if (!opt.startRegPos.compare("E")) { 
			rotSamplingLevel = registerPos[4];
			rotamerAdderPosition = 4;
		}
		else if (!opt.startRegPos.compare("F")) { 
			rotSamplingLevel = registerPos[5];
			rotamerAdderPosition = 5;
		}
		else if (!opt.startRegPos.compare("G")) { 
			rotSamplingLevel = registerPos[6];
			rotamerAdderPosition = 6;
		}
		else {
			cerr << "Incorrect register position entry" << endl;
			exit(0);
		}

		for (uint i = 0; i < sys.positionSize();i++){
			Position &pos = sys.getPosition(i);

			rotSamplingLevel = registerPos[rotamerAdderPosition];

			//if (!sysRot.loadRotamers(&pos, "BALANCED-200", sys.getPosition(i).getResidueName(), 0, rotLevel[sys.getPosition(i).getResidueName()][rotSamplingLevel])) {
			if (!sysRot.loadRotamers(&pos, sys.getPosition(i).getResidueName(), 0, rotLevel[sys.getPosition(i).getResidueName()][rotSamplingLevel])) {
				cerr << "ERROR 2: Cannot load rotamers " << sys.getPosition(i).getResidueName() << endl;
				exit(2);
			}

			if (rotamerAdderPosition < registerPos.size()-1) {
				rotamerAdderPosition++;
			}
			else {
				rotamerAdderPosition = 0;
			}

		}

		char c[1000];
		sprintf(c, "coil-%06u.pdb", k);
		//char d[1000];
		//sprintf(d, "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f", r0_vec[coilStates[k][0]], w0_vec[coilStates[k][1]], r1_vec[coilStates[k][2]], w1_vec[coilStates[k][3]], phi1_vec[coilStates[k][4]], rpr_vec[coilStates[k][5]], pitch_vec[coilStates[k][6]]);
		sys.writePdb(c);
		cout << "Written PDB " << c << endl;

		// Rotamer Addition Time
		time (&endRotTime);
		rotTime = difftime (endRotTime, startRotTime);
		cout << "Rotamer Addition Time: " << rotTime << " seconds" << endl;

		// Calculating Energies
		time (&startRotTime);
		spm.setSystem(&sys);
		spm.calculateEnergies(); 
		time (&endRotTime);
		rotTime = difftime (endRotTime, startRotTime);
		cout << "Energy Calc time: " << rotTime << " seconds" << endl;

		// Run Optimization
		spm.setRunDEE(opt.runDEE);
		spm.setRunEnum(opt.runEnum);
		spm.setRunSCMF(opt.runSCMF);
		spm.setRunSCMFBiasedMC(opt.runMCO);
		spm.setVerbose(false);
		spm.seed(opt.seed);

		spm.runOptimizer();

		// Final Output
		vector<unsigned int> MCOfinal = spm.getBestSCMFBiasedMCState();

		cout << "MCO accepted state: ";
		for (int i = 0; i < MCOfinal.size(); i++) {
			cout << MCOfinal[i] << ",";
		}
		cout << endl;
		cout << "Energy: " << setiosflags(ios::fixed) << setprecision(10) << spm.getStateEnergy(MCOfinal) << endl; 

		// Save the energy, MCO state and the Enumerator state
		savedEnergyVector.push_back(spm.getStateEnergy(MCOfinal));
		savedMCOState.push_back(MCOfinal);
		energyVectorIndex.push_back(k);

		time (&endTime);
		diffTime = difftime (endTime, startTime);
		cout << "Seed value: " << spm.getSeed() << endl;
		cout << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
	}

	// Sort Energy Low to High
	MslTools::quickSortWithIndex(savedEnergyVector, energyVectorIndex);
	vector<unsigned int> variablePos = sys.getVariablePositions();

	ofstream resultsFile_fs;
	resultsFile_fs.open(opt.resultsFile.c_str());
	if (resultsFile_fs.fail()) {
		cerr << "ERROR: Cannot write " << opt.resultsFile<< endl;
		exit(1);
	}

	char c[1000];
	// index of the variable positions
	resultsFile_fs << "#";
	for (int i = 0; i < variablePos.size(); i++) {
		sprintf(c, "%4u", variablePos[i]);
		resultsFile_fs << c;
	}
	resultsFile_fs << endl;
	
	for (int i=0; i < savedEnergyVector.size(); i++) {
		// energy
		sprintf(c, "%25.4f", savedEnergyVector[i]);
		resultsFile_fs << c;

		// coiled coil parameters
		sprintf(c, "%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", r0_vec[coilStates[energyVectorIndex[i]][0]], r1_vec[coilStates[energyVectorIndex[i]][1]], w1_vec[coilStates[energyVectorIndex[i]][2]], phi1_vec[coilStates[energyVectorIndex[i]][3]], rpr_vec[coilStates[energyVectorIndex[i]][4]], pitch_vec[coilStates[energyVectorIndex[i]][5]]);
		resultsFile_fs << c;

		// rotamer index of each variable position
		string line;
		for (int j = 0; j < savedMCOState[energyVectorIndex[i]].size(); j++) {
			sprintf(c, "%4u", savedMCOState[energyVectorIndex[i]][j]);
			line += c;
		}
		resultsFile_fs << line << endl;

		/*
		opt.resultsFile += "NO. " + MslTools::intToString(i) + "\n";
		opt.resultsFile += "Energy: " + savedEnergyVector[i] + "\n";
		opt.resultsFile += "MCO state: ";
		for (int j = 0; j < savedMCOState[energyVectorIndex[i]].size(); j++) {
			opt.resultsFile += savedMCOState[energyVectorIndex[i]][j] + ",";
		}
		opt.resultsFile += "\n";
		opt.resultsFile += "r0: " + setiosflags(ios::fixed) + setprecision(2) + r0_vec[coilStates[energyVectorIndex[i]][0]];
		opt.resultsFile += "\tr1: " + r1_vec[coilStates[energyVectorIndex[i]][1]];
		opt.resultsFile += "\tw1: " + w1_vec[coilStates[energyVectorIndex[i]][2]];
		opt.resultsFile += "\tphi1: " + phi1_vec[coilStates[energyVectorIndex[i]][3]];
		opt.resultsFile += "\trpr: " + rpr_vec[coilStates[energyVectorIndex[i]][4]];
		opt.resultsFile += "\tpitch: " + pitch_vec[coilStates[energyVectorIndex[i]][5]] + "\n";
		*/
	}

	//resultsFile_fs << opt.resultsFile;
	resultsFile_fs.close();

	// Key File Output
	resultsFile_fs.open(opt.resultsFileKey.c_str());
	if (resultsFile_fs.fail()) {
		cerr << "ERROR: Cannot write " << opt.resultsFileKey << endl;
		exit(1);
	}
	resultsFile_fs << "Key for the result file:" << endl;
	resultsFile_fs << "" << endl;
	resultsFile_fs << "First line (commented out): indeces of the variable positions" << endl;
	resultsFile_fs << "" << endl;
	resultsFile_fs << "Energy line fields (space separated):" << endl;
	resultsFile_fs << "1) Energy" << endl;
	resultsFile_fs << "2) r0" << endl;
	resultsFile_fs << "3) r1" << endl;
	resultsFile_fs << "4) w1" << endl;
	resultsFile_fs << "5) phi1" << endl;
	resultsFile_fs << "6) rpr" << endl;
	resultsFile_fs << "7) pitch" << endl;
	resultsFile_fs << "8) Rotameric state of the first variable position" << endl;
	resultsFile_fs << "9) Rotameric state of the second variable position" << endl;
	resultsFile_fs << "10) ... etc (up to the number of variable positions)" << endl;
	resultsFile_fs << "" << endl;

	resultsFile_fs.close();

	// End of Program comments
	cout << endl;
	time (&endTotalTime);
	totalTime = difftime(endTotalTime, startTotalTime);
	cout << "Total Time: " << totalTime << endl;
	cout << "Finish." << endl;
}


/*****************************************
 *       OTHER FUNCTIONS
 ****************************************/
void saveMin(double _boundE, vector<unsigned int> _stateVec, vector<double> & _minBound, vector<vector<unsigned int> > & _minStates, int _maxSaved, Options & _opt) {

	// case the list is emtpy
	if (_minBound.size() == 0) {
		_minBound.push_back(_boundE);
		_minStates.push_back(_stateVec);
		return;
	} else if (_minBound.size() == _maxSaved) {
		if (_boundE >= _minBound.back()) {
			return;
		}
	}
	
	// make sure that we have not yet saved the state
	vector<vector<unsigned int> >::iterator v;

	// ... else try to fit it in the middle
	v = _minStates.begin();
	vector<double>::iterator b;
	int counter = 0;
	for (b = _minBound.begin(); b != _minBound.end(); b++, v++, counter++) {
		double Enew = _boundE;
		double Eold = *b;
		if (Enew <= Eold) {
			if (Enew == Eold) {
				// make sure that we have not yet saved the state
				bool identical = true;
				for (unsigned int i=0; i<_stateVec.size(); i++) {
					if (_stateVec[i] != (*v)[i]) {
						identical = false;
						break;
					}
				}

				if (identical) {
					return;
				}
			}
			_minBound.insert(b, _boundE);
			_minStates.insert(v, _stateVec);
			// trim the list if oversized
			if (_minBound.size() > _maxSaved) {
				b = _minBound.end()-1;
				v = _minStates.end()-1;
				_minBound.erase(b);
				_minStates.erase(v);
			}
			return;
		}
	}

	// ... or stick it at the end if the list is short
	_minBound.push_back(_boundE);
	_minStates.push_back(_stateVec);

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

	opt.required.push_back("sequence");
	opt.required.push_back("startPos");
	opt.required.push_back("startRegPos");

	opt.required.push_back("r0_start");
	opt.required.push_back("r0_end");
	opt.required.push_back("r0_step");
	opt.required.push_back("r1_start");
	opt.required.push_back("r1_end");
	opt.required.push_back("r1_step");
	opt.required.push_back("w1_start");
	opt.required.push_back("w1_end");
	opt.required.push_back("w1_step");
	opt.required.push_back("phi1_start");
	opt.required.push_back("phi1_end");
	opt.required.push_back("phi1_step");
	opt.required.push_back("rpr_start");
	opt.required.push_back("rpr_end");
	opt.required.push_back("rpr_step");
	opt.required.push_back("pitch_start");
	opt.required.push_back("pitch_end");
	opt.required.push_back("pitch_step");

	opt.required.push_back("nRes");
	opt.required.push_back("symmetry");
	opt.required.push_back("N");

	//opt.required.push_back("testNo");
	opt.allowed.push_back("setState");
	opt.allowed.push_back("rotamerStates");

	opt.allowed.push_back("rotamerSamplingSize");
	opt.allowed.push_back("seed");
	opt.required.push_back("rotamerLibrary"); //filename of rotamer library file
	opt.allowed.push_back("outputdir");
	opt.allowed.push_back("outputfile"); //filename of output file
	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("v"); // -v is equivalent to --version
	opt.allowed.push_back("help"); // --help
	opt.allowed.push_back("h"); // -h is equivalent to --help
	opt.allowed.push_back("configfile");

	opt.allowed.push_back("runDEE");
	opt.allowed.push_back("runEnum");
	opt.allowed.push_back("runSCMF");
	opt.allowed.push_back("runMCO");

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
	 *  OUTPUT DIR AND FILES
	 *****************************************/
	opt.outputdir = OP.getString("outputdir");
	if (OP.fail() || opt.outputdir == "") {
		opt.outputdir = defaults.outputdir;
	}
	
	opt.outputfile = OP.getString("outputfile");
	if (OP.fail()) {
		opt.outputfile = defaults.outputfile;
	}
	if (opt.outputfile != "" && opt.outputfile != "STDOUT") {
		opt.coutRedirected = true;
		opt.cout_fs = new ofstream();
		opt.cerr_fs = new ofstream();
	} else {
		opt.coutRedirected = false;
		opt.cout_fs = NULL;
		opt.cerr_fs = NULL;
	}


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

	opt.sequence= OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages = "sequence not specified";
		opt.errorFlag = true;
	}
	opt.startPos = OP.getStringVector("startPos");
	if (OP.fail()) {
		opt.errorMessages = "startPos not specified";
		opt.errorFlag = true;
	}
	opt.startRegPos = OP.getString("startRegPos");
	if (OP.fail()) {
		opt.errorMessages = "startRegPos not specified";
		opt.errorFlag = true;
	}
	opt.r0_start = OP.getDouble("r0_start");
	if (OP.fail()) {
		opt.errorMessages = "r0 not specified";
		opt.errorFlag = true;
	}
	opt.r0_end = OP.getDouble("r0_end");
	if (OP.fail()) {
		opt.errorMessages = "r0 not specified";
		opt.errorFlag = true;
	}
	opt.r0_step = OP.getDouble("r0_step");
	if (OP.fail()) {
		opt.errorMessages = "r0 not specified";
		opt.errorFlag = true;
	}
	opt.r1_start = OP.getDouble("r1_start");
	if (OP.fail()) {
		opt.errorMessages = "r1 not specified";
		opt.errorFlag = true;
	}
	opt.r1_end = OP.getDouble("r1_end");
	if (OP.fail()) {
		opt.errorMessages = "r1 not specified";
		opt.errorFlag = true;
	}
	opt.r1_step = OP.getDouble("r1_step");
	if (OP.fail()) {
		opt.errorMessages = "r1 not specified";
		opt.errorFlag = true;
	}
	opt.w1_start = OP.getDouble("w1_start");
	if (OP.fail()) {
		opt.errorMessages = "w1 not specified";
		opt.errorFlag = true;
	}
	opt.w1_end = OP.getDouble("w1_end");
	if (OP.fail()) {
		opt.errorMessages = "w1 not specified";
		opt.errorFlag = true;
	}
	opt.w1_step = OP.getDouble("w1_step");
	if (OP.fail()) {
		opt.errorMessages = "w1 not specified";
		opt.errorFlag = true;
	}
	opt.phi1_start = OP.getDouble("phi1_start");
	if (OP.fail()) {
		opt.errorMessages = "phi1 not specified";
		opt.errorFlag = true;
	}
	opt.phi1_end = OP.getDouble("phi1_end");
	if (OP.fail()) {
		opt.errorMessages = "phi1 not specified";
		opt.errorFlag = true;
	}
	opt.phi1_step = OP.getDouble("phi1_step");
	if (OP.fail()) {
		opt.errorMessages = "phi1 not specified";
		opt.errorFlag = true;
	}
	opt.rpr_start = OP.getDouble("rpr_start");
	if (OP.fail()) {
		opt.errorMessages = "rpr not specified";
		opt.errorFlag = true;
	}
	opt.rpr_end = OP.getDouble("rpr_end");
	if (OP.fail()) {
		opt.errorMessages = "rpr not specified";
		opt.errorFlag = true;
	}
	opt.rpr_step = OP.getDouble("rpr_step");
	if (OP.fail()) {
		opt.errorMessages = "rpr not specified";
		opt.errorFlag = true;
	}
	opt.pitch_start = OP.getDouble("pitch_start");
	if (OP.fail()) {
		opt.errorMessages = "pitch not specified";
		opt.errorFlag = true;
	}
	opt.pitch_end = OP.getDouble("pitch_end");
	if (OP.fail()) {
		opt.errorMessages = "pitch not specified";
		opt.errorFlag = true;
	}
	opt.pitch_step = OP.getDouble("pitch_step");
	if (OP.fail()) {
		opt.errorMessages = "pitch not specified";
		opt.errorFlag = true;
	}
	opt.nRes = OP.getInt("nRes");
	if (OP.fail()) {
		opt.errorMessages = "number of residues not specified";
		opt.errorFlag = true;
	}
	opt.symmetry = OP.getString("symmetry");
	if (OP.fail()) {
		opt.errorMessages = "symmetry not specified";
		opt.errorFlag = true;
	}
	opt.N = OP.getInt("N");
	if (OP.fail()) {
		opt.errorMessages = "N not specified";
		opt.errorFlag = true;
	}
	opt.rotamerLibrary = OP.getString("rotamerLibrary");
	if (OP.fail()) {
		opt.errorMessages = "rotamerLibrary not specified";
		opt.errorFlag = true;
	}
	opt.outputfile = OP.getString("outputfile");
	if (OP.fail()) {
		opt.errorMessages = "outputfile not specified";
		opt.errorFlag = true;
	}

	/*
	opt.testNo = OP.getInt("testNo");
	if (OP.fail()) {
		opt.errorMessages = "testNo not specified";
		opt.errorFlag = true;
	}
	*/

	opt.rotamerSamplingSize = OP.getString("rotamerSamplingSize");
	opt.seed = OP.getInt("seed");
	opt.runDEE = OP.getBool("runDEE");
	opt.runEnum = OP.getBool("runEnum");
	opt.runSCMF = OP.getBool("runSCMF");
	opt.runMCO = OP.getBool("runMCO");
	opt.setState = OP.getBool("setState");
	opt.rotamerStates = OP.getUnsignedIntVector("rotamerStates");

	opt.rerunConf = "########################################################\n";
	opt.rerunConf += "#  This configuration file was automatically generated,\n";
	opt.rerunConf += "#  it will rerun this job with the same options. Run as:\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#  Run as:\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#    % " + programName + " --configfile " + opt.rerunConfFile + "\n";
	opt.rerunConf += "#\n";
	opt.rerunConf += "#  Job started on " + (string)ctime(&start_time);
	opt.rerunConf += "#  on host " + opt.host + ", path " + opt.pwd + "\n";
	if (opt.seed == 0) {
		opt.rerunConf += "#  seed " + MslTools::intToString(opt.seed) + " (time based)\n";
	} else {
		opt.rerunConf += "#  seed " + MslTools::intToString(opt.seed) + "\n";
	} 
	opt.rerunConf += "########################################################\n";
	opt.rerunConf += "\n";
	opt.rerunConf += OP.getConfFile();

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

	string sequence;
	vector<string> startPos;
	string CCparametersFile; //filename of CC params file
	string rotamerLibrary; //filename of rotamer library file
	string outputFile; //filename of output file

void help(Options defaults) {
       	cout << "Run  as:" << endl;
	cout << " % coiledCoilBuilder --sequence <sequence> --startPos <list of start pos> --startRegPos <starting heptad position> --r0 <superhelical radius> --w0 <superhelical frequency> --r1 <alpha helical radius> --w1 <alpha helical frequency> --phi1 <alpha helical phase> --rpr <rise per residue> --pitch <superhelical pitch> --nRes <number of residues> --symmetry <symmetry operator> --N <symmetry N value> --rotamerLibrary <rotLib.txt> --outputFile <nameOfOutPutFile.pdb>" << endl;
	cout << endl;
}
