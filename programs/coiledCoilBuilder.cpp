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


using namespace MSL;
using namespace std;

string programName = "coiledCoilBuilder";
string programDescription = "This program generates Coiled Coils based on a given set of parameters";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.2";
string programDate = "2 June 2010";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

//THINGS TO DO:
/* Option Parser:
	- sequence
	- CC parameters
	- output
	- rotamer library
  */
struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/

	// CC GEOMETRY AND SEQUENCE VARIABLES
	string sequence;
	vector<string> startPos;
	double r0;
	double w0;
	double alpha;
	double r1;
	double w1;
	double phi1;
	double rpr;
	double pitch;
	int nRes;
	string symmetry;
	int N;

	// RUN SETTINGS
	string rotamerSamplingSize;

	// DEAD END ELIMINATION SETTINGS
	bool runDEE;

	// ENUMERATION SETTINGS
	bool runEnum;

	// SELF CONSISTENT MEAN FIELD SETTINGS
	bool runSCMF;

	// MONTECARLO SETTINGS
	bool runMCO;


	// FILE VARIABLES
	string rotamerLibrary; //filename of rotamer library file
	string outputFile; //filename of output file
	string configfile;


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
void saveMin(double _boundE, vector<unsigned int> _stateVec, vector<double> & _minBound, vector<vector<unsigned int> > & _minStates, int _maxSaved, Options & _opt);


/******************************************
 *  
 *  =======  MAIN METHOD  =======
 *
 ******************************************/

int main(int argc, char *argv[]) {
	System sys;
	CoiledCoils cc;
	OptionParser op;
	SelfPairManager spm;
	vector<string> paramList;
	vector<unsigned int> rots;
	vector<string> stringVect;

	cout << "Start..." << endl;

	time_t startTime, endTime;
	double diffTime;

	time (&startTime);


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



	// Map containing list of amino acids and their None, Low, Medium and High rotamer counts
	map<string, map<string, int> > rotLevel;
	map<string, int> tempRotLevel;
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


	//Read in Sequence
	cout << "Reading Polymer sequence" << endl;
	PolymerSequence seq(opt.sequence);

	CharmmSystemBuilder CSB(sys, "/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp");
	CSB.buildSystem(seq);

	// Build Coiled Coil
	cout << "Building Coiled Coil" << endl;
	for (int i = 0; i < opt.startPos.size(); i++) {
		paramList.push_back(opt.startPos[i]);
	}

	if(!cc.primarySequenceToCoiledCoil(opt.r0, opt.rpr, opt.pitch, opt.r1, opt.w1, opt.phi1, 0.0, opt.nRes, opt.symmetry, opt.N, sys, paramList)){
		cerr << "Error.256" << endl;
		exit(0);
	}


	// Add Side Chains
	sys.buildAtoms();

	//Build rotamers
	time_t startRotTime, endRotTime;
	double rotTime;

	time (&startRotTime);

	cout << "Adding Rotamers" << endl;
	SystemRotamerLoader sysRot(sys, opt.rotamerLibrary);

	for (uint i = 0; i < sys.positionSize();i++){
		Position &pos = sys.getPosition(i);

		//string chainId = pos.getChainId();
		//int resNum     = pos.getResidueNumber();

		if (!sysRot.loadRotamers(&pos, "BALANCED-200", sys.getPosition(i).getResidueName(), 0, rotLevel[sys.getPosition(i).getResidueName()][opt.rotamerSamplingSize])) {
			cerr << "Cannot load rotamers " << sys.getPosition(i).getResidueName() << endl;
			exit(1);
		}
				
	}
	spm.setSystem(&sys);
	spm.calculateEnergies();


	time (&endRotTime);
	rotTime = difftime (endRotTime, startRotTime);
	cout << "Rotamer Addition Time: " << rotTime << " seconds" << endl;

	spm.setRunDEE(opt.runDEE);
	spm.setRunEnum(opt.runEnum);
	spm.setRunSCMF(opt.runSCMF);
	spm.setRunMC(opt.runMCO);
	spm.setVerbose(false);

	spm.runOptimizer();

	//vector<vector<unsigned int> > aliveRots = spm.getDEEAliveRotamers();
	//for (int i = 0; i < aliveRots.size(); i++) {
	//	for (int j = 0; j < aliveRots[i].size(); j++){
	//		cout << aliveRots[i][j] << ",";
	//	}
	//	cout << endl;
	//}
	vector<unsigned int> MCOfinal = spm.getMCOstate();

	cout << "MCO accepted state: ";
	for (int i = 0; i < MCOfinal.size(); i++) {
		cout << MCOfinal[i] << ",";
	}
	cout << endl;
	cout << "Energy: " << spm.getStateEnergy(MCOfinal) << endl; 

	sys.setActiveRotamers(spm.getMCOstate());
	sys.writePdb(opt.outputFile);

	time (&endTime);
	diffTime = difftime (endTime, startTime);

	cout << "Output file written to: " << opt.outputFile << endl;
	cout << "Total Time: " << diffTime << " seconds" << endl;
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
	opt.required.push_back("r0");
	opt.required.push_back("w0");
	opt.required.push_back("alpha");
	opt.required.push_back("r1");
	opt.required.push_back("w1");
	opt.required.push_back("phi1");
	opt.required.push_back("rpr");
	opt.required.push_back("pitch");
	opt.required.push_back("nRes");
	opt.required.push_back("symmetry");
	opt.required.push_back("N");
	opt.allowed.push_back("rotamerSamplingSize");
	opt.required.push_back("rotamerLibrary"); //filename of rotamer library file
	opt.required.push_back("outputFile"); //filename of output file
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

	// Print Configuration File / Commmand Line Options
	OP.printConfFile();

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
	opt.outputFile = OP.getString("outputFile");
	if (OP.fail()) {
		opt.errorMessages = "outputFile not specified";
		opt.errorFlag = true;
	}

	opt.rotamerSamplingSize = OP.getString("rotamerSamplingSize");
	opt.runDEE = OP.getBool("runDEE");
	opt.runEnum = OP.getBool("runEnum");
	opt.runSCMF = OP.getBool("runSCMF");
	opt.runMCO = OP.getBool("runMCO");

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
	cout << " % coiledCoilBuilder --sequence <sequence> --startPos <list of start pos> --r0 <superhelical radius> --w0 <superhelical frequency> --alpha <pitch angle> --r1 <alpha helical radius> --w1 <alpha helical frequency> --phi1 <alpha helical phase> --rpr <rise per residue> --pitch <superhelical pitch> --nRes <number of residues> --symmetry <symmetry operator> --N <symmetry N value> --rotamerLibrary <rotLib.txt> --outputFile <nameOfOutPutFile.pdb>" << endl;
	cout << endl;
}
