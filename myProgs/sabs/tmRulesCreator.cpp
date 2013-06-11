#include <iostream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "Transforms.h"
#include "AtomSelection.h"

using namespace MSL;
using namespace std;

string programName = "tmRulesCreator";
string programDescription = "This program checks to see if a residue is compatible at a given position";
string programAuthor = "Benjamin Keymar Mueller";
string programVersion = "1.0.0";
string programDate = "13 February 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {

	string backbonePdb;
	string modelFile;
	string residuesToTest;
	int testPos;

	string topFile;
	string parFile;
	string rotLibFile;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

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

	string configfile;
};

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);

string convertToPolymerSequence(string _seq, int _startNum) {
	// convert a 1 letter _sequence like AIGGG and startNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
	}
	ps = ":{" + MslTools::intToString(_startNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

bool hydrogenBondCheck(vector<string> _parsedGeoInformation, double & _xShiftStart) {
	unsigned int numOfHbonds = (_parsedGeoInformation.size() / 5) - 1;

	if(_parsedGeoInformation.size() < 5) {
		return false;
	}

	// by default lets assign dout based on dmin as follows 
	double dmin = MslTools::toDouble(_parsedGeoInformation[4]); // dmin
	/*
	if(dmin < 7.2) {
		_xShiftStart = 7.4;
	} else if(dmin >= 7.6) {
		_xShiftStart = 7.8;
	} else {
		_xShiftStart = dmin + 0.2;
	}
	*/
	_xShiftStart = dmin + 0.50;

	// of course, if there are hydrogen bonds, then 
	//cout << _sys << endl
	/*;
	for(uint k=6; k < _parsedGeoInformation.size(); k+=5) {
		string donorPos = _parsedGeoInformation[k-1] + "," + _parsedGeoInformation[k]; // Chain,Pos

		string acceptorChain = "";
		if(_parsedGeoInformation[k-1] == "A") {
			acceptorChain = "B";
		}
		else {
			acceptorChain = "A";
		}
		string acceptorPos = acceptorChain + "," + _parsedGeoInformation[k+1];

		unsigned int xShiftStartPosition = 9 + (5 * (numOfHbonds - 4));
		_xShiftStart = MslTools::toDouble(_parsedGeoInformation[xShiftStartPosition]);

	}
	*/
	if(numOfHbonds >= 4) {

		unsigned int xShiftStartPosition = 9 + (5 * (numOfHbonds - 4));
		_xShiftStart = MslTools::toDouble(_parsedGeoInformation[xShiftStartPosition]) - 0.1;
	}

	return true;
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
	
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);

	// if helices are moving symmetrically copy A coor to B and rotate B around the Z axis
	c2Symmetry(_chainA, _chainB);
}

int main(int argc, char *argv[]) {

	time(&startTime);	

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
  
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

	// Read Rotamer Library File
	SystemRotamerLoader sysRot;
	if(!sysRot.readRotamerLibraryFile(opt.rotLibFile)) {
		cerr << "Cannot load rotamer file: " << opt.rotLibFile << endl;
		exit(1);
	}

	// Read in Gly-69 to use as backbone coordinate template
	System gly69;
	if(!gly69.readPdb(opt.backbonePdb)) {
		cout << "Unable to read " << opt.backbonePdb << endl;
		exit(0);
	}

	AtomPointerVector& glyAPV = gly69.getAtomPointers();
	ifstream file;
	file.open(opt.modelFile.c_str());
	if(!file.is_open()) {
		cerr << "Unable to open " << file << endl;
		exit(0);
	}
	string geometryLine;
	vector<string> geometries;
	while(file) {
		getline(file, geometryLine);
		MslTools::uncomment(geometryLine);
		if(geometryLine.length() > 5) {
			geometries.push_back(geometryLine);
		}
	}
	file.close();
	
	for(int geomNum = 0; geomNum < geometries.size(); geomNum++) {
		// Parse parameter file line
		vector<string> parsedGeoInformation = MslTools::tokenize(geometries[geomNum], " ");

		if (parsedGeoInformation.size() % 5 != 0) {
			cerr << "Geometry line is of incompatible length" << endl;
			exit(1);
		}
		// index axialRot crossingAngle zShift xShift ChainDonor DonorResNum AcceptorResNum DonorAtom DistBondBreaks
		// 00001 75 -35.0 1.0 6.8 A 7 8 HA2 9.6 
		// vector length = 5, 0 hbonds
		// vector length = 10, 1 hbonds
		// vector length = 15, 2 hbonds
		double crossingAngle = MslTools::toDouble(parsedGeoInformation[2]);
		double axialRotation = MslTools::toDouble(parsedGeoInformation[1]);
		double zShift = MslTools::toDouble(parsedGeoInformation[3]);
		double xShift = 0.0;
		hydrogenBondCheck( parsedGeoInformation, xShift); // Assign xShift

		// Reference points for Helices
		CartesianPoint ori(0.0,0.0,0.0);
		CartesianPoint zAxis(0.0,0.0,1.0);
		CartesianPoint xAxis(1.0,0.0,0.0);

		string rejectedRes = "";
		//cout <<  "xShift " <<  xShift << endl;

		for (uint i = 0; i < opt.residuesToTest.size(); i++) {

			string residueToChangeTo = opt.residuesToTest.substr(i,1);

			string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
			sequence.replace(opt.testPos-1, 1, residueToChangeTo);
			//cout << "sequence being tried: " << sequence << endl;

			sequence = convertToPolymerSequence(sequence,1); // so that the 4th residue will be the middle one (35th) on the GLY 69 backbone
			PolymerSequence PS(sequence); 

			// Create system with sequence - a string with the following format
			// A:{startingResNum} ALA ILE ...\n
			// B:{startingResNum} ALA ILE ...

			/******************************************************************************
			 *                     === DECLARE SYSTEM ===
			 ******************************************************************************/
			System sys;
			CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
			sysRot.setSystem(sys);
			sysRot.defineRotamerSamplingLevels();

			if(!CSB.buildSystem(PS)) {
				cerr << "Unable to build system from " << sequence << endl;
				exit(0);
			}
			// Set up chain A and chain B atom pointer vectors
			AtomPointerVector &chainA = sys.getChain("A").getAtomPointers();
			AtomPointerVector &chainB = sys.getChain("B").getAtomPointers();

			// Objects used for transformations
			Transforms trans; 
			trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)


			/******************************************************************************
			 *                     === COPY BACKBONE COORDINATES ===
			 ******************************************************************************/
			sys.assignCoordinates(glyAPV,false);		
			sys.buildAtoms();

			/******************************************************************************
			 *                     === INITIAL VARIABLE SET UP ===
			 ******************************************************************************/
			EnergySet* Eset = sys.getEnergySet();

			// Set all terms active, besides Charmm-Elec
			Eset->setAllTermsInactive();
			Eset->setTermActive("CHARMM_VDW", true);

			/******************************************************************************
			 *                  === LOAD ROTAMERS ===
			 ******************************************************************************/
			Position &posA = sys.getPosition(opt.testPos - 1);
			Position &posB = sys.getPosition(69 + opt.testPos - 1); // Hardcoded 69 BEWARE
			

			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				if (!sysRot.loadRotamers(&posA, posA.getResidueName(), "SL95.00")) { 
					cerr << "Cannot load rotamers for " << posA.getResidueName() << endl;
					exit(0);
				}
			}

			if (posB.getResidueName() != "GLY" && posB.getResidueName() != "ALA" && posB.getResidueName() != "PRO") {
				if (!sysRot.loadRotamers(&posB, posB.getResidueName(), "SL95.00")) { 
					cerr << "Cannot load rotamers for " << posB.getResidueName() << endl;
					exit(0);
				}
			}

			/******************************************************************************
			 *                     === COPY BACKBONE COORDINATES ===
			 ******************************************************************************/
			posA.wipeAllCoordinates();
			posB.wipeAllCoordinates();

			sys.assignCoordinates(glyAPV,false);		
			sys.buildAllAtoms();
			sys.saveAltCoor("original");

			/******************************************************************************
			 *                     === MOVE TO TEST GEOMETRY ===
			 ******************************************************************************/
			// Transform helices to monomer position
			transformation(chainA, chainB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, 500, trans);

			SelfPairManager spm;
			double monomerE = 0.0;
			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				spm.setSystem(&sys);
				spm.setRunEnum(true);
				spm.setEnumerationLimit(4000000);
				spm.setRunDEE(false);
				spm.setRunSCMF(false);
				spm.setVerbose(false);

				spm.calculateEnergies();
				spm.runOptimizer();
				monomerE = (spm.getMinBound())[0];
			}
			else {
				monomerE = sys.calcEnergy();
			}

			/******************************************************************************
			 *                     === MOVE TO TEST GEOMETRY ===
			 ******************************************************************************/
			// Transform helices to monomer position
			sys.applySavedCoor("original");
			transformation(chainA, chainB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);

			double bestE = 0.0;
			if (posA.getResidueName() != "GLY" && posA.getResidueName() != "ALA" && posA.getResidueName() != "PRO") {
				spm.calculateEnergies();
				spm.runOptimizer();

				bestE = (spm.getMinBound())[0];
			}
			else {
				bestE = sys.calcEnergy();
			}

			if(bestE > monomerE + 10.0) {
				rejectedRes += residueToChangeTo;
			}
			//cout << "monomerE: " << monomerE << " bestE: " << bestE << endl; 
		}
		if(rejectedRes != "") {
			char tmp[1000];
			sprintf(tmp,"%s,A%2d:![%s],B%2d:![%s]",parsedGeoInformation[0].c_str(),opt.testPos,rejectedRes.c_str(),opt.testPos,rejectedRes.c_str());
			cout << tmp << endl; 
		}
	}

}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);	
	}

	// Rotation matrix for 180 degrees
	Matrix m(3,3,0.0);
	m[0][0] = -1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = -1.0;
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;

	// Rotate chain B around Z axis
	Transforms trans; 
	trans.rotate(_apvB, m);
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

	opt.required.push_back("backbonePdb");
	opt.required.push_back("modelFile");

	opt.required.push_back("testPos");
	opt.required.push_back("residuesToTest");

	opt.required.push_back("topFile");
	opt.required.push_back("parFile");
	opt.required.push_back("rotLibFile");

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
	if (OP.fail()) {
		opt.help = OP.getBool("h");
	}

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


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.backbonePdb = OP.getString("backbonePdb");
	if (OP.fail()) {
		opt.errorMessages = "backbonePdb file not specified";
		opt.errorFlag = true;
	}
	opt.residuesToTest= OP.getString("residuesToTest");
	if (OP.fail()) {
		opt.errorMessages = "residuesToTest not specified";
		opt.errorFlag = true;
	}
	opt.modelFile= OP.getString("modelFile");
	if (OP.fail()) {
		opt.errorMessages = "modelFile not specified";
		opt.errorFlag = true;
	}
	opt.testPos = OP.getInt("testPos");
	if (OP.fail()) {
		opt.errorMessages = "testPos not specified";
		opt.errorFlag = true;
	}
	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		opt.errorMessages = "topFile not specified";
		opt.errorFlag = true;
	}
	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		opt.errorMessages = "parFile not specified";
		opt.errorFlag = true;
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		opt.errorMessages = "rotLibFile not specified";
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
	cout << "This program runs as:" << endl;
	cout << endl;
}
