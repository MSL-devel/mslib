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

string programName = "tmHeteroRulesCreator";
string programDescription = "This program checks to see if a residue is compatible at a given position";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "31 May 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {

	string backbonePdb;
	string modelGeometry;
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
	unsigned int numOfHbonds = ((_parsedGeoInformation.size()-7) / 5);

	if(_parsedGeoInformation.size() < 7) {
		return false;
	}

	// by default lets assign dout based on dmin as follows 
	double dmin = MslTools::toDouble(_parsedGeoInformation[6]); // dmin
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

		unsigned int xShiftStartPosition = 11 + (5 * (numOfHbonds - 4));
		_xShiftStart = MslTools::toDouble(_parsedGeoInformation[xShiftStartPosition]) - 0.1;
	}

	return true;
}
void transformation(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShiftA,double _zShiftB,double _crossingAngle, double _axialRotateA,double _axialRotateB,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis) {
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftACP(0.0, 0.0, _zShiftA);
	_trans.translate(_chainA, zShiftACP);

	CartesianPoint zShiftBCP(0.0, 0.0, _zShiftB);
	_trans.translate(_chainB, zShiftBCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotateA, _origin, _zAxis);
	_trans.rotate(_chainB, _axialRotateB, _origin, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _origin, _xAxis);

	_trans.rotate(_chainB, (_crossingAngle/2.0), _origin, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((_xShift/2.0) * -1.0, 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);

	interDistVect.setCoor((_xShift/2.0) * -1.0, 0.0, 0.0);
	_trans.translate(_chainB, interDistVect);

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
	trans.rotate(_chainB, m);

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

	// Parse parameter file line
	vector<string> parsedGeoInformation = MslTools::tokenize(opt.modelGeometry, " ");

	if ((parsedGeoInformation.size() - 7) % 5 != 0) {
		cerr << "Geometry line is of incompatible length" << endl;
		exit(1);
	}
	// index axialRotA axialRotB crossingAngle zShiftA zShiftB xShift ChainDonor DonorResNum AcceptorResNum DonorAtom DistBondBreaks
	// 184487 10 -32 -11 -0.45 0.45 7.2 A 35 31 HA1 7.3 B 42 42 HA1 7.5 B 46 46 HA2 7.5 A 28 27 HA2 8.0 B 35 35 HA2 8.1 A 39 38 HA2 8.2
	// vector length = 7, 0 hbonds
	// vector length = 12, 1 hbonds
	// vector length = 17, 2 hbonds
	double crossingAngle = MslTools::toDouble(parsedGeoInformation[3]);
	double axialRotationA = MslTools::toDouble(parsedGeoInformation[1]);
	double axialRotationB = MslTools::toDouble(parsedGeoInformation[2]);
	double zShiftA = MslTools::toDouble(parsedGeoInformation[4]);
	double zShiftB = MslTools::toDouble(parsedGeoInformation[5]);
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
//		void transformation(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShiftA,double _zShiftB,double _crossingAngle, double _axialRotateA,double _axialRotateB,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis) {
		transformation(chainA, chainB, zShiftA, zShiftB, crossingAngle, axialRotationA, axialRotationB, 500,trans,ori,zAxis,xAxis);

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
		transformation(chainA, chainB, zShiftA, zShiftB, crossingAngle, axialRotationA, axialRotationB, xShift, trans, ori, zAxis, xAxis);

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
	opt.required.push_back("modelGeometry");

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
	opt.modelGeometry= OP.getString("modelGeometry");
	if (OP.fail()) {
		opt.errorMessages = "modelGeometry not specified";
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
