
#include <iostream>
#include <stdio.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"


using namespace MSL;
using namespace std;

string programName = "genHeteroUniverse";
string programDescription = "This program creates coiled coils with different configurations and measures energy. Finds the Dmin and also the Douts for each hydrogen bond at Dmin.";
string programAuthor = "Sabareesh Subramaniam,Benjamin K Mueller";
string programVersion = "1.0.0";
string programDate = "12 April 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {
	string pdbFile;
	string helicalAxisPdbFile;
	string output;

	double xShiftStart;
	double zShiftAStart;
	double zShiftBStart;
	double axialRotAStart;
	double axialRotBStart;
	double crossingAngleStart;
	//vector<int> rotCount;	
	double xShiftEnd;
	double zShiftAEnd;
	double zShiftBEnd;
	double axialRotAEnd;
	double axialRotBEnd;
	double crossingAngleEnd;
	//vector<int> rotCount;	
	double xShiftSteps;
	double zShiftASteps;
	double zShiftBSteps;
	double axialRotASteps;
	double axialRotBSteps;
	double crossingAngleSteps;
	//vector<int> rotCount;	

	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;

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

void transformCoiledCoil(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShiftA,double _zShiftB,double _crossingAngle, double _axialRotateA,double _axialRotateB,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, AtomPointerVector& _axisA, AtomPointerVector& _axisB) {
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
	_trans.rotate(_axisA,  (_crossingAngle/2.0), _origin, _xAxis);

	_trans.rotate(_chainB, (_crossingAngle/2.0), _origin, _xAxis);
	_trans.rotate(_axisB,  (_crossingAngle/2.0), _origin, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((_xShift/2.0) * -1.0, 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	interDistVect.setCoor((_xShift/2.0) * -1.0, 0.0, 0.0);
	_trans.translate(_chainB, interDistVect);
	_trans.translate(_axisB, interDistVect);

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
	trans.rotate(_axisB, m);

}


unsigned int getNumberOfInterHelicalHbonds(vector<Interaction*>& hbondInteractions) {
	unsigned int numHbonds = 0;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			numHbonds++;
		}
	}
	return numHbonds;
}
map<string,double> getInterHelicalHbondInfo(vector<Interaction*>& hbondInteractions) {
	map<string,double> info;

	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			info[atoms[0]->getAtomId() + ":" + atoms[2]->getAtomId()] = e;
		}
	}
	return info;
}

void printInfo(double axialRotateA, double axialRotateB,double crossingAngle,double zShiftA,double zShiftB, vector<double> xShifts,vector<double> energies,vector<map<string,double> > hBonds) {
	int minIdx = 0; // into xShifts, energies,hBonds
	for(int i = 1; i < energies.size(); i++) {
		if(energies[i] < energies[minIdx]) {
			minIdx = i;
		}
	}

	// get the hbonds at minIdx and see where they get broken
	map<string,double> hmin = hBonds[minIdx];
	map<string,double> douts;

	for(map<string,double>::iterator it = hmin.begin(); it != hmin.end(); it++) {
		string name = it->first;
		douts[name] = xShifts[minIdx];
		for(int i = minIdx+1; i < hBonds.size(); i++) {
			if(hBonds[i].find(name) == hBonds[i].end()) {
				douts[name] = xShifts[i];
				break;
			}
		}
	}
	char name[1000];
	sprintf(name,"model_%02.0f_%02.0f_%02.0f_%03.2f_%03.2f_%02.1f",axialRotateA,axialRotateB,crossingAngle,zShiftA,zShiftB,xShifts[minIdx]);
	cout << name << " " << hBonds[minIdx].size();
	for(int i = minIdx + 1; i < xShifts.size(); i++) {
		char str[100];
		sprintf(str,"%02.1f",xShifts[i]);
		for(map<string,double>::iterator it = hmin.begin(); it != hmin.end(); it++) {
			if(douts[it->first] == xShifts[i]) {
				cout << " " << it->first << " " << str;
			}
		}
	}
	cout << endl;

}

int main(int argc, char *argv[]) {

	time(&startTime);	
	
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

	/******************************************************************************
	 *                     === SYSTEM SETUP ===
	 ******************************************************************************/
	
	// Declare System
	System sys;

	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile, opt.solvFile);
	//CharmmSystemBuilder CSB(sys, "/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp", "/library/mslib/toppar/charmm/solvpar22.inp");
	CSB.setSolvent("CHEX");
	CSB.setVdwRescalingFactor(1.00);

	// Read in PDB File
	//if(!CSB.buildSystemFromPDB("/data00/bkmueller/gpaProject/createHelixAandB.pdb")) {
	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
		cerr << "Unable to build system" << endl;
		exit(0);
	}

	// Build system
	sys.buildAllAtoms();

	// Add hydrogen bonds
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	//HydrogenBondBuilder hb(sys, "/data00/bkmueller/dataFiles/hbondlist_nonCanon_adj.txt");
	hb.buildInteractions(50);

	// Redirect Output
	string filename = opt.output;
	freopen (filename.c_str(), "w", stdout);

	// Set up APVs
	AtomPointerVector &chainA = sys.getChain("A").getAtomPointers();
	AtomPointerVector &chainB = sys.getChain("B").getAtomPointers();



	/******************************************************************************
	 *                     === MONTE CARLO SET UP ===
	 ******************************************************************************/

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxisPdbFile);
	//helicalAxis.readPdb("/data00/bkmueller/gpaProject/mutationRunWithPhiPsi/helicalAxis.pdb");

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/

	EnergySet* Eset = sys.getEnergySet();
	AtomSelection sel(sys.getAtomPointers());

	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	
	sys.saveCoor("initialState");
	helicalAxis.saveCoor("initialState");

	EnergySet* pESet = sys.getEnergySet();
	vector<Interaction*> hbondInteractions = (*(pESet->getEnergyTerms()))["SCWRL4_HBOND"];

	//cout << "Number of hbond interactions "<< hbondInteractions.size() << endl;
	// Reference points to set up Helical starting postions
	CartesianPoint origin(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);


	vector<double> xShifts;
	for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
		xShifts.push_back(xShift);
	}


	for(double axialRotateA = opt.axialRotAStart; axialRotateA < opt.axialRotAEnd; axialRotateA += opt.axialRotASteps) {
		for(double axialRotateB = opt.axialRotBStart; axialRotateB < opt.axialRotBEnd; axialRotateB += opt.axialRotBSteps) {
			for(double crossingAngle = opt.crossingAngleStart; crossingAngle < opt.crossingAngleEnd; crossingAngle += opt.crossingAngleSteps ) {
				for(double zShiftA = opt.zShiftAStart; zShiftA < opt.zShiftAEnd; zShiftA += opt.zShiftASteps) {
					for(double zShiftB = opt.zShiftBStart; zShiftB < opt.zShiftBEnd; zShiftB += opt.zShiftBSteps) {
						vector<map<string,double> >  hBonds;
						vector<double> energies;
						for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
							sys.applySavedCoor("initialState");
							helicalAxis.applySavedCoor("initialState");
							transformCoiledCoil(chainA,chainB,zShiftA,zShiftB,crossingAngle,axialRotateA,axialRotateB,xShift,trans,origin,zAxis,xAxis,axisA,axisB);
							double thisEnergy = sys.calcEnergy();
							//unsigned int numHbonds = getNumberOfInterHelicalHbonds(hbondInteractions);
							hBonds.push_back(getInterHelicalHbondInfo(hbondInteractions));
							energies.push_back(thisEnergy);
							//cout << "axialRotationA: " << axialRotateA  << "axialRotationB: " << axialRotateB  <<" crossingAngle: " << crossingAngle << " zShiftA: " << zShiftB << " zShiftB: " << zShiftB <<  " xShift: " << xShift <<  endl;
						}
						printInfo(axialRotateA,axialRotateB,crossingAngle,zShiftA,zShiftB,xShifts,energies,hBonds);
					}
				}
			}
		}
	}



	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

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

	opt.required.push_back("pdbFile");
	opt.required.push_back("helicalAxisPdbFile");
	opt.required.push_back("output");

	opt.required.push_back("xShiftStart");
	opt.required.push_back("xShiftEnd");
	opt.required.push_back("xShiftSteps");

	opt.required.push_back("zShiftAStart");
	opt.required.push_back("zShiftAEnd");
	opt.required.push_back("zShiftASteps");

	opt.required.push_back("zShiftBStart");
	opt.required.push_back("zShiftBEnd");
	opt.required.push_back("zShiftBSteps");

	opt.required.push_back("axialRotAStart");
	opt.required.push_back("axialRotAEnd");
	opt.required.push_back("axialRotASteps");

	opt.required.push_back("axialRotBStart");
	opt.required.push_back("axialRotBEnd");
	opt.required.push_back("axialRotBSteps");

	opt.required.push_back("crossingAngleStart");
	opt.required.push_back("crossingAngleEnd");
	opt.required.push_back("crossingAngleSteps");

	//opt.required.push_back("rotCount");

	opt.required.push_back("topFile");
	opt.required.push_back("parFile");
	opt.required.push_back("solvFile");
	opt.required.push_back("hBondFile");

	//opt.equivalent.push_back(vector<string>());
	//opt.equivalent.back().push_back("v");
	//opt.equivalent.back().push_back("version");
	//opt.equivalent.push_back(vector<string>());
	//opt.equivalent.back().push_back("h");
	//opt.equivalent.back().push_back("help");

	//opt.defaultArgs.push_back("configfile");


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


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages = "pdb file not specified";
		opt.errorFlag = true;
	}
	opt.helicalAxisPdbFile = OP.getString("helicalAxisPdbFile");
	if (OP.fail()) {
		opt.errorMessages = "helicalAxisPdbFile file not specified";
		opt.errorFlag = true;
	}
	opt.output = OP.getString("output");
	if (OP.fail()) {
		opt.errorMessages = "output not specified";
		opt.errorFlag = true;
	}
	opt.xShiftStart = OP.getDouble("xShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "xShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.xShiftEnd = OP.getDouble("xShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "xShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.xShiftSteps = OP.getDouble("xShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "xShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.zShiftAStart = OP.getDouble("zShiftAStart");
	if (OP.fail()) {
		opt.errorMessages = "zShiftAStart not specified";
		opt.errorFlag = true;
	}
	opt.zShiftAEnd = OP.getDouble("zShiftAEnd");
	if (OP.fail()) {
		opt.errorMessages = "zShiftAEnd not specified";
		opt.errorFlag = true;
	}
	opt.zShiftASteps = OP.getDouble("zShiftASteps");
	if (OP.fail()) {
		opt.errorMessages = "zShiftASteps not specified";
		opt.errorFlag = true;
	}

	opt.zShiftBStart = OP.getDouble("zShiftBStart");
	if (OP.fail()) {
		opt.errorMessages = "zShiftBStart not specified";
		opt.errorFlag = true;
	}
	opt.zShiftBEnd = OP.getDouble("zShiftBEnd");
	if (OP.fail()) {
		opt.errorMessages = "zShiftBEnd not specified";
		opt.errorFlag = true;
	}
	opt.zShiftBSteps = OP.getDouble("zShiftBSteps");
	if (OP.fail()) {
		opt.errorMessages = "zShiftBSteps not specified";
		opt.errorFlag = true;
	}

	opt.axialRotAStart = OP.getDouble("axialRotAStart");
	if (OP.fail()) {
		opt.errorMessages = "axialRotAStart not specified";
		opt.errorFlag = true;
	}
	opt.axialRotAEnd = OP.getDouble("axialRotAEnd");
	if (OP.fail()) {
		opt.errorMessages = "axialRotAEnd not specified";
		opt.errorFlag = true;
	}
	opt.axialRotASteps = OP.getDouble("axialRotASteps");
	if (OP.fail()) {
		opt.errorMessages = "axialRotASteps not specified";
		opt.errorFlag = true;
	}

	opt.axialRotBStart = OP.getDouble("axialRotBStart");
	if (OP.fail()) {
		opt.errorMessages = "axialRotBStart not specified";
		opt.errorFlag = true;
	}
	opt.axialRotBEnd = OP.getDouble("axialRotBEnd");
	if (OP.fail()) {
		opt.errorMessages = "axialRotBEnd not specified";
		opt.errorFlag = true;
	}
	opt.axialRotBSteps = OP.getDouble("axialRotBSteps");
	if (OP.fail()) {
		opt.errorMessages = "axialRotBSteps not specified";
		opt.errorFlag = true;
	}

	opt.crossingAngleStart = OP.getDouble("crossingAngleStart");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleStart not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleEnd = OP.getDouble("crossingAngleEnd");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleEnd not specified";
		opt.errorFlag = true;
	}
	opt.crossingAngleSteps = OP.getDouble("crossingAngleSteps");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleSteps not specified";
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
	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		opt.errorMessages = "solvFile not specified";
		opt.errorFlag = true;
	}
	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		opt.errorMessages = "hBondFile not specified";
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
	cout << " % genHomoUniverse --pdbFile <pdbfile> --helicalAxisPdbFile <filename> --output <filename>" << endl;
	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
	cout << " --zShiftAStart <double> --zShiftAEnd <double> --zShiftASteps <double> " << endl;
	cout << " --zShiftBStart <double> --zShiftBEnd <double> --zShiftBSteps <double> " << endl;
	cout << " --axialRotAStart <double> --axialRotAEnd <double> --axialRotASteps <double> " << endl;
	cout << " --axialRotBStart <double> --axialRotBEnd <double> --axialRotBSteps <double> " << endl;
	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
	cout << endl;
}
