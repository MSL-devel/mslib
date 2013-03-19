#include <iostream>
#include <fstream>

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
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"
#include "EZpotentialBuilder.h"


using namespace MSL;
using namespace std;

string programName = "predictHelixOligomer";
string programDescription = "This program repacks a dimer of given crossingAngle, axialRotation, zShift with xShift from dout to dmin and prints the best model";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "0.0.0";
string programDate = "19 March 2013";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

static SysEnv ENV;

struct Options {
	// Required
	string tmSequence;

	// optional
	string backboneCrd;
	string logFile;

	string helixGeoLine;

	string hisProtState;

	string topFile;
	string parFile;
	string hBondFile;
	string rotLibFile;

	bool verbose;
	int greedyCycles;
	int seed;

	// protein information (optional)
	string uniprotName;
	string uniprotAccession;

	bool printTermEnergies;

	int startResNum;

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

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "backboneCrd " << _op.backboneCrd << endl;
	_fout << "logFile " << _op.logFile << endl;

	_fout << "tmSequence " << _op.tmSequence << endl;

	_fout << "helixGeoLine " << _op.helixGeoLine << endl;
	_fout << "hisProtState " << _op.hisProtState << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "uniprotName " << _op.uniprotName << endl;
	_fout << "uniprotAccession " << _op.uniprotAccession << endl;
	_fout << "startResNum " << _op.startResNum << endl;

	_fout << "printTermEnergies " << _op.printTermEnergies << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	
	/* Faster code 
	for (uint i=0; i < _apvA.size(); i++) {
			_apvB[i]->copyAllCoor(*_apvA[i]);
			vector<CartesianPoint*>& bCoors = _apvB[i]->getAllCoor();

			for(uint j = 0; j < bCoors.size(); j++) {
				bCoors[j]->setX(0 - bCoors[j]->getX());
				bCoors[j]->setY(0 - bCoors[j]->getY());
			}
					
		}
	*/

	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
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

void setRotamerLevelByPosition(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	for (uint i=0; i < _apvA.size(); i++) {
		if(_apvA[i]->getParentResidue()->getResidueName() != "GLY" && _apvA[i]->getParentResidue()->getResidueName() != "ALA" && _apvA[i]->getParentResidue()->getResidueName() != "PRO") {
			if (_apvA[i]->getName() == "CB") {
				unsigned int counter = 0;
				for (uint j=0; j < _apvB.size(); j++) {
					if (_apvB[j]->getName() == "CA") {
						if (_apvA[i]->distance(*_apvB[j]) < 10.0) {
							counter++;
						}
					}
				}
				if (counter >= 7) { _apvA[i]->getParentResidue()->setRotamerSamplingLevel("SL95.00"); }
				else if(counter > 4) { _apvA[i]->getParentResidue()->setRotamerSamplingLevel("SL85.00"); }
				else if (counter > 2) { _apvA[i]->getParentResidue()->setRotamerSamplingLevel("SL80.00"); }
				else { _apvA[i]->getParentResidue()->setRotamerSamplingLevel("SL70.00"); }
			}
		}
	}
	for (uint i=0; i < _apvB.size(); i++) {
		if(_apvB[i]->getParentResidue()->getResidueName() != "GLY" && _apvB[i]->getParentResidue()->getResidueName() != "ALA" && _apvB[i]->getParentResidue()->getResidueName() != "PRO") {
			if (_apvB[i]->getName() == "CB") {
				unsigned int counter = 0;
				for (uint j=0; j < _apvA.size(); j++) {
					if (_apvA[j]->getName() == "CA") {
						if (_apvB[i]->distance(*_apvA[j]) < 10.0) {
							counter++;
						}
					}
				}
				if (counter >= 7) { _apvB[i]->getParentResidue()->setRotamerSamplingLevel("SL95.00"); }
				else if(counter > 4) { _apvB[i]->getParentResidue()->setRotamerSamplingLevel("SL85.00"); }
				else if (counter > 2) { _apvB[i]->getParentResidue()->setRotamerSamplingLevel("SL80.00"); }
				else { _apvB[i]->getParentResidue()->setRotamerSamplingLevel("SL70.00"); }
			}
		}
	}
}

string convertToPolymerSequence(string _seq, int _startResNum, string _hisProtState) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " " + _hisProtState;
		} else {
			ps = ps + " " + resName;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}


void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
	
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);


	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV) {

	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());
	//fout << x << " " << y << " " << pt << " " << caApV.size() << endl;

	for(int i = 0; i < _apV.size(); i++) {
		CartesianPoint& pt = _apV[i]->getCoor();
		pt.setZ(pt.getZ() +  zShift);
	}

}

vector<string> getInterHelicalHbonds(EnergySet* & _ESet) {
	unsigned int numHbonds = 0;
	// Why are we doing this?
	//_ESet->setAllTermsInactive();
	//_ESet->setTermActive("SCWRL4_HBOND", true);
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];

	vector<string> hbonds;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			//_fout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
			hbonds.push_back(atoms[0]->getAtomOfIdentityId() + ":" + atoms[2]->getAtomOfIdentityId() + "=" + MslTools::doubleToString(atoms[0]->distance(*atoms[2])));
			numHbonds++;
		}
	}
	// Why are we doing this?
	//_ESet->setAllTermsActive();
	//_ESet->setTermActive("CHARMM_ELEC", false);
	return hbonds;
}

map<string,double> getEnergyByTerm(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
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

	ofstream fout;
	fout.open(opt.logFile.c_str());

	if(!fout.is_open()) {
		cerr << "Unable to open " << opt.logFile << endl;
		exit(0);
	}

	printOptions(opt, fout);

	// Import sequence, determine length
	unsigned int sequenceLength = opt.tmSequence.length(); 
	// cant do sequences of size less than 4
	if(sequenceLength < 4) {
		cerr << "Sequence " << opt.tmSequence << "is too small (should be >= 4 AA long)" << endl; 
		exit(0);
	}

	string modelledTMSeq ="" ; // replace P with A
	string prolineMask = "";
	for(int i = 0; i < opt.tmSequence.length(); i++) {
		char AA = opt.tmSequence[i];
		if(AA == 'P') {
			modelledTMSeq += "A";
			prolineMask += "1";
		} else {
			modelledTMSeq += AA;
			prolineMask += "0";
		}
	}
	

	//cout << modelledTMSeq << endl;

	// A:30 B:30 55 -8 -0.75 9 9 A,45,GLN,HE22;B,45,GLN,OE1

	vector<string> geoLine = MslTools::tokenizeAndTrim(opt.helixGeoLine);
	int threadStart = MslTools::toInt(geoLine[0].substr(2));


	string sequence = convertToPolymerSequence(modelledTMSeq,threadStart,opt.hisProtState); 
	PolymerSequence PS(sequence); 

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << sequence << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd);
	if(!cRead.read()) {
		fout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();
	AtomPointerVector& glyAPV = cRead.getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(opt.seed);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	hb.buildInteractions(30);


	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/

	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";
	

	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	helicalAxis.saveCoor("originState");

	// Declare SelfPairManager and Set Seed
	SelfPairManager spm;
	spm.seed(RNG1.getSeed()); 

	spm.setOnTheFly(true);
	spm.setVerbose(opt.verbose);

	/******************************************************************************
	 *                  === LOAD ROTAMERS & SET-UP SPM ===
	 ******************************************************************************/
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!sysRot.loadRotamers(&pos, pos.getResidueName(), "SL95.00")) { 
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

//	spm.saveEnergiesByTerm(true);
	spm.setSystem(&sys);
	sys.saveAltCoor("originState");

//	fout <<  "Monomer Energy " << monomerEnergy << endl;
//	cout <<  "Monomer Energy " << computeMonomerEnergy(sys, trans, RNG1, sysRot, spm) << endl;
//	cout << spm.getSummary(spm.getMinStates()[0]) << endl;

	/******************************************************************************
	 *                  === READ IN RULES ===
	 ******************************************************************************/


	// We will need a format converter to convert to PDB names
	FormatConverter fc;

	// transform using parameters in the geoLine
	// A:30 B:30 55 -8 -0.75 9 9 A,45,GLN,HE22;B,45,GLN,OE1

	double crossingAngle = MslTools::toDouble(geoLine[2]);
	double axialRotation = MslTools::toDouble(geoLine[3]);
	double zShift = MslTools::toDouble(geoLine[4]);
	double dmin = MslTools::toDouble(geoLine[5]);
	double dout = MslTools::toDouble(geoLine[6]);
	double xShift = dout;
	transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);

	double bestEnergy = MslTools::doubleMax;
	double bestXShift = dout;

	while(xShift >= dmin ) {
		// repack
		//sys.writePdb("/tmp/pdb_" + MslTools::doubleToString(xShift) + ".pdb");
		setRotamerLevelByPosition(apvChainA,apvChainB);
		spm.calculateEnergies();
		spm.runGreedyOptimizer(opt.greedyCycles);
		double currentEnergy = spm.getMinBound()[0];
		fout << "XSHIFT " << xShift << " ENERGY " << currentEnergy << endl;
		if(currentEnergy < bestEnergy) {
			sys.setActiveRotamers(spm.getMinStates()[0]);
			sys.calcEnergy();
			sys.saveCoor("savedBestState");
			bestXShift = dmin;
			bestEnergy = currentEnergy;
		}
		sys.applySavedCoor("originState");
		helicalAxis.applySavedCoor("originState");
		xShift -= 0.1;
		transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
	}

	sys.applySavedCoor("savedBestState");
	double finalEnergy = sys.calcEnergy();
	vector<string> hbonds = getInterHelicalHbonds(Eset);
	if(hbonds.size() > 0) {
		fout << "HBONDS ";
	}
	for(int i = 0; i < hbonds.size(); i++) {
		fout << hbonds[i] << " ";
	}
	fout << endl;

	if(opt.printTermEnergies) {
		map<string,double> eByTerm = getEnergyByTerm(Eset);
		fout << "Best Energy " << finalEnergy << endl;
		for(map<string,double>::iterator it = eByTerm.begin(); it != eByTerm.end(); it++) {
			fout << it->first << " " << it->second << endl;	
		}
	}

	fout << "GEOM " << opt.helixGeoLine << " " << bestXShift << endl;
	// renumber using startResNum, moveToCentreOfMass and print pdbs and 
	// TODO even .inp files for pymol
	AtomPointerVector& allAtoms = sys.getAtomPointers();
	moveZCenterOfCAMassToOrigin(allAtoms);

	chainA.renumberChain(opt.startResNum);
	chainB.renumberChain(opt.startResNum);

	string outName = "struct_" + geoLine[0].substr(2) + "_" + geoLine[2] + "_" + geoLine[3] + "_" + geoLine[4] + "_" + MslTools::doubleToString(bestXShift) + ".pdb"; 
	if(!sys.writePdb(outName)) {
		cerr << "Unable to write " << outName << endl;
		exit(0);
	}

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	fout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
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

	opt.required.push_back("tmSequence");

	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("logFile");

	opt.required.push_back("helixGeoLine");
	opt.allowed.push_back("hisProtState");

	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("hBondFile");
	opt.allowed.push_back("rotLibFile");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	opt.allowed.push_back("uniprotName");
	opt.allowed.push_back("uniprotAccession");
	opt.allowed.push_back("startResNum");

	opt.allowed.push_back("printTermEnergies");

	opt.allowed.push_back("configfile");

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

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.tmSequence = OP.getString("tmSequence");
	if (OP.fail()) {
		opt.errorMessages += "tmSequence (1 letter aa) not specified\n";
		opt.errorFlag = true;
	}

	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.warningMessages += "backboneCrd file not specified using /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd\n";
		opt.warningFlag = true;
		opt.backboneCrd = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd";
	}
	opt.helixGeoLine = OP.getString("helixGeoLine");
	if (OP.fail()) {
		opt.errorMessages += "helixGeoLine not specified\n";
		opt.errorFlag = true;
	}

	opt.uniprotName = OP.getString("uniprotName");
	if (OP.fail()) {
		opt.uniprotName = "PROTEIN_UNK";
		opt.warningMessages += "uniprotName not specified using " + opt.uniprotName + "\n";
		opt.warningFlag = true;
	}
	opt.uniprotAccession = OP.getString("uniprotAccession");
	if (OP.fail()) {
		opt.uniprotAccession = "P00000";
		opt.warningMessages += "uniprotAccession not specified using " + opt.uniprotAccession + "\n";
		opt.warningFlag = true;
	}

	opt.logFile = OP.getString("logFile");
	if (OP.fail()) {
		opt.logFile = opt.uniprotAccession + ".log";
		opt.warningMessages += "logFile not specified using " + opt.logFile + "\n";
		opt.warningFlag = true;
	}

	opt.hisProtState = OP.getString("hisProtState");	
	if(OP.fail()) {
		opt.hisProtState = "HSE";
		opt.warningMessages += "hisProtState not specified using HSE\n";
		opt.warningFlag = true;
	}
	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(ENV.isDefined(envVar)) {
			opt.topFile = ENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}

	}
	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(ENV.isDefined(envVar)) {
			opt.parFile = ENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_PAR";
		if(ENV.isDefined(envVar)) {
			opt.hBondFile = ENV.getEnv(envVar);
			opt.warningMessages += "hBondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hBondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_ROTLIB";
		if(ENV.isDefined(envVar)) {
			opt.rotLibFile = ENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}

	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using 1\n";
		opt.warningFlag = true;
		opt.startResNum = 1;
	}
	opt.printTermEnergies = OP.getBool("printTermEnergies");
	if (OP.fail()) {
		opt.printTermEnergies = true;
		opt.warningMessages += "printTermEnergies not specified using true\n";
		opt.warningFlag = true;
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
	cout << " % predictHelixOligomer " << endl;
	cout << "   --tmSequence <1-letter sequence> " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --backboneCrd <backboneCrd> --logFile <logFile> --helixGeoLine <geometryLine> hisProtState <HSD|HSE|HSP>" << endl;
	cout << "   --topFile <file> --parFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --uniprotName <name> --uniportAccession <name> --printTermEnergies <true/false>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << endl;
}

