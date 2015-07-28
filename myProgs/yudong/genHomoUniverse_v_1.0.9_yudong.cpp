
#include <iostream>
#include <stdio.h>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "SysEnv.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "HelixGenerator.h"
#include "PDBReader.h"
#include "AtomBondBuilder.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

string programName = "genHomoUniverse";
string programDescription = "This program creates helix oligomers of specified geometries and measures interhelical hydrogen bond energy. Finds the Dmin and also the Douts for each hydrogen bond at Dmin.";
string programAuthor = "Sabareesh Subramaniam, Samson Condon, and Samantha Andderson";
string programVersion = "1.0.7";
string programDate = "20 July 2015";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

struct Options {
	string pdbFile;
	string sequence;
	int startResNum;
	int numHelices;
	bool antiparallel;
	string helicalAxisPdbFile;
	string output;
	string outputDir;

	double xShiftStart;
	double zShiftStart;
	double axialRotStart;
	double crossingAngleStart;
	//vector<int> rotCount;	
	double xShiftEnd;
	double zShiftEnd;
	double axialRotEnd;
	double crossingAngleEnd;
	//vector<int> rotCount;	
	double xShiftSteps;
	double zShiftSteps;
	double axialRotSteps;
	double crossingAngleSteps;
	double axZShift; // to move the z along with the axialRotation along the trapezoid
	double zAxShift; // to move the axialRotation along with the zShift along the trapezoid
	//vector<int> rotCount;	
	bool unitCellCoord; // true if using unit cell coordinates w' and Z' instead of Cartesian coordinates w and Z
	string topFile;
	string parFile;
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

void printOptions(Options& _opt) {
	cout << "Program            " << programName << " v." << programVersion << "," << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl; 
	cout << "pdbFile            " <<  _opt.pdbFile << endl;
	cout << "sequence           " <<  _opt.sequence << endl;
	cout << "helicalAxisPdbFile " <<  _opt.helicalAxisPdbFile << endl;
	cout << "output             " <<  _opt.output << endl;
	cout << "topFile            " <<  _opt.topFile << endl;
	cout << "parFile            " <<  _opt.parFile << endl;
	cout << "hBondFile          " <<  _opt.hBondFile << endl;
	cout << "numHelices         " <<  _opt.numHelices << endl;
	cout << "antiparallel       " <<  _opt.antiparallel << endl;
	cout << "xShiftStart        " <<  _opt.xShiftStart << endl;
	cout << "zShiftStart        " <<  _opt.zShiftStart << endl;
	cout << "axialRotStart      " <<  _opt.axialRotStart << endl;
	cout << "crossingAngleStart " <<  _opt.crossingAngleStart << endl;
	cout << "xShiftEnd          " <<  _opt.xShiftEnd << endl;
	cout << "zShiftEnd          " <<  _opt.zShiftEnd << endl;
	cout << "axialRotEnd        " <<  _opt.axialRotEnd << endl;
	cout << "crossingAngleEnd   " <<  _opt.crossingAngleEnd << endl;
	cout << "xShiftSteps        " <<  _opt.xShiftSteps << endl;
	cout << "zShiftSteps        " <<  _opt.zShiftSteps << endl;
	cout << "axialRotSteps      " <<  _opt.axialRotSteps << endl;
	cout << "crossingAngleSteps " <<  _opt.crossingAngleSteps << endl;
	cout << "axZShift           " <<  _opt.axZShift << endl; 
	cout << "zAxShift           " <<  _opt.zAxShift << endl;
	cout << "unitCellCoord      " <<  _opt.unitCellCoord << endl;
}


Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);

//void transformCoiledCoil(AtomPointerVector& _chainA, AtomPointerVector& _chainB,double _zShift,double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, AtomPointerVector& _axisA, AtomPointerVector& _axisB);
void transformCoiledCoil(System &_sys, double _zShift, double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, bool _antiparallel);

string convertToPolymerSequence(string _seq, int _startResNum);

unsigned int getNumberOfInterHelicalHbonds(vector<Interaction*>& hbondInteractions);

map<string,double> getInterHelicalHbondInfo(vector<Interaction*>& hbondInteractions);

string printInfo(double axialRotate,double crossingAngle,double zShift,double relAxRot,double relZShift,double xShift);

bool createHelixdeNovo (System &_sys, CharmmSystemBuilder &_csb, HelixGenerator &_hg, PolymerSequence &_ps, string _sequence, int _startResNum);
bool createHelixFromPDB (System &_sys, PDBReader _pdbReader, string _pdbFile, CharmmSystemBuilder &_csb, PolymerSequence &_ps, string _sequence, int _startResNum);
bool duplicateHelix (System &_sys, CharmmSystemBuilder &_csb, int _helices);

bool orientHelix (System &_sys, Transforms &_tr, int _startResNum);
void applyRotationalSymmetry(System &_sys, Transforms &_tr, CartesianPoint _axis);
void applyAntiparallelSymmetry(System &_sys, Transforms &_tr, CartesianPoint _axis);
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB);
/*========================================== BEGIN MAIN ==========================================*/
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

	printOptions(opt);

	/******************************************************************************
	 *                     === SYSTEM SETUP ===
	 ******************************************************************************/
	
	// Declare System
	System sys;
	HelixGenerator hg("/exports/home/scondon/mslib/trunk/tables/PeptidesBBQTable.txt");
	PDBReader pdbRead;
	//hg.setHelixParameters(3.78, 100.0/180.0*M_PI, 51.2/180.0*M_PI);
	PolymerSequence PS;

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Set up Charmm System Builder
	CharmmSystemBuilder CSB(sys, opt.topFile, opt.parFile);
	CSB.setBuildNoTerms();
	CSB.setBuildTerm("CHARMM_VDW");
	CSB.setBuildTerm("SCWRL4_HBOND");

	if (opt.pdbFile!="") {
		cout << "Creating model helix from PDB File " << opt.pdbFile <<endl;
		createHelixFromPDB(sys, pdbRead, opt.pdbFile, CSB, PS, opt.sequence, opt.startResNum);
	} else {
		cout << "Creating model helix de novo" << endl;
		createHelixdeNovo(sys, CSB, hg, PS, opt.sequence, opt.startResNum);
	}

	//orientHelix(sys, trans, opt.startResNum);
	duplicateHelix(sys, CSB, opt.numHelices);
	sys.buildAllAtoms();
	CSB.updateNonBonded();

	
//	// Read in PDB File

//	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
//		cerr << "Unable to build system" << endl;
//		exit(0);
//	}


	// Build system
	//sys.buildAllAtoms();
//	cout << sys << endl;

	// Add hydrogen bonds
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	//HydrogenBondBuilder hb(sys, "/data00/bkmueller/dataFiles/hbondlist_nonCanon_adj.txt");
	hb.buildInteractions(30);
	cout << sys << endl;
	sys.calcEnergy();
	sys.printEnergySummary();
	sys.writePdb("initialState.pdb");
	cout << "==================***************************======================================" << endl;

	// Redirect Output
	string filename = opt.outputDir +"/"+ opt.output;
	if(!freopen (filename.c_str(), "w", stdout)) {
		cerr << "ERROR_freopen: Did not open output file " << opt.output << endl;
		exit(0);
	}
	// Set up APVs
	//AtomPointerVector &chainA = sys.getChain("A").getAtomPointers();
	//AtomPointerVector &chainB = sys.getChain("B").getAtomPointers();



	/******************************************************************************
	 *                     === MONTE CARLO SET UP ===
	 ******************************************************************************/

	System helicalAxis;
	helicalAxis.readPdb(opt.helicalAxisPdbFile);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/

	AtomSelection sel(sys.getAtomPointers());
	sel.select("chainA,chain A");
	sel.select("chainB,chain B");
	
	/******************************************************************************
	 *                     === INITIAL STARTING POSITION ===
	 ******************************************************************************/
	
	sys.saveCoor("initialState");
	helicalAxis.saveCoor("initialState");

	EnergySet* pESet = sys.getEnergySet();
	vector<Interaction*> hbondInteractions = (*(pESet->getEnergyTerms()))["SCWRL4_HBOND"];


	sys.saveEnergySubset("interHelical","chainA", "chainB");
	// Reference points to set up Helical starting postions
	CartesianPoint origin(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);
	CartesianPoint yAxis(0.0,1.0,0.0);


	vector<double> xShifts;
	for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
		xShifts.push_back(xShift);
	}

	int counter = 0;

	//Axial Rotation Loop
	for(double axialRotate = opt.axialRotStart; axialRotate < opt.axialRotEnd; axialRotate += opt.axialRotSteps) {
		//sys.applySavedCoor("initialState");
	//	AtomPointerVector chainA = sys.getChain("A").getAtomPointers();

		double zShiftStart = opt.zShiftStart + axialRotate * opt.axZShift; // opt.axZShift = -0.015, converts relative Z' to absolute Z
		double zShiftEnd = opt.zShiftEnd + axialRotate * opt.axZShift; // opt.axZShift = -0.015
		double relZ = opt.zShiftStart;

		//Z Shift Loop
		for( double zShift = zShiftStart; zShift < zShiftEnd; zShift += opt.zShiftSteps ,relZ += opt.zShiftSteps) {
			double shiftedAxialRotate = axialRotate + (zShift - axialRotate * opt.axZShift)  * opt.zAxShift; //convert relative w' to absolute w
		sys.applySavedCoor("initialState");
		AtomPointerVector chainA = sys.getChain("A").getAtomPointers();

			if (opt.antiparallel){
				applyAntiparallelSymmetry(sys, trans, yAxis);
			}
			AtomPointerVector chainB = sys.getChain("B").getAtomPointers();				

			//===== Axial Rotation ======
			trans.rotate(chainA, shiftedAxialRotate, origin, zAxis);
			if (opt.antiparallel){
				trans.rotate(chainB, -shiftedAxialRotate, origin, zAxis);
			}

			//====== Z Shift (Crossing Point) ======
			CartesianPoint zShiftCP(0.0, 0.0, zShift);
			trans.translate(chainA, zShiftCP);
			if (opt.antiparallel){
				CartesianPoint zShiftCPB (0.0, 0.0, -zShift);
				trans.translate(chainB, zShiftCPB);
			}
			sys.saveCoor("AxialZState");

			//Crossing Angle Loop
			for(double crossingAngle = opt.crossingAngleStart; crossingAngle < opt.crossingAngleEnd; crossingAngle += opt.crossingAngleSteps ) {
				vector<map<string,double> >  hBonds;
				vector<double> energies;	
				
				sys.applySavedCoor("AxialZState");
				//====== Local Crossing Angle ======
				trans.rotate(chainA, (crossingAngle/2.0), origin, xAxis);
				if (opt.antiparallel){
					trans.rotate(chainB, (-crossingAngle/2.0), origin, xAxis);
				}

				sys.saveCoor("crossingAngleState");

				//X Shift Loop
				double bestx = 0;
				double bestE = std::numeric_limits<double>::max();
				for(double xShift = opt.xShiftStart; xShift < opt.xShiftEnd; xShift += opt.xShiftSteps) {
					sys.applySavedCoor("crossingAngleState");

					//====== X shift (Interhelical Distance) =======
					CartesianPoint interDistVect;
					interDistVect.setCoor(xShift * -1.0/2.0, 0.0, 0.0);
					trans.translate(chainA, interDistVect);
				
					//====== Rotational Symmetry Z Axis=======
					if (!opt.antiparallel){
						applyRotationalSymmetry(sys, trans, zAxis);
					} else {
						CartesianPoint oppInterDistVect;
						oppInterDistVect.setCoor(xShift * 1/2, 0.0, 0.0);
						trans.translate(chainB, oppInterDistVect);
					}

					//double thisEnergy = sys.calcEnergy();
//					thisEnergy = sys.calcEnergyOfSubset("interHelical");
//					*****CODE FOR TESTING THE ENERGY SUBSETS*****
//					//map<string, map<string, vector<Interaction*> > >::iterator interHelicalSubset = pESet->getEnergyTermsSubsets()->find("interHelical");
//					for (map<string, vector<Interaction*> >::iterator k=interHelicalSubset->second.begin(); k!=interHelicalSubset->second.end(); k++) {
//						// for all the terms
//						for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
//							// for all the interactions
//								cout << (*l)->toString() << endl;
//						}
//					}
//					exit(0);
//

					//if (thisEnergy < bestE) {
					//	bestx = xShift;
					//	bestE = thisEnergy;
					//	sys.saveCoor("BestX");
					//}
					
					//unsigned int numHbonds = getNumberOfInterHelicalHbonds(hbondInteractions);
					//hBonds.push_back(getInterHelicalHbondInfo(hbondInteractions));
					//energies.push_back(thisEnergy);
				

			        	counter++;
			        	//if (counter ==1 || counter < 30) {
				        string name = opt.outputDir + "/" + printInfo(axialRotate,crossingAngle,relZ,shiftedAxialRotate,zShift,xShift) + ".pdb";
				        //sys.applySavedCoor("BestX");
			        	//transformCoiledCoil(sys,zShift,crossingAngle,shiftedAxialRotate,bestx,trans,origin,zAxis,xAxis,opt.antiparallel);
			       		sys.writePdb(name);
			        	//} else {
			        	//	string name = printInfo(axialRotate,crossingAngle,relZ,shiftedAxialRotate,zShift,xShifts,energies,hBonds);
			        	//}
			        }
                        }
		}
	}



	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;

}
/*========================================== END MAIN ==========================================*/

void applyRotationalSymmetry(System &_sys, Transforms &_tr, CartesianPoint _axis=CartesianPoint(0.0, 0.0, 1.0)) {
	double angle = 0.0;
	//number of chains to apply symmetry: 2 chains for C2, 3 for C3, etc
	unsigned int N = _sys.chainSize();
	AtomPointerVector chainA = _sys.getChain(0).getAtomPointers();
	for (unsigned int i = 1; i < N; i++) {
		angle += 360.0/N;
		AtomPointerVector nextChain = _sys.getChain(i).getAtomPointers();
		//Set coordinates of the next chain in the system to the first chain
		for (uint j = 0; j < chainA.size(); j++) {
			nextChain[j]->setCoor(chainA[j]->getCoor());
		}
		Matrix rotMat = CartesianGeometry::getRotationMatrix(angle, _axis);
		_tr.rotate(nextChain, rotMat);

	}

}

void applyAntiparallelSymmetry(System &_sys, Transforms &_tr, CartesianPoint _axis=CartesianPoint(0.0, 1.0, 0.0)) {
	//number of chains to apply symmetry: 2 chains for C2, 3 for C3, etc
	double xAngle = 180.0;
	unsigned int N = _sys.chainSize();
	//Even number of monomers only
	if (N%2 == 0){
		AtomPointerVector chainA = _sys.getChain(0).getAtomPointers();
		for (unsigned int i = 1; i < N; i+=2) {
			AtomPointerVector nextChain = _sys.getChain(i).getAtomPointers();
			//Set coordinates of the next chain in the system to the first chain
			for (uint j = 0; j < chainA.size(); j++) {
				nextChain[j]->setCoor(chainA[j]->getCoor());
			}
			Matrix rotMat = CartesianGeometry::getRotationMatrix(xAngle, _axis);
			_tr.rotate(nextChain, rotMat);
		}
	} 
}


void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->setCoor(_apvA[i]->getCoor());	
	}
	// Rotation matrix for 180 degrees
	//      ____________
	//     |            |
	//     | -1   0   0 |
	//     |            |
	//     |  0  -1   0 |
	//     |            |
	//     |  0   0   1 |
	//     |____________|
	
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

void transformCoiledCoil(System &_sys, double _zShift,double _crossingAngle, double _axialRotate,double _xShift,Transforms& _trans,CartesianPoint& _origin, CartesianPoint& _zAxis, CartesianPoint& _xAxis, bool _antiparallel) {
	AtomPointerVector chainA = _sys.getChain("A").getAtomPointers();
	//===== Apply Antiparallel Symmetry ======
	if (_antiparallel){
		CartesianPoint yAxis (0.0,1.0,0.0);
		applyAntiparallelSymmetry (_sys, _trans, yAxis);
		AtomPointerVector chainB = _sys.getChain("B").getAtomPointers();
		CartesianPoint zShiftCPB (0.0, 0.0, -_zShift);
		_trans.translate(chainB, zShiftCPB);
		_trans.rotate(chainB, -_axialRotate, _origin, _zAxis);
		_trans.rotate(chainB, (-_crossingAngle/2.0), _origin, _xAxis);
		CartesianPoint oppInterDistVect;
		oppInterDistVect.setCoor(_xShift * 1/2, 0.0, 0.0);
		_trans.translate(chainB, oppInterDistVect);
	}


	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(chainA, _axialRotate, _origin, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(chainA, (_crossingAngle/2.0), _origin, _xAxis);
	//_trans.rotate(_axisA,  (_crossingAngle/2.0), _origin, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor(_xShift * -1.0/2.0, 0.0, 0.0);
	_trans.translate(chainA, interDistVect);
	//_trans.translate(_axisA, interDistVect);
	
	//Apply Symmetry
	if (!_antiparallel){
		applyRotationalSymmetry(_sys, _trans, _zAxis);
		//c2Symmetry(_chainA, _chainB);
		//c2Symmetry(_axisA, _axisB);
	} 
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
		//cout << hbondInteractions[i]->toString() << endl;
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



string printInfo(double axialRotate,double crossingAngle,double zShift,double relAxRot,double relZShift,double xShift) {
	char name[100];
	sprintf(name,"model_%03.3f_" "%03.1f_" "%03.3f_" "%03.3f_" "%03.3f_" "%03.1f", axialRotate, crossingAngle, zShift, relAxRot, relZShift, xShift);
	return string(name);

}



//function to create idealized alpha helices of defined sequence
bool createHelixdeNovo (System &_sys, CharmmSystemBuilder &_csb, HelixGenerator &_hg, PolymerSequence &_ps, string _sequence, int _startResNum) {
	string polymerSequence = convertToPolymerSequence (_sequence, _startResNum);
	_ps.setSequence(polymerSequence);
	_csb.buildSystem(_ps);

	AtomPointerVector helix;
	_hg.generateHelix(helix, _sequence.length()+4); //extra turn added to prevent helix "unwinding" when building from PS
	_sys.assignCoordinates(helix, false);

	_sys.buildAllAtoms();
	return 1;
}

bool createHelixFromPDB (System &_sys, PDBReader _pdbReader, string _pdbFile, CharmmSystemBuilder &_csb, PolymerSequence &_ps, string _sequence, int _startResNum) {
	string polymerSequence = convertToPolymerSequence (_sequence, _startResNum);
	_ps.setSequence(polymerSequence);
	_csb.buildSystem(_ps);

	_pdbReader.open(_pdbFile);
	if (!_pdbReader.read()) {
		cout << "PDB READER ERROR: CANNOT READ " << _pdbFile << endl;
		exit(0);
	}
	_pdbReader.close();
	_sys.assignCoordinates(_pdbReader.getAtomPointers(), false);
	_sys.buildAllAtoms();
	return 1;

}
//function to orient a helix on the Z axis and with the middle alpha carbon on the X axis
bool orientHelix (System &_sys, Transforms &_tr, int _startResNum) {
	AtomPointerVector helix = _sys.getChain("A").getAtomPointers();
	int chainMidpoint = _sys.getChain("A").positionSize() / 2 + 1 + _startResNum;	//The alpha carbon at the chain midpoint will be set to lie on the X-axis
	string midpointCAID = "A " + MslTools::intToString(chainMidpoint) + " CA";
	CartesianPoint O(0.0,0.0,0.0);
	CartesianPoint X(1.0,0.0,0.0);
	CartesianPoint Y(0.0,1.0,0.0);
	CartesianPoint Z(0.0,0.0,1.0);
	
	CartesianPoint geoCenter = helix.getGeometricCenter();
	//Orient the helix on the Z-axis
	_tr.orient(_sys.getAtomPointers(), geoCenter, X, O, Z);
	_tr.translate(_sys.getAtomPointers(), O - geoCenter);
	
	//Place the middle alpha carbon on the X-axis
	CartesianPoint midCA = _sys.getAtom(midpointCAID).getCoor();
	_tr.orient(_sys.getAtomPointers(), midCA, X, O, Z);
	CartesianPoint distance = CartesianPoint(midCA.getX(), midCA.getY(), 0.0) - midCA;
	_tr.translate (_sys.getAtomPointers(), distance);
	return 1;
}

bool duplicateHelix (System &_sys, CharmmSystemBuilder &_csb, int _helices=1) {
	AtomBondBuilder Abb;
	for (unsigned int i = 1; i < _helices; i++) {
		_sys.duplicateChain("A");
		Abb.buildConnections(_sys.getChain(i).getAtomPointers());
	}
	return 1;
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


	opt.allowed.push_back("pdbFile");
	opt.required.push_back("sequence");
	opt.allowed.push_back("numHelices");
	opt.allowed.push_back("antiparallel");
	opt.allowed.push_back("startResNum");

	opt.required.push_back("output");
	opt.allowed.push_back("outputDir");

	opt.required.push_back("xShiftStart");
	opt.required.push_back("xShiftEnd");
	opt.required.push_back("xShiftSteps");

	opt.required.push_back("zShiftStart");
	opt.required.push_back("zShiftEnd");
	opt.required.push_back("zShiftSteps");

	opt.required.push_back("axialRotStart");
	opt.required.push_back("axialRotEnd");
	opt.required.push_back("axialRotSteps");

	opt.required.push_back("crossingAngleStart");
	opt.required.push_back("crossingAngleEnd");
	opt.required.push_back("crossingAngleSteps");
	opt.required.push_back("axZShift");
	opt.required.push_back("zAxShift");

	//opt.required.push_back("rotCount");
	//
	opt.required.push_back("helicalAxisPdbFile");
	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("hBondFile");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("unitCellCoord");

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
		opt.warningMessages = "pdb file not specified";
		opt.warningFlag = true;
	}
	opt.sequence = OP.getString("sequence");
	if (OP.fail()) {
		opt.errorMessages = "polymer sequence not specified";
		opt.errorFlag = true;
	}
	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.startResNum = 1;
		opt.warningMessages += "starting residue number not specified; defaulting to 1";
		opt.warningFlag = true;
	}
	opt.numHelices = OP.getInt("numHelices");
	if (OP.fail()) {
		opt.numHelices = 2;
		opt.warningMessages += " number of helices not specified, defaulting to 2";
		opt.warningFlag = true;
	}
	opt.antiparallel = OP.getBool("antiparallel");
	if (OP.fail()){
		opt.antiparallel = false;
		opt.warningMessages += "antiparall not indicated, defaulting to parallel";
		opt.warningFlag = true;
	} else {
		if (opt.numHelices%2==1){
			opt.errorMessages = "cannot be antiparallel with an odd number of helices";
			opt.errorFlag = true;
		}
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
	opt.outputDir = OP.getString("outputDir");
	if (OP.fail()) {
		opt.outputDir="./";
		opt.warningMessages = "output directory not specified";
		opt.warningFlag = true;
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

	opt.zShiftStart = OP.getDouble("zShiftStart");
	if (OP.fail()) {
		opt.errorMessages = "zShiftStart not specified";
		opt.errorFlag = true;
	}
	opt.zShiftEnd = OP.getDouble("zShiftEnd");
	if (OP.fail()) {
		opt.errorMessages = "zShiftEnd not specified";
		opt.errorFlag = true;
	}
	opt.zShiftSteps = OP.getDouble("zShiftSteps");
	if (OP.fail()) {
		opt.errorMessages = "zShiftSteps not specified";
		opt.errorFlag = true;
	}

	opt.axialRotStart = OP.getDouble("axialRotStart");
	if (OP.fail()) {
		opt.errorMessages = "axialRotStart not specified";
		opt.errorFlag = true;
	}
	opt.axialRotEnd = OP.getDouble("axialRotEnd");
	if (OP.fail()) {
		opt.errorMessages = "axialRotEnd not specified";
		opt.errorFlag = true;
	}
	opt.axialRotSteps = OP.getDouble("axialRotSteps");
	if (OP.fail()) {
		opt.errorMessages = "axialRotSteps not specified";
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
	opt.axZShift = OP.getDouble("axZShift");
	if (OP.fail()) {
		opt.errorMessages = "axZShift not specified";
		opt.errorFlag = true;
	}

	opt.zAxShift = OP.getDouble("zAxShift");
	if (OP.fail()) {
		opt.errorMessages = "zAxShift not specified";
		opt.errorFlag = true;
	}

	opt.crossingAngleSteps = OP.getDouble("crossingAngleSteps");
	if (OP.fail()) {
		opt.errorMessages = "crossingAngleSteps not specified";
		opt.errorFlag = true;
	}

	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(envVar)) {
			opt.topFile = SYSENV.getEnv(envVar);
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
		if(SYSENV.isDefined(envVar)) {
			opt.parFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(SYSENV.isDefined(envVar)) {
			opt.hBondFile = SYSENV.getEnv(envVar);
			opt.warningMessages += "hBondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hBondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.unitCellCoord = OP.getBool("unitCellCoord");
	if (OP.fail()) {
		opt.unitCellCoord = false;
		opt.warningMessages += "Unit Cell coordinates not specified using --unitCellCoord. Assuming Cartesian coordinates for axRot and Z\n";
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
	cout << " % genHomoUniverse --pdbFile <pdbfile> --helicalAxisPdbFile <filename> --output <filename>" << endl;
	cout << " --xShiftStart <double> --xShiftEnd <double> --xShiftSteps <double> " << endl;
	cout << " --zShiftStart <double> --zShiftEnd <double> --zShiftSteps <double> " << endl;
	cout << " --axialRotStart <double> --axialRotEnd <double> --axialRotSteps <double> " << endl;
	cout << " --crossingAngleStart <double> --crossingAngleEnd <double> --crossingAngleSteps <double> " << endl;
	cout << endl;
}


string convertToPolymerSequence(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
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
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps;
}
