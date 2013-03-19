#include <iostream>
#include <fstream>

#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "CRDReader.h"
#include "SysEnv.h"

using namespace MSL;
using namespace std;

string programName = "filterOligomerByConstraints";
string programDescription = "This program repacks oligomeric helices from a set starting position based on constraints";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "0.0.0";
string programDate = "19 March 2013";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

static SysEnv SYSENV;
void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {

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
string convertToPolymerSequence(string _seq, int _startResNum, string _protState) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " " + _protState;
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

void repackSideChains(System & _sys, SelfPairManager & _spm, bool _greedy, int _greedyCycles) {

	_spm.calculateEnergies();
	_spm.runGreedyOptimizer(_greedyCycles);
}
double computeMonomerEnergy(System & _sys, Transforms & _trans, RandomNumberGenerator & _RNG, SystemRotamerLoader & _sysRot, SelfPairManager & _spm) {

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector &chainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector &chainB = _sys.getChain("B").getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES 500 A APART ===
	 ******************************************************************************/
	double xShift = 500.0;

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*xShift/2.0), 0.0, 0.0);
	_trans.translate(chainA, interDistVect);

	c2Symmetry(chainA, chainB);

	// Repack side chains
	repackSideChains(_sys, _spm, true, 1);

	return _spm.getMinBound()[0];

}


int main(int argc, char *argv[]) {

	if(argc != 6) {
		cerr << "Usage: filterOligomerByConstraints <startThread> <endThread> <hisProtonation state> <hbondList> <sequence>" << endl;
		exit(0);
	}

	int start = atoi(argv[1]);
	int end = atoi(argv[2]);
	string protState = string(argv[3]);

	/********************************************************************************
	*
	*Let's start the helices at a 12A separation, and move them closer 0.2A at the time (0.1A if the program is very fast) until things clash.
	*
	*The output should be the grid points in which the His is hydrogen bonded either as NE2 donor or as OE1 acceptor, and for each grid point the dout (where the hydrogen bonding starts) and din (where things become clashy).
	*
	*
	******************************************************************************/
	

	time(&startTime);	

	string modelledTMSeq = string(argv[5]) ; // replace P with A
	
	/**********************************************************************************
	*
	*    printProteinOutFile
	*
	**********************************************************************************/


	//cout << modelledTMSeq << endl;

	string sequence = convertToPolymerSequence(modelledTMSeq,1,protState); 
	PolymerSequence PS(sequence); 

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,SYSENV.getEnv("MSL_CHARMM_TOP"),SYSENV.getEnv("MSL_CHARMM_PAR"), SYSENV.getEnv("MSL_CHARMM_SOLV"));
	CSB.setSolvent("CHEX");

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << sequence << endl;
		exit(0);
	} else {
		//cout << "Sequence Built " << sequence << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();


	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open("/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd");
	if(!cRead.read()) {
		cout << "Unable to read /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd" << endl;
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
	RNG1.setSeed(1);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, "/data01/sabs/msl_working/mslib/trunk/library/EBL_11-2011_CHARMM22.txt");
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, "/data01/sabs/msl_working/mslib/trunk/toppar/scwrl4hb/par_hbond_2.txt");
	hb.buildInteractions(-1.0);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();

	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsInactive();
	Eset->setTermActive("CHARMM_VDW", true);
	Eset->setTermActive("SCWRL4_HBOND", true);

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
	spm.setVerbose(false);

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
	double monomerEnergy = computeMonomerEnergy(sys,trans,RNG1,sysRot,spm);
	cout << "monomerEnergy: " << monomerEnergy << endl;

	vector<Interaction*> qbonds;
	
	vector<Interaction*> hbondInteractions = (*(Eset->getEnergyTerms()))["SCWRL4_HBOND"];

	map<string,vector<string> > hbonds;
	// read the hbondlist file and populate the map from donor to acceptor
	ifstream file;
	file.open(argv[4]);
	if(!file.is_open()) {
		cerr << "Unable to open " << file << endl;
		exit(0);
	}

	string tmp;
	while(file) {
		getline(file,tmp);
		tmp = MslTools::uncomment(tmp);
		if(tmp.length() > 4) { // there will be atleast 4 commas
			vector<string> toks = MslTools::tokenizeAndTrim(tmp);
			if(hbonds.find(toks[0]) != hbonds.end()) {
				hbonds[toks[0]].push_back(toks[1]);
			} else {
				hbonds[toks[0]] = vector<string>();
				hbonds[toks[0]].push_back(toks[1]);
			}
			//cout << "LOOK FOR : " << toks[0] << " " << toks[1] << endl;
		}

	}



//	cout << hbondInteractions.size() << endl;

	// look for interhelical HIS-HIS or HIS-TRP or HIS-THR sidechain hydrogen bonds
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms = hbondInteractions[i]->getAtomPointers();
		//cout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		// ok they are interhelical
		string donor = atoms[0]->getAtomId();
		string acceptor = atoms[2]->getAtomId();
		if(hbonds.find(donor) != hbonds.end()) {
			vector<string> list = hbonds[donor];
			for(int j = 0; j < list.size(); j++) {
				if(list[j] == acceptor) {
					qbonds.push_back(hbondInteractions[i]);
					//cout << "FOUND " << donor <<  " " << acceptor << endl;
					break;
				}
			}
		}
	}

	cout << qbonds.size() << endl;


	/******************************************************************************
	 *              === LOOP OVER ALL POSSIBLE INTERFACE POSITIONS ===
	 ******************************************************************************/
	for(int j = start ; j <= end ; j++) {

		//bool computedMonomerEnergy = false;
		chainA.renumberChain(j);
		chainB.renumberChain(j);

		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);		

		sys.buildAllAtoms();
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
		sys.saveAltCoor("currentTMFragment");

	
		/******************************************************************************
		 *              === LOOP OVER ALL grid points ===
		 ******************************************************************************/
		for(double zShift = 0; zShift < 6.01; zShift += 0.6) {
			for(double crossingAngle = -55; crossingAngle < 56; crossingAngle += 5.0) {
				for(double axialRotation = -100; axialRotation < 1; axialRotation += 20) {

					double prevEnergy = MslTools::floatMax;

					vector<double> dmins(qbonds.size(),MslTools::floatMax);
					vector<double> douts(qbonds.size(),0);

					//cout <<  " A:" << j << " B:" << j << " zShift: " << zShift << " crossingAngle: " << crossingAngle << " axialRot: " << axialRotation << endl;
					for(double xShift = 12.0; xShift >= 6.0; xShift -= 0.2) {

						helicalAxis.applySavedCoor("originState");
						sys.applySavedCoor("currentTMFragment");

						// Transform helices to initial starting position
						transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
						
						/*
						string filename = "/tmp/new/pdb_" + MslTools::doubleToString(xShift) + ".pdb";
						if(!sys.writePdb(filename)) {
							cerr << "Unable to write " << filename << endl;
							exit(0);
						}
						continue;
						*/
						// Optimizatize Initial Starting Position
						repackSideChains(sys, spm, true, 1);

						sys.setActiveRotamers(spm.getMinStates()[0]);
						double currentEnergy = spm.getMinBound()[0];
						//cout << " xShift: " << xShift << " " << currentEnergy << endl;
					
						for(int l = 0; l < qbonds.size(); l++) {
							Interaction* inter = qbonds[l];
							if(inter->getEnergy() < 0) {
								dmins[l] = xShift;
								if(xShift > douts[l]) {
									douts[l] = xShift;
								}
							}
						}
						if (currentEnergy > prevEnergy && currentEnergy > (monomerEnergy + 10.0)) {
							break;

						}
						prevEnergy = currentEnergy;
						
					} // xshift loop
					
					for(int l = 0; l < qbonds.size(); l++) {
						vector<Atom*> atoms = qbonds[l]->getAtomPointers();
						string name = atoms[0]->getAtomOfIdentityId() + ";" + atoms[2]->getAtomOfIdentityId();
						if(douts[l] > 0) {
							cout <<  "HBOND A:" << j << " B:" << j << " " << name;
							cout << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
							cout << " dmin: " << dmins[l] << " dout: " << douts[l] << endl;
						}
					}
				} // axial loop
			} // crossingAngle loop
		} // zShift loop
	} // Threading Loop End


	time(&endTime);
	diffTime = difftime (endTime, startTime);
	cout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
}

