#include "repackSideChains.h"

void loadRotamers(System& _sys, Options& _opt) {
//	cout << "Include WT " << _opt.includeCR << endl;
	SystemRotamerLoader sysRot;
	sysRot.setSystem(_sys);
	sysRot.defineRotamerSamplingLevels();

	if(!sysRot.readRotamerLibraryFile(_opt.rotlibFile)) {
		cerr << "Unable to read " << _opt.rotlibFile << endl;
		exit(0);
	}

	vector<Position*> & positions = _sys.getPositions();

	for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++) {
		Position &pos = **p;
		string resName = pos.getResidueName();

		if(_opt.fixedPosMap.find(pos.getPositionId()) != _opt.fixedPosMap.end()) {
			cout << "Not repacking " << pos.getPositionId() << endl;
			continue;
		}

	//	cout << "Loading " << _opt.numRots[resName] << " at " << pos.getPositionId() << endl;
	//	cout << pos.getTotalNumberOfRotamers() << endl;
		
		if(_opt.numRots.find(resName) != _opt.numRots.end()) {
			if(_opt.numRots[resName] < 1) {
				continue;
			}
			if (!sysRot.loadRotamers(&pos, resName,_opt.numRots[resName],"",_opt.includeCR)) {
				cerr << "Cannot load rotamers " << pos.getPositionId() << "," << resName << endl;
				exit(0);
			} 
		} else if (_opt.rotLevel != "") {
			if (!sysRot.loadRotamers(&pos, resName,_opt.rotLevel,"",_opt.includeCR)) {
				cerr << "Cannot load rotamers " << pos.getPositionId() << "," << resName << endl;
				exit(0);
			}
		} else {
			// retain the crystal structure
			cerr << "Not repacking " << pos.getPositionId() << endl; 
			continue;
		}
	}

}

void rotlibRepack (System& _sys, Options& _opt) {

	SidechainOptimizationManager scom;
	scom.setSystem(&_sys);
	if(_opt.useTimeToSeed) {
		scom.seed(time(NULL));
	} else {
		scom.seed(_opt.seed);
	}
	scom.setOnTheFly(_opt.onTheFly);

	scom.calculateEnergies();
	scom.setVerbose(_opt.verbose);

	if(!_opt.runGreedy) {
		scom.setRunDEE(_opt.runGoldsteinSingles,_opt.runGoldsteinPairs);
		scom.setEnumerationLimit(10000);
		scom.setMCOptions(_opt.startT,_opt.endT,_opt.nCycles,_opt.shape,_opt.maxReject,_opt.deltaSteps,_opt.minDeltaE);

		scom.setRunSCMF(_opt.runSCMF);
		scom.setRunSCMFBiasedMC(_opt.runSCMFBiasedMC);
		scom.setRunUnbiasedMC(_opt.runUnbiasedMC);
		scom.runOptimizer();
	} else {
		scom.runGreedyOptimizer(_opt.greedyCycles);
	}


	vector<unsigned int> best = scom.getMinStates()[0];
	_sys.setActiveRotamers(best);
	
}



/******************************************
 *  
 *  =======  MAIN METHOD  =======
 *
 ******************************************/

int main(int argc, char *argv[]) {
	Options opt = parseOptions(argc, argv);
	if(opt.errorFlag) {
		cout << opt.OPerrors ;
		cout << opt.errorMessages ;
		exit(0);
	}

	if(opt.warningFlag) {
		cout << opt.warningMessages;
	}
	
	cout << "PDB: " << opt.pdbFile << endl;

	if(opt.rotLevel == "") {
		cout << "Number of Rotamers" << endl;
		for(map<string,int>::iterator it = opt.numRots.begin(); it != opt.numRots.end(); it++) {
			cout << it->first << "\t" << it->second << endl;
		}
	}
	
	time_t start,end;
	time(&start);

	System sys;
	CharmmSystemBuilder CSB(sys,opt.charmmTopFile,opt.charmmParFile);
	CSB.setBuildNonBondedInteractions(false);
	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
		cerr << "Unable to build" << opt.pdbFile << endl;
		exit(0);
	}

	// Add Side Chains
	sys.buildAllAtoms();

	CSB.updateNonBonded(opt.cuton,opt.cutoff,opt.cutnb);

	HydrogenBondBuilder HBB(sys,opt.hbondParFile);
	HBB.buildInteractions(opt.cuthb); // 

	EnergySet *eSet = sys.getEnergySet();
	for(int i = 0; i < opt.excludeTerms.size(); i++) {
		eSet->setTermActive(opt.excludeTerms[i],false);
	}


	cout << "Crystal Energy " << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();

	loadRotamers(sys,opt);
	rotlibRepack(sys,opt);

	time(&end);
	cout << "Repack took " << difftime(end,start) << " seconds" << endl;

	cout << "Repack Energy " << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();

	if(opt.outputPDBFile != "" ) {
		if(!sys.writePdb(opt.outputPDBFile)) {
			cerr << "Unable to write repacked structure to " << opt.outputPDBFile << endl;
		}
	}
}	




