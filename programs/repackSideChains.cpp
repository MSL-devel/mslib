#include "repackSideChains.h"

vector<double> chi1s;
vector<double> chi2s;
vector<double> altchi2s;
vector<double> sasa;

bool hasChi2(string resName) {
	if(resName == "ALA" || resName == "CYS" || resName == "GLY" || resName == "PRO" || resName == "SER" || resName == "THR" || resName == "VAL") {
		return false;
	} 
	return true;
}

bool hasAltChi2(string resName) {
	if(resName == "ASP" || resName == "PHE"  || resName == "TYR") {
		return true;
	} else {
		return false;
	}
}


double getChi1 (Residue& _res) {
	string name = _res.getResidueName();
	vector<string> atomNames;
	atomNames.push_back("N");
	atomNames.push_back("CA");
	atomNames.push_back("CB");
	if(name == "ILE" || name == "VAL") {
		atomNames.push_back("CG1");
	} else if(name == "CYS") {
		atomNames.push_back("SG");
	} else if(name == "SER") {
		atomNames.push_back("OG");
	} else if(name == "THR") {
		atomNames.push_back("OG1");
	} else {
		atomNames.push_back("CG");
	}

	vector<Atom*> atoms;
	for(int i = 0; i < atomNames.size(); i++) {
		if(_res.atomExists(atomNames[i])) {
			atoms.push_back(&_res.getLastFoundAtom());
		} else {
			return 1000.0;
		}
	}	

	return (atoms[0]->dihedral(*atoms[1],*atoms[2],*atoms[3]));
}

double getChi2 (Residue & _res) {
	string name = _res.getResidueName();
	if(!hasChi2(name)) {
		return 1000.0;
	}
	vector<string> atomNames;
	atomNames.push_back("CA");
	atomNames.push_back("CB");
	if(name == "ILE" ) {
		atomNames.push_back("CG1");
	} else {
		atomNames.push_back("CG");
	}
	if(name == "LEU" || name == "PHE" || name == "TRP" || name == "TYR") {
		atomNames.push_back("CD1");
	} else if(name == "ASN" || name == "ASP") {
		atomNames.push_back("OD1");
	} else if(name == "HSD" || name == "HSE" || name == "HSP") {
		atomNames.push_back("ND1");
	} else if(name == "MET") {
		atomNames.push_back("SD");
	} else {
		atomNames.push_back("CD");
	}

	vector<Atom*> atoms;
	for(int i = 0; i < atomNames.size(); i++) {
		if(_res.atomExists(atomNames[i])) {
			atoms.push_back(&_res.getLastFoundAtom());

		} else {
			return 1000.0;
		}
	}	
	return(atoms[0]->dihedral(*atoms[1],*atoms[2],*atoms[3]));
}

double getAltChi2(Residue& res) {
	
	if(!hasAltChi2(res.getResidueName())) {
		return 1000.0;
	}

	Atom* CA = NULL;
	Atom* CB = NULL;
	Atom* CG = NULL;
	Atom* atom4 = NULL;
	if(!res.atomExists("CA")) {
		return 1000.0;		
	} else {
		CA = &(res.getLastFoundAtom());
	}
	if(!res.atomExists("CB")) {
		return 1000.0;		
	} else {
		CB = &(res.getLastFoundAtom());
	}
	if(!res.atomExists("CG")) {
		return 1000.0;		
	} else {
		CG = &(res.getLastFoundAtom());
	}

	if(res.atomExists("CD2") || (res.atomExists("OD2") && res.getResidueName() == "ASP")) {
		atom4 = &(res.getLastFoundAtom());
	} else {
		return 1000.0;
	}
	return CA->dihedral(*CB,*CG,*atom4);
}

void loadRotamers(System& _sys, Options& _opt) {
//	cout << "Include WT " << _opt.includeCR << endl;
	SystemRotamerLoader sysRot;
	sysRot.setSystem(_sys);

	if(!sysRot.readRotamerLibraryFile(_opt.rotlibFile)) {
		cerr << "Unable to read " << _opt.rotlibFile << endl;
		exit(0);
	}

	vector<unsigned> varPos;
	vector<Position*> & positions = _sys.getPositions();

	for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++) {
		Position &pos = **p;
		Residue & res = pos.getCurrentIdentity();
		string resName = pos.getResidueName();

		if(_opt.numRots.find(resName) == _opt.numRots.end() || _opt.numRots[resName] < 1) {
			continue;
		}
	//	cout << "Loading " << _opt.numRots[resName] << " at " << pos.getPositionId() << endl;
	//	cout << pos.getTotalNumberOfRotamers() << endl;
		varPos.push_back(p-positions.begin());
		if(_opt.printStats) {
			chi1s.push_back(getChi1(res));
			chi2s.push_back(getChi2(res));
			altchi2s.push_back(getAltChi2(res));
			sasa.push_back(pos.getSasa());
		}
		if (!sysRot.loadRotamers(&pos, resName,_opt.numRots[resName],"",_opt.includeCR)) {
			cerr << "Cannot load rotamers " << pos.getPositionId() << "," << resName << endl;
			exit(0);
		} 
	}
	_sys.setVariablePositions(varPos);

}

void rotlibRepack (System& _sys, Options& _opt) {

	SidechainOptimizationManager socm;
	socm.setSystem(&_sys);
	if(_opt.useTimeToSeed) {
		socm.seed(time(NULL));
	} else {
		socm.seed(_opt.seed);
	}
	socm.setOnTheFly(_opt.onTheFly);

	socm.calculateEnergies();

	socm.setRunDEE(_opt.runGoldsteinSingles,_opt.runGoldsteinPairs);
	socm.setEnumerationLimit(10000);
	socm.setMCOptions(_opt.startT,_opt.endT,_opt.nCycles,_opt.shape,_opt.maxReject,_opt.deltaSteps,_opt.minDeltaE);

	socm.setRunSCMF(_opt.runSCMF);
	socm.setRunSCMFBiasedMC(_opt.runSCMFBiasedMC);
	socm.setRunUnbiasedMC(_opt.runUnbiasedMC);
	socm.setVerbose(_opt.verbose);

	socm.runOptimizer();

	vector<unsigned int> best = socm.getMinStates()[0];
	_sys.setActiveRotamers(best);
	
}


void printStats(System& _sys) {
	vector<unsigned>  varPos = _sys.getVariablePositions();
	cout << "Output Format : " << endl;
	cout << "\tResidueID\torigChi1\torigChi2\torigAltChi2\torigSasa\trecdChi1\trecdChi2\trecdSasa" << endl;
	for(int i = 0; i < varPos.size(); i++) {
		Residue& res = _sys.getPosition(varPos[i]).getCurrentIdentity();
		cout << "STAT " << res.getIdentityId() << " " << chi1s[i] << " " << chi2s[i] << " " << altchi2s[i] << " " << sasa[i] << " ";
						cout << getChi1(res) << " " << getChi2(res) << " " << getAltChi2(res) << " " << res.getSasa() << endl;
	}
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
	cout << "Number of Rotamers" << endl;

	for(map<string,int>::iterator it = opt.numRots.begin(); it != opt.numRots.end(); it++) {
		cout << it->first << "\t" << it->second << endl;
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

	SasaCalculator sc(sys.getAtomPointers());
	if(opt.printStats) {
		sc.calcSasa();
	}

	cout << "Crystal Energy " << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();

	loadRotamers(sys,opt);
	rotlibRepack(sys,opt);

	time(&end);
	cout << "Repack took " << difftime(end,start) << " seconds" << endl;

	cout << "Repack Energy " << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();

	if(opt.printStats) {
		sc.calcSasa();
		printStats(sys);
	}
	if(opt.outputPDBFile != "" ) {
		if(!sys.writePdb(opt.outputPDBFile)) {
			cerr << "Unable to write repacked structure to " << opt.outputPDBFile << endl;
		}
	}
}	




