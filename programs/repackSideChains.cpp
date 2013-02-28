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

	if (_opt.rotLevel != "") {
		unsigned int numAla = sysRot.getRotamerLibrary()->getLevel(_opt.rotLevel, "ALA");
		if (numAla == 0) {
			cout << _opt.rotLevel << " does not contain ALA, using wt rotamer alone" << endl;
			_opt.numRots["ALA"] = 0;
		}
		unsigned int numGly = sysRot.getRotamerLibrary()->getLevel(_opt.rotLevel, "GLY");
		if (numGly == 0) {
			cout << _opt.rotLevel << " does not contain GLY, using wt rotamer alone" << endl;
			_opt.numRots["GLY"] = 0;
		}
		unsigned int numPro = sysRot.getRotamerLibrary()->getLevel(_opt.rotLevel, "PRO");
		if (numPro == 0) {
			cout << _opt.rotLevel << " does not contain PRO, using wt rotamer alone" << endl;
			_opt.numRots["PRO"] = 0;
		}
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
	/*
	cout << "FixedE " << scom.getFixedEnergy() << endl;
		vector<vector<double> > & se = scom.getSelfEnergy();
		cout << "SelfE " << endl;
		for(int i = 0; i < se.size(); i++) {
			for(int j = 0; j < se[i].size(); j++) {
				cout << se[i][j] << endl;
			}
		}
		vector<vector<vector<vector<double> > > > & pe = scom.getPairEnergy();
		cout << "PairE " << endl;
		for(int i = 0; i < pe.size(); i++) {
			for(int j = 0; j < pe[i].size(); j++) {
				for(int k = 0; k < pe[i][j].size(); k++) {
					for(int l = 0; l < pe[i][j][k].size(); l++) {
						cout << pe[i][j][k][l] << endl;
					}
				}
			}

		}
		exit(0);
		*/

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
	if(opt.solvFile != "") {
		if(!CSB.readSolvation(opt.solvFile)) {
			cerr << "Unable to read solvation file " << opt.solvFile << endl;
			exit(0);
		}
		CSB.setSolvent(opt.solvent);
	}
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




