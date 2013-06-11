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
#include <string>
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "Timer.h"
#include "myEnergyOptimizations.h"
#include "MonteCarloManager.h"


using namespace std;

using namespace MSL;

// Global objects
Timer t;
double startTime = 0.0;

void printRots (SelfPairManager& scom,vector<Position*>& masterPositions) {
	vector<unsigned int> numRots = scom.getNumberOfRotamers();
	for(int i = 0; i < numRots.size(); i++) {
		cout << "NUMROTS " << masterPositions[i]->getResidueName() << ":" << numRots[i] << endl;
	}
}


int main(int argc, char *argv[]){

	// Option Parser
	cout << "Setup Options"<<endl;
	MonteCarloOptions opt = setupMonteCarloOptions(argc, argv);
	cout << "Create System"<<endl;
	// Create a system from the structural input options
	startTime = t.getWallTime();
	System sys;
	createSystem(opt.energyOpt.structOpt, sys,&opt.energyOpt);
	cout << "Built system after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	startTime = t.getWallTime();

	//EnergySet *pESet = sys.getEnergySet();

	SelfPairManager scom;
	map<string,double> sequenceEnergy;

	scom.setSystem(&sys);
	scom.setOnTheFly(true);
	//scom.saveInteractionCounts(true);
	scom.setVerbose(false);
	scom.calculateEnergies();
	scom.seed(opt.randomSeed);

	cout << "Calc Fixed/Self energies after "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	startTime = t.getWallTime();
	scom.runGreedyOptimizer(1);
	vector<unsigned int> rotamerState = (scom.getMinStates())[0];
	sys.setActiveRotamers(rotamerState);
	/*
	cout << "Initial AtomSize " << sys.atomSize() << endl; 
	if(!sys.writePdb("initial.pdb")) {
		cerr << "Unable to write initial.pdb" << endl;
		exit(0);
	}
	*/
	double e = sys.calcEnergy();
	string sequence = PolymerSequence::toOneLetterCode(sys.getChain(0).getAtomPointers());
	fprintf(stdout,"SEQ %s ENERGY %8.3f %8.3f\n",sequence.c_str(),e,(scom.getMinBound())[0]);
	//cout << "SPM Count " << scom.getStateInteractionCount(rotamerState) << " SYS Count " << pESet->getTotalNumberOfInteractionsCalculated() << endl;;
	if(sequenceEnergy.find(sequence) == sequenceEnergy.end() || sequenceEnergy[sequence] > e ) {
		sequenceEnergy[sequence] = e;
	}
	vector<unsigned int> masterIndices = sys.getMasterPositions();
	vector<Position*> masterPositions;
	vector<vector<string> > identities;
	identities.resize(masterIndices.size());
	for(int i = 0; i < masterIndices.size(); i++) {
		masterPositions.push_back(&(sys.getPosition(masterIndices[i])));
		string thisResName = masterPositions.back()->getResidueName();
		int numIds = masterPositions.back()->identitySize();
		identities[i].resize(numIds);
		for(int j = 0; j < numIds; j++) {
			identities[i][j] = (masterPositions.back()->getResidue(j)).getResidueName();
		}
		masterPositions.back()->hideAllIdentitiesButOne(thisResName);
//		vector<Position*> linked = masterPositions[i]->getLinkedPositions();
	//	cout << "AtomSize " << sys.atomSize() << " Done with " << masterPositions[i]->getPositionId() << endl;
	}

	/*
	for(int i = 0; i < identities.size(); i++) {
		cout << masterPositions[i]->getPositionId() << "  ";
		for(int j = 0; j < identities[i].size();j++) {
			cout << identities[i][j] << " ";
		}
		cout << endl;
	}
	*/
	RandomNumberGenerator rng;
	rng.setSeed(opt.randomSeed);

	MonteCarloManager MCMngr(10000,999,100,MonteCarloManager::CONSTANT,100,opt.deltaSteps,0);
	MCMngr.setRandomNumberGenerator(&rng);
	MCMngr.setEner(e);

	//printRots(scom,masterPositions);
	
	/*
	scom.setSystem(&sys);
	vector<Position*> positions = sys.getPositions();
	cout << "Listing all positions" << endl;
	for(int i = 0; i < positions.size(); i++) {
		cout << positions[i]->getPositionId() << "," << positions[i]->getResidueName() << endl;
	}
	for(int i = 0; i < masterPositions.size(); i++) {
		vector<Position*> linked = masterPositions[i]->getLinkedPositions();
		cout << "REAL " << masterPositions[i]->getPositionId() << "," << masterPositions[i]->getResidueName() <<  " " << linked[0]->getPositionId() << "," << linked[0]->getResidueName() << endl;
	}
	
	scom.calculateEnergies();
	scom.runGreedyOptimizer(1);
	scom.debug();
	double repackEnergy = scom.getMinBound()[0];
	rotamerState = (scom.getMinStates())[0];
	sys.setActiveRotamers(rotamerState);

	cout << "Energy after hiding stuff " << endl;
	if(!sys.writePdb("afterHiding.pdb")) {
		cerr << "Unable to write afterHiding.pdb" << endl;
		exit(0);
	}
	e = sys.calcEnergy();
	printRots(scom,masterPositions);
	fprintf(stdout,"SEQ %s ENERGY %8.3f %8.3f\n",sequence.c_str(),e,repackEnergy);
	cout << "SPM Count " << scom.getStateInteractionCount(rotamerState) << " SYS Count " << pESet->getTotalNumberOfInteractionsCalculated() << endl;;
	exit(0);
	*/

	int seqNum = 0;
	while(!MCMngr.getComplete()) {
		int posIndex = rng.getRandomInt(0,masterPositions.size()-1);
		Position* selPos = masterPositions[posIndex]; 
		string origIdentity = (selPos->getResidueName());
		int identityIndex = rng.getRandomInt(0,identities[posIndex].size()-1);
		string selIdentity = identities[posIndex][identityIndex];
		if(origIdentity == selIdentity) {
			cout << "Skip " << posIndex << " " << origIdentity << " " << selIdentity << " " << identityIndex << endl;
			continue;
		}
		cout << "Mutating " << selPos->getPositionId() << " " << origIdentity << " " << selIdentity << endl;
		if(!selPos->hideAllIdentitiesButOne(selIdentity)) {
			cerr << "Unable to hide-all-but " << selIdentity << " " << selPos->getPositionId() << endl;
			exit(0);
		}
		scom.setSystem(&sys);
		scom.calculateEnergies();
		scom.runGreedyOptimizer(1);
		double repackEnergy = scom.getMinBound()[0];
		vector<unsigned int> rotamerState = (scom.getMinStates())[0];
		sys.setActiveRotamers(rotamerState);
		cout << "STATE ";
		for(int i = 0; i < rotamerState.size(); i++) {
			cout << rotamerState[i] << ",";
		}
		cout << endl;

		string sequence = PolymerSequence::toOneLetterCode(sys.getChain(0).getAtomPointers());
		//printRots(scom,masterPositions);
		double e = sys.calcEnergy();
		//cout << "SPM Count " << scom.getStateInteractionCount(rotamerState) << " SYS Count " << pESet->getTotalNumberOfInteractionsCalculated() << endl;;
		if(sequenceEnergy.find(sequence) == sequenceEnergy.end() || sequenceEnergy[sequence] > e ) {
			sequenceEnergy[sequence] = e;
		}
		if(!MCMngr.accept(e)) {
			fprintf(stdout,"REJECTED SEQ %s ENERGY %8.3f %8.3f\n",sequence.c_str(),e,repackEnergy);
			if(!selPos->hideAllIdentitiesButOne(origIdentity)) {
				cerr << "Unable to restore " << origIdentity << " " << selPos->getPositionId() << endl;
				exit(0);
			}
		} else {
			fprintf(stdout,"ACCEPTED S%04d %s ENERGY %8.3f %8.3f\n",seqNum, sequence.c_str(),e,repackEnergy);
			char pdbName[50];
			sprintf(pdbName,"S%04d.pdb",seqNum);
			seqNum++;
			if(!sys.writePdb(string(pdbName))) {
				cerr << "Unable to write " << pdbName << endl;
				exit(0);
			}
		}
		
	}

	for(map<string,double>::iterator it = sequenceEnergy.begin(); it != sequenceEnergy.end(); it++) {
		cout << "SEQ " << it->first << " Energy " << it->second << endl;
	}

}
