/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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
// MSL Includes
#include "MslTools.h"
#include "energyOptimizations.h"
#include "Timer.h"
#include "OptionParser.h"
#include "MonteCarloOptimization.h"
#include "DeadEndElimination.h"


// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace std;
using namespace MslTools;

// Need a global MC object
MonteCarloOptimization mc;
Timer t;
double startTime =  MslTools::doubleMax;

int main(int argc, char *argv[]) {

	MonteCarloOptions opt = setupMonteCarloOptions(argc, argv);



	startTime = t.getWallTime();

	// Setup a signal catching function
	signal(SIGINT, cleanExit);
	signal(SIGTERM, cleanExit);


	// Dead End Elimination
	if (opt.DEE){
		DeadEndElimination dee;

		// Read energy table..
		dee.readEnergyTable(opt.energyTable);

		// Run Goldstein Singles
		dee.runSimpleGoldsteinSingles();

		// Get Eliminated Rotamers
		vector<vector<bool> > inputMask;
		inputMask = dee.getMask();

		// Inform about the number of eliminated rotamers..
		cout << "Eliminated Rotamers: "<<dee.getEliminatedCounter()<<" out of "<<dee.getTotalNumberRotamers()<<endl;
		stringstream ss;
		bool globalMinimumFound = true;
		for (uint i  = 0; i < inputMask.size();i++){

			int rotCount = 0;
			for (uint j = 0 ; j < inputMask[i].size();j++){
				if (inputMask[i][j]){
					rotCount++;
					ss << j<<":";

				}
			}

			// If at least 1 position has more than 1 rotamer we do not have a GMEC.
			if (rotCount != 1){
				globalMinimumFound = false;
			}
			fprintf(stdout,"POSITION %4d has %5d rotamers left\n",i,rotCount);
		}
		if (globalMinimumFound){

			cout << "GMEC: "<<ss.str()<<endl;
			exit(0);
		}

		// Add Mask to MonteCarlo Object
		mc.setInputRotamerMasks(inputMask);
	}



	cout << "Read Energy Table"<<endl;
	mc.readEnergyTable(opt.energyTable);

	// Setup MC
	cout << "SetupMC"<<endl;
	mc.setNumberOfCycles(opt.numCycles);
	mc.setNumberOfStoredConfigurations(opt.numStoredConfigurations);

	mc.setAnnealSchedule(MonteCarloOptimization::LIN_TEMP_ANNEAL,opt.annealStart,opt.annealEnd);

	if (opt.annealType == "EXPONENTIAL"){
		mc.setAnnealSchedule(MonteCarloOptimization::EXP_TEMP_ANNEAL,opt.annealStart,opt.annealEnd);
	}

	if (opt.annealType == "LINEAR_CYCLES"){
		mc.setAnnealSchedule(MonteCarloOptimization::SAWTOOTH_TEMP_ANNEAL,opt.annealStart,opt.annealEnd,opt.numberOfAnnealCycles);
	}

	if (opt.annealType == "EXPONENTIAL_CYCLES"){
		mc.setAnnealSchedule(MonteCarloOptimization::EXPCYCLE_TEMP_ANNEAL,opt.annealStart,opt.annealEnd,opt.numberOfAnnealCycles);
	}


	
	mc.setInitializationState(MonteCarloOptimization::LOWESTSELF);

	if (opt.initAlgorithm == "RANDOM"){
		mc.setInitializationState(MonteCarloOptimization::RANDOM);
	}
	if (opt.initAlgorithm == "QUICKSCAN"){
		mc.setInitializationState(MonteCarloOptimization::QUICKSCAN);
	}
	if (opt.initAlgorithm == "USERDEF"){
		mc.setInitializationState(MonteCarloOptimization::USERDEF, opt.initConfiguration);
	}


	// Random Seed
	mc.setRandomSeed(opt.randomSeed);



	// Run MC
	cout << "Run MC"<<endl;
	mc.runMC();


	fprintf(stdout, "Random Seed Used: %d\n",mc.getRandomSeed());



	

	cout << "After "<<mc.getCurrentStep()<<" steps and "<<(t.getWallTime()-startTime)<<" seconds"<<endl;




	// Either print rotamer selections + energies out, or generate PDBs and print out
	if (opt.structureConfig == ""){

		// Print results..
		cout << "Sampled Energies: "<<endl;
		mc.printSampledConfigurations();
	} else {

		// Create a system from the structural input options
		System sys;
		createSystem(opt.structOpt, sys);

		// Get priority queue of resulting conformations
		priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > &conformations = mc.getSampledConformations();

		int solution = 1;
		while (!conformations.empty()){
			cout << "Working on solution conformation "<<solution<<conformations.top().first <<" "<<conformations.top().second<<endl;
			
			vector<string> toks = MslTools::tokenize(conformations.top().second);
			vector<int> rotamerState;
			for (uint i = 0; i < toks.size();i++){
				rotamerState.push_back(MslTools::toInt(toks[i]));
			}
			conformations.pop();

			// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
			changeRotamerState(opt.structOpt,sys,rotamerState);

			// Write out PDB
			char name[80];
			sprintf(name, "winnerMC-%04d.pdb",solution++);
			sys.writePdb(name);
		
		}
				
		
		
				
	}
	
}






void cleanExit(int sig) {

	cout << "*************SIGNAL CAUGHT*****************\n";
	cout << "\tSIG"<<sig<<endl;

	fprintf(stdout, "Random Seed Used: %d\n",mc.getRandomSeed());
	cout << "Best-Sampled Energies after "<<mc.getCurrentStep()<<" steps and "<<(t.getWallTime()-startTime)<<" seconds."<<endl;
	mc.printSampledConfigurations();


	cout << "GoodBye."<<endl;
	exit(0);
}
