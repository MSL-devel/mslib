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

using namespace MSL;


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
		vector<unsigned int> lastRotamerIndex(inputMask.size(),0);
		double energy = 0.0;
		for (uint i  = 0; i < inputMask.size();i++){

			int lastIndex = 0;
			int rotCount = 0;
			for (uint j = 0 ; j < inputMask[i].size();j++){
				if (inputMask[i][j]){
					rotCount++;
					ss << j<<":";
					lastIndex = j;
				}
			}

			// If at least 1 position has more than 1 rotamer we do not have a GMCC.
			if (rotCount != 1){
				globalMinimumFound = false;
			}
			fprintf(stdout,"POSITION %4d has %5d rotamers left\n",i,rotCount);
		}
		if (globalMinimumFound){

			cout << "GMEC: "<<ss.str()<<endl;
			if (opt.energyOpt.structureConfig != ""){
				// Create a system from the structural input options
				System sys;
				createSystem(opt.energyOpt.structOpt, sys);


				// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
				sys.setActiveRotamers(lastRotamerIndex);

				
				fprintf(stdout,"Energy GMEC: %8.3f\n",sys.getEnergySet()->calcEnergy());

				// Write out PDB
				char name[80];
				sprintf(name, "winnerMC-GMEC.pdb");
				sys.writePdb(name);
		
			}
			exit(0);
		}

		// Add Mask to MonteCarlo Object
		mc.setInputRotamerMasks(inputMask);
	}



	cout << "Read Energy Table"<<endl;
	mc.readEnergyTable(opt.energyTable);


	// Setup MC
	cout << "SetupMC"<<endl;
	
	mc.setNumberOfStoredConfigurations(opt.numStoredConfigurations);

	
	mc.setInitializationState(opt.initAlgo);


	// Random Seed
	mc.seed(opt.randomSeed);



	// Run MC
	cout << "Run MC with "<<mc.getNumPositions()<<" positions and "<<mc.getStateEnergy()<<" total energy."<<endl;
	mc.runMC(opt.annealStart,opt.annealEnd,opt.numCycles,opt.annealShape,opt.maxRejections,opt.deltaSteps,opt.minDeltaE);


	fprintf(stdout, "Random Seed Used: %d\n",mc.getSeed());

	cout << "After "<<(t.getWallTime()-startTime)<<" seconds"<<endl;

	// Either print rotamer selections + energies out, or generate PDBs and print out
	if (opt.energyOpt.structureConfig == ""){

		// Print results..
		cout << "Sampled Energies: "<<endl;
		mc.printSampledConfigurations();
	} else {

		// Create a system from the structural input options
		System sys;
		createSystem(opt.energyOpt.structOpt, sys);

		// Get priority queue of resulting conformations
		priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > &conformations = mc.getSampledConformations();

		int solution = 1;
		while (!conformations.empty()){
		        cout << "Working on solution conformation "<<solution<<" energy "<<conformations.top().first <<" "<<conformations.top().second<<endl;
			
			
			vector<string> toks = MslTools::tokenize(conformations.top().second,":");
			vector<unsigned int> rotamerState;
			for (uint i = 0; i < toks.size();i++){
			  rotamerState.push_back((unsigned int)MslTools::toInt(toks[i]));
			}
			conformations.pop();

			// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
			sys.setActiveRotamers(rotamerState);

			
			string sysString = "";
			for (uint c = 0; c < sys.chainSize();c++){
			  sysString += MslTools::stringf("%1s: %s\n",sys.getChain(c).getChainId().c_str(),PolymerSequence::toOneLetterCode(sys.getAtomPointers(),"CA").c_str());
			}

			fprintf(stdout,"Energy %04d: %8.3f\n%s\n",solution,sys.getEnergySet()->calcEnergy(),sysString.c_str());

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

	fprintf(stdout, "Random Seed Used: %d\n",mc.getSeed());
	cout << "Best-Sampled Energies after "<<(t.getWallTime()-startTime)<<" seconds."<<endl;
	mc.printSampledConfigurations();


	cout << "GoodBye."<<endl;
	exit(0);
}
