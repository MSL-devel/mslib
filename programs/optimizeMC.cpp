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
#include "Timer.h"
#include "OptionParser.h"
#include "MonteCarloOptimization.h"
#include "DeadEndElimination.h"
#include "MslTools.h"
#include "optimizeMC.h"

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

	Options opt = setupOptions(argc, argv);



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

	// Print results..
	cout << "Sampled Energies: "<<endl;
	mc.printSampledConfigurations();
	

	cout << "After "<<mc.getCurrentStep()<<" steps and "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	
}




Options setupOptions(int theArgc, char * theArgv[]){

    // Create the options
    Options opt;
    

    // Parse the options
    OptionParser OP;
    OP.readArgv(theArgc, theArgv);
    OP.setRequired(opt.required);    
    OP.setAllowed(opt.optional);
    //    OP.setShortOptionEquivalent(opt.equivalent);
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"


    if (OP.countOptions() == 0){
        cout << "Usage:" << endl;
        cout << endl;
        cout << "optimizeMC CONF\n";
        exit(0);
    }

    opt.configfile = OP.getString("configfile");
    
    if (opt.configfile != "") {
        OP.readFile(opt.configfile);
        if (OP.fail()) {
            string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
            cerr << "ERROR 1111 "<<errorMessages<<endl;
        }
    }

    if (OP.getBool("help")){
	    
	    cout << "# Options for optimizeMC\n\n";
	    cout << "# Energy Table\n";
	    cout << "energyTable energy.txt\n\n";
	    cout << "# Annealing Schedule\n";
	    cout << "#     can be LINEAR, EXPONENTIAL, LINEAR_CYCLES, EXPONENTIAL_CYCLES\n";
	    cout << "annealScheduleType EXPONENTIAL\n";
	    cout << "annealScheduleStartTemp 1000\n";
	    cout << "annealScheduleEndTemp   1 \n\n";
	    cout << "numberOfAnnealCycles    1 \n\n";
	    cout << "# Number of MC Cycles \n";
	    cout << "numberOfCycles 10000\n\n";
	    cout << "# Initial Rotamer Configuration set via initializationAlgorithm\n";
	    cout << "#     can be RANDOM LOWESTSELF QUICKSCAN USERINPUT\n";
	    cout << "initializationAlgorithm LOWESTSELF\n\n";
	    cout << "# Initial Configuration set by user, RotamerNumberForPosition1:RotamerNumberForPosition2:\n";
	    cout << "#   Notice the last character is a ':'\n";
	    cout << "# initializationConfiguration 0:0:0:0:\n\n";
	    cout << "# Number of Rotamer Configurations to report , top 100 lowest energy scoring found\n";
	    cout << "numberOfStoredConfigurations 100\n\n";
	    cout << "# Random Seed to use, -1 means create a time-based seed\n";
	    cout << "randomSeed 838201\n\n";
	    exit(0);
        
    }

    opt.debug = OP.getBool("debug");

    opt.energyTable = OP.getString("energyTable");
    if (OP.fail()){
	    cerr << "ERROR 1111 energyTable not specified."<<endl;	
	    exit(1111);
    }


    opt.annealType = OP.getString("annealScheduleType");
    if (OP.fail()){
	    opt.annealType = "LINEAR";
    }

    opt.annealStart = OP.getDouble("annealScheduleStartTemp");
    if (OP.fail()){
	    opt.annealStart = 1000;
    }

    opt.annealEnd = OP.getDouble("annealScheduleEndTemp");
    if (OP.fail()){
	    opt.annealEnd = 1;
    }

    opt.numberOfAnnealCycles = OP.getInt("numberOfAnnealCycles");
    if (OP.fail()){
	    opt.numberOfAnnealCycles = 1;
    }

    opt.numCycles = OP.getInt("numberOfCycles");
    if (OP.fail()){
	    opt.numCycles = 10000;
    }

    opt.initAlgorithm = OP.getString("initializationAlgorithm");
    if (OP.fail()){
	    opt.initAlgorithm = "LOWESTSELF";
    }
    opt.initConfiguration = OP.getString("initializationConfiguration");
    if (OP.fail()){
	    opt.initConfiguration = "";
    }

    opt.numStoredConfigurations = OP.getInt("numberOfStoredConfigurations");
    if (OP.fail()){
	    opt.numStoredConfigurations = 50;
    }
    
    opt.randomSeed = OP.getInt("randomSeed");
    if (OP.fail()){
	    opt.randomSeed = -1;
    }

    
    opt.DEE  = OP.getBool("DEE");
    if (OP.fail()){
	    opt.DEE = false;
    }

    cout << OP<<endl;
    return opt;
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
