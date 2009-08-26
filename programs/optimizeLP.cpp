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
#include "LinearProgrammingOptimization.h"
#include "MslTools.h"
#include "optimizeLP.h"

// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace std;
using namespace MslTools;

// Need a global MC object
LinearProgrammingOptimization lp;
Timer t;
double startTime =  MslTools::doubleMax;

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);


	startTime = t.getWallTime();

	// Setup a signal catching function
	signal(SIGINT, cleanExit);
	signal(SIGTERM, cleanExit);


	cout << "Read Energy Table"<<endl;
	lp.readEnergyTable(opt.energyTable);

	// Setup LP
	cout << "SetupLP"<<endl;
	lp.analyzeEnergyTable();
	lp.createLP();

	// Solve LP
	lp.solveLP();

	// Print Results
	lp.printMe();
	cout << "Optimized Energy: "<<lp.getTotalEnergy()<<"\t"<<lp.getRotString()<<endl;

	cout << "Process took "<<(t.getWallTime()-startTime)<<" seconds"<<endl;
	
	
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
        cout << "optimizeLP CONF\n";
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
	    
	    cout << "# Options for optimizeLP\n\n";
	    cout << "# Energy Table\n";
	    cout << "energyTable energy.txt\n\n";
	    exit(0);
        
    }

    opt.debug = OP.getBool("debug");

    opt.energyTable = OP.getString("energyTable");
    if (OP.fail()){
	    cerr << "ERROR 1111 energyTable not specified."<<endl;	
	    exit(1111);
    }


    cout << OP<<endl;
    return opt;
}

void cleanExit(int sig) {

	cout << "*************SIGNAL CAUGHT*****************\n";
	cout << "\tSIG"<<sig<<endl;

	cout << "Print Current Flags"<<endl;
	lp.printMe();
	lp.getTotalEnergy();

	cout << "GoodBye."<<endl;
	exit(0);
}
