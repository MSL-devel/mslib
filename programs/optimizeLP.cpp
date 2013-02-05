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
#include "Timer.h"
#include "OptionParser.h"
#include "LinearProgrammingOptimization.h"
#include "MslTools.h"
#include "energyOptimizations.h"

// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace std;

using namespace MSL;


// Need a global MC object
LinearProgrammingOptimization lp;
LinearProgrammingOptions opt;
Timer t;
double startTime =  MslTools::doubleMax;

int main(int argc, char *argv[]) {

	opt = setupLinearProgrammingOptions(argc, argv);


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
	

	// If strucure configuration has been specified, output solution PDBs
	if (opt.structureConfig != ""){
		
		// Create a system from the structural input options
		System sys;
		createSystem(opt.structOpt, sys);

		vector<unsigned int> rotamerState = lp.getRotamerSelection();

		// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
		sys.setActiveRotamers(rotamerState);
		sys.writePdb("winnerLP.pdb");
	}
}




void cleanExit(int sig) {

	cout << "*************SIGNAL CAUGHT*****************\n";
	cout << "\tSIG"<<sig<<endl;

	cout << "Print Current Flags"<<endl;
	lp.printMe();
	lp.getTotalEnergy();

	// If strucure configuration has been specified, output solution PDBs
	if (opt.structureConfig != ""){
		
		// Create a system from the structural input options
		System sys;
		createSystem(opt.structOpt, sys);

		vector<unsigned int> rotamerState = lp.getRotamerSelection();

		// Helper function takes structOptions, a System and a rotamer state , putting system into given rotamer state.
		sys.setActiveRotamers(rotamerState);
		sys.writePdb("winnerLP.pdb");
	}

	cout << "GoodBye."<<endl;
	exit(0);
}
