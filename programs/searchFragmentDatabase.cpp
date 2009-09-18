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
#include "System.h"
#include "Timer.h"
#include "ChiStatistics.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "searchFragmentDatabase.h"

// STL Includes
#include<iostream>
using namespace std;




int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	Timer t;
	double start = t.getWallTime();


	// Read a list of PDBs into a single atom vector.
	AtomVector fragDB;
	fragDB.load_checkpoint(opt.database);


	for (uint i = 0 ; i < fragDB.size();i++){

		
	}


	PDBWriter pout;
	pout.open("/tmp/test.pdb");
	pout.write(fragDB);
	pout.close();

	cout <<"Done. took: "<<(t.getWallTime() - start)<<" seconds."<<endl<<endl;
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "searchFragmentDatabase --database FRAG_DB \n";
		exit(0);
	}

	opt.database = OP.getString("database");
	if (OP.fail()){
		cerr << "ERROR 1111 database not specified.\n";
		exit(1111);
	}

	return opt;
}



