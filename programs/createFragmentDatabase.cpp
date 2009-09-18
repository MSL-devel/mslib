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
#include "createFragmentDatabase.h"

// STL Includes
#include<iostream>
using namespace std;


int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	Timer t;
	double start = t.getWallTime();

	cout << "READ LIST"<<endl;
	vector<string> pdbs;  
	ifstream fs;

	fs.open(opt.list.c_str());
	if (fs.fail()){
		cerr<<"Cannot open file "<<opt.list<<endl;
		exit(1);
	}

	while(true){
		string line;
		getline(fs, line);

		if(fs.fail()){
			//no more lines to read, quite the while.
			break;
		}

		if(line==""){
			continue;
		}
		pdbs.push_back(line);
	}

	fs.close();


	// Read a list of PDBs into a single atom vector.
	AtomVector results;
	for (uint i = 0; i < pdbs.size();i++){
		
		cout << "Opening "<<pdbs[i]<<endl;
		PDBReader pin;
		pin.open(pdbs[i]);
		pin.read();
		pin.close();


		AtomVector &tmp = pin.getAtoms();
		cout << "\tfound "<<tmp.size()<<" atoms in file."<<endl;
		for (uint a = 0; a < tmp.size();a++){

			if (tmp(a).getName() == "CA"){
				results.push_back(new Atom(tmp(a)));
				results.back()->setSegID(MslTools::getFileName(pdbs[i]));
			}
		}
		
	       
	}

	cout << "Time to write checkpoint file"<<opt.database<< " with "<<results.size()<<" atoms."<<endl;

	// Write out binary checkpoint file.
	results.save_checkpoint(opt.database);

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
		cout << "createFragmentDatabase --list LIST_OF_PDBS --database out.fragdb\n";
		exit(0);
	}

	opt.list = OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}

	opt.database = OP.getString("database");
	if (OP.fail()){
		cerr << "ERROR 1111 database not specified.\n";
		exit(1111);
	}

	return opt;
}



