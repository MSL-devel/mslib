/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
 Sabareesh Subramaniam, Ben Mueller

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

using namespace MSL;



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
	AtomPointerVector results;
	for (uint i = 0; i < pdbs.size();i++){
		
		cout << "Opening "<<pdbs[i]<<endl;
		PDBReader pin;
		pin.open(pdbs[i]);
		pin.read();
		pin.close();

		
		AtomPointerVector &tmp = pin.getAtomPointers();
		int totalNumber = results.size();
		map<int, bool> nearALoopMap;
		for (uint a = 0; a < tmp.size();a++){

		  bool isALoop = false;
		  bool isNearALoop = false;
		  bool isChainTerminal = false;

		  if (a > 1 && a < tmp.size()-1 && tmp(a).getChainId() == tmp(a-1).getChainId() && tmp(a).getChainId() == tmp(a+1).getChainId()){ isChainTerminal = true; }



		  if (tmp(a).getSegID() == "LLLL" || tmp(a).getSegID() == "TTTT" || tmp(a).getSegID() == "SSSS"){
		    isALoop = true;
		    nearALoopMap[a] = true;
		  }

		  // Look ahead 4 residues for loop residues
		  for (uint a2 = a+1; a2 < a+5; a2++){
		    if (a2 < tmp.size()){
		      if (tmp(a).getChainId() == tmp(a2).getChainId()){
			if (tmp(a2).getSegID()  == "LLLL" || tmp(a2).getSegID() == "TTTT" || tmp(a2).getSegID() == "SSSS"){


			  nearALoopMap[a2] = true;
			  isNearALoop = true;
			}
		      } else {
			break;
		      }
		    }
		  }

		  // If near an N-term loop
		  if (!isChainTerminal && (a < 4 || nearALoopMap[a-4] || nearALoopMap[a-3] ||nearALoopMap[a-2] ||nearALoopMap[a-1])){
		    isNearALoop = true;
		  }


		  if (tmp(a).getName() == "CA" && (isALoop || isNearALoop)){ 
				results.push_back(new Atom(tmp(a)));
				results.back()->setSegID(MslTools::getFileName(pdbs[i]));
			}
		}
		
		cout << "\tfound "<<tmp.size()<<" atoms in file. Saved "<<results.size()-totalNumber<<" atoms. In Total: "<<results.size()<<endl;		
	       
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



