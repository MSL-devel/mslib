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
#include "System.h"
#include "Timer.h"
#include "RegEx.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "createFragmentDatabase.h"

// STL Includes
#include<iostream>
using namespace std;

using namespace MSL;

#include "MslOut.h"
static MslOut MSLOUT("createFragmentDatabase");


int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	Timer t;
	double start = t.getWallTime();

	MSLOUT.stream() << "READ LIST"<<endl;
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

	// Regular expression object
	RegEx re;
	re.setStringType(RegEx::SegID);  // Use SegID

	// Read a list of PDBs into a single atom vector.
	int totalMatches = 0;
	AtomPointerVector results;
	for (uint i = 0; i < pdbs.size();i++){
		
		MSLOUT.stream() << "Opening "<<pdbs[i]<<endl;
		System sys;
		sys.readPdb(pdbs[i]);
		

		for (uint c = 0; c < sys.chainSize();c++){
		  Chain &ch = sys.getChain(c);

		  if (opt.regex == ""){
		    for (uint a = 0; a < ch.atomSize();a++){
		      if (!opt.allAtom && ch.getAtom(a).getName() != "CA") continue;

		      Atom *tmp = new Atom(ch.getAtom(a));
		      tmp->setSegID(MslTools::getFileName(pdbs[i]));
		      results.push_back(tmp);
		      tmp = NULL;
		    }

		  } else {
		    vector<pair<int,int> > matches = re.getResidueRanges(ch,opt.regex);
		    MSLOUT.stream() << "Num Matches: "<<matches.size()<<endl;

		    for (uint m = 0; m < matches.size();m++){
			
		      bool sequential = true;
		      for (uint r = matches[m].first; r < matches[m].second;r++){
			if (r < matches[m].second-1 && ch.getResidue(r).getResidueNumber()+1 != ch.getResidue(r+1).getResidueNumber() && ch.getResidue(r).getResidueNumber()   != ch.getResidue(r+1).getResidueNumber()){
			MSLOUT.stream() << "Not sequential: "<<ch.getResidue(r).toString()<<" "<<ch.getResidue(r+1).toString()<<endl;
			    sequential = false;
			    break;
			}
			
		      }
		      
		      if (!sequential) { MSLOUT.stream() << "not sequential\n"; continue;}
		      totalMatches++;

		      for (uint r = matches[m].first; r < matches[m].second;r++){
			if (opt.allAtom){
			  for (uint a = 0; a < ch.getResidue(r).atomSize();a++){
			    Atom *tmp = new Atom(ch.getResidue(r).getAtom(a));
			    tmp->setSegID(MslTools::getFileName(pdbs[i]));
			    results.push_back(tmp);
			    tmp = NULL;
			  }
			} else {
			  if (ch.getResidue(r).atomExists("CA")){
			    Atom *tmp = new Atom(ch.getResidue(r)("CA"));
			    tmp->setSegID(MslTools::getFileName(pdbs[i]));
			    results.push_back(tmp);
			  }
			} // IF-ELSE opt.allAtom

		      } // END FOR matches[m] (a range of residues)

		    } // END FOR matches.size()

		  } // IF-ELSE regex == ""

		} // END FOR chain.size()
	
		MSLOUT.stream() << "\tsize: "<<results.size()<<endl;
	} // END FOR pdbs.size()

	if (opt.allAtom){
	  results.setName("allatom");
	} else {
	  results.setName("ca-only");
	}
	// Write out binary checkpoint file.
	results.save_checkpoint(opt.database);

	MSLOUT.stream() <<"Done. Got "<<totalMatches<<" total matches. took: "<<(t.getWallTime() - start)<<" seconds."<<endl<<endl;
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
		cout << "createFragmentDatabase --list LIST_OF_PDBS --database out.fragdb [--regex 'EELLLEE' --allAtoms]\n";
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

	opt.regex = OP.getString("regex");
	if (OP.fail()){
	  opt.regex = "";
	}

	opt.allAtom = OP.getBool("allAtoms");
	if (OP.fail()){
	  opt.allAtom = false;
	}
	cout << OP<<endl;
	return opt;
}



