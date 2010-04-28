/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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


#include <string>
#include <sstream>
#include <vector>
#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>


#include "runQuench.h"
#include "Quench.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"
#include "SystemRotamerLoader.h"
#include "PairwiseEnergyCalculator.h"
#include "AtomicPairwiseEnergy.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"

using namespace std;

using namespace MSL;



int main(int argc, char *argv[]){

	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	// Need to check for HIS residues.
	for (uint i = 0 ; i < initialSystem.positionSize();i++){
	  Residue &res = initialSystem.getPosition(i).getCurrentIdentity();
	  if (res.getResidueName() == "HIS" &&
	      res.atomExists("HE2") &&
	      res.atomExists("HD1")){
		fprintf(stdout,"Residue %10s is a protonated histidine change name to HSP\n",res.toString().c_str());
		res.setResidueName("HSP");
		continue;
	  }

	  if (res.getResidueName() == "HIS" &&
	      res.atomExists("HE2")){
		fprintf(stdout,"Residue %10s is a epsilon neutral histidine change name to HSE\n",res.toString().c_str());
		res.setResidueName("HSE");
	        continue;
	  }

	  if (res.getResidueName() == "HIS" &&
	      res.atomExists("HD1")){
		fprintf(stdout,"Residue %10s is a delta   neutral histidine change name to HSD\n",res.toString().c_str());
		res.setResidueName("HSD");
	        continue;
	  }

	  if (res.getResidueName() == "HIS"){
	    cout << "Residue: "<<res.toString()<< " is HIS but hydrogen configuration is not specified\n";
	    exit(2222);
	  }
	}

	Quench quencher(opt.topfile, opt.parfile, opt.rotlib);

	// Set the number of rotamers you want for large and small side chains, respectively
	quencher.setVariableNumberRotamers(opt.largeRotNum,opt.smallRotNum);

	// Set which positions you want to repack (optional)
	if (opt.positions.size() > 0) {	
		vector<int> variablePositions;
		for (uint i = 0; i < opt.positions.size(); i++) {
			//vector<string> pos = MslTools::tokenize(opt.positions[i],"_");
			variablePositions.push_back(initialSystem.getPositionIndex(opt.positions[i]));
		}
		System sys = quencher.runQuench(initialSystem,variablePositions);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << opt.outfile << endl;
		}
		writer.close();
	}else {

	  // Set which positions by automatically finding positions with large numbers of sidechain clashes, then quench those sidechains + other neighboring sidechains
	  if (opt.autoFindPositions){

	    set<int> variablePositions ;
	    set<int> ignorePositions ;  // Used for disulfide bonds. Otherwise CYS-CYS will be considered clashing, or they could be included in the "neighbor search".  This will keep the conformation of the CYS-CYS, forcing other nearby residues to accomodate it.

	    int clashTolerance = 1;
	    for (uint i = 0; i < initialSystem.positionSize();i++){
	      Position &pos1 = initialSystem.getPosition(i);

	      int numClashes = 0;
	      for (uint j = 0; j < initialSystem.positionSize();j++){
		if ( abs((int)j-(int)i) <= 2) continue; // skip direct neighbors

		Position &pos2 = initialSystem.getPosition(j);
		

		for (uint a1 = 0; a1 < pos1.atomSize();a1++){
		  if (pos1.getAtom(a1).getName()[0] == 'H') continue;
		  for (uint a2 = 0; a2 < pos2.atomSize();a2++){

		    if (pos2.getAtom(a2).getName()[0] == 'H') continue;

		    // Disulfides
		    if (pos1.getAtom(a1).getName() == "SG" &&
			pos2.getAtom(a2).getName() == "SG" &&
	                (pos1.getAtom(a1).distance(pos2.getAtom(a2)) < 2.2)){
	    		         ignorePositions.insert(i);
				 continue;
		    }
		    
		    if (pos1.getAtom(a1).distance(pos2.getAtom(a2)) < 2.2){ 	
		      numClashes++;
		      cout << "Clash between :"<<pos1.getAtom(a1).toString() << " and "<<pos2.getAtom(a2).toString()<<endl;
		    }

		  }
		}

		// Don't try more positions if already over clash limit
 		if (numClashes >= clashTolerance) break;

	      }

	      if (numClashes >= clashTolerance) {
		variablePositions.insert(i);
		
		// Also find neighobrs within 8 angstroms
		vector<int> neighbors = pos1.getCurrentIdentity().findNeighbors(8);
		for (uint n = 0; n < neighbors.size();n++){
		  variablePositions.insert(neighbors[n]);
		}

	      }

	    }


	    // Take the set difference between variable and ignore positions
	    vector<int> uniquePositions;
	    std::set_difference(variablePositions.begin(),variablePositions.end(),
				ignorePositions.begin(),ignorePositions.end(),
				std::inserter(uniquePositions,uniquePositions.end()));

	    if (uniquePositions.size()  < 2){
	      cerr << "ERROR 2222 AutoFindPositions code did not find any clashing residues\n";
	      exit(2222);
	    }
	
	    


	    cout << "Positions that clashed + their neighbors are: \n";
	    for (uint v = 0; v < uniquePositions.size();v++){
	      cout << initialSystem.getPosition(uniquePositions[v]).getCurrentIdentity().toString()<<endl;
	    }
	    cout <<endl;

	    System sys = quencher.runQuench(initialSystem,uniquePositions);

	    // Rename HIS residues.
	    /*
	    for (uint i = 0 ; i < sys.positionSize();i++){
	    	  Residue &res = sys.getPosition(i).getCurrentIdentity();
		  if (res.getResidueName() == "HSP" ||
		      res.getResidueName() == "HSD" ||
		      res.getResidueName() == "HSE"){

			fprintf(stdout,"Residue %10s is being renamed to HIS\n",res.toString().c_str());
			res.setResidueName("HIS");
			continue;
		  }
	    }
	    */
	    cout << "Write post quench pdb " << opt.outfile << endl;
	    sys.writePdb(opt.outfile);

	  } else {
		System sys = quencher.runQuench(initialSystem);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtomPointers())) {
			cerr << "Problem writing " << opt.outfile << endl;
		}
		writer.close();
	  }
	}


	cout << "Done."<<endl;
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
		cout << "runQuench --pdb PDB [--topfile TOPFILE --parfile PARFILE --rotlib ROTLIB --outfile OUTPDB --positions A_73 A_74 A_75 --smallRotNum NUM --largeRotNum NUM --autoFinPositions]\n";
		exit(0);
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		cerr << "WARNING no topfile specified, using default /library/charmmTopPar/top_all22_prot.inp\n";
		opt.topfile = "/library/charmmTopPar/top_all22_prot.inp";
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		cerr << "WARNING no parfile specified, using default /library/charmmTopPar/par_all22_prot.inp\n";
		opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
	}

	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "WARNING no rotlib specified, using default /library/rotlib/balanced/rotlib-balanced-200.txt\n";
		opt.rotlib = "/library/rotlib/balanced/rotlib-balanced-200.txt";
	}

	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		cerr << "WARNING no outfile specified, using default /tmp/currentConformation.pdb\n";
		opt.outfile = "/tmp/currentConformation.pdb";
	}

	opt.positions = OP.getStringVectorJoinAll("positions");

	opt.largeRotNum = OP.getInt("largeRotNum");
	if (OP.fail()){
		cerr << "WARNING largeRotNum not specified, using default 50\n";
		opt.largeRotNum = 50;
	}

	opt.smallRotNum = OP.getInt("smallRotNum");
	if (OP.fail()){
		cerr << "WARNING smallRotNum not specified, using default 5\n";
		opt.smallRotNum = 5;
	}
	if (opt.largeRotNum < 0 || opt.smallRotNum < 0) { cerr << "Need a positive rotamer number" << endl; exit(15); }


	opt.autoFindPositions = OP.getBool("autoFindPositions");
	if (OP.fail()){
	  opt.autoFindPositions = false;
	}


	return opt;
}
