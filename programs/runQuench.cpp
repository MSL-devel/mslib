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
#include "runQuench.h"

int main(int argc, char *argv[]){
	// Option Parser
	Options opt = setupOptions(argc,argv);

	if(opt.errorFlag) {
		cout << opt.OPerrors ;
		cout << opt.errorMessages ;
		exit(0);
	}

	if(opt.warningFlag) {
		cout << opt.warningMessages;
	}

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	Quench quencher(opt.topfile, opt.parfile, opt.rotlib);

	// Set the number of rotamers you want for large and small side chains, respectively
	if (opt.rotLevel != "") {
		quencher.setRotamerLevel(opt.rotLevel);
	} else {
		quencher.setVariableNumberRotamers(opt.largeRotNum,opt.smallRotNum);
	}

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

