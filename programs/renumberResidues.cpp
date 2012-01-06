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
#include "MslTools.h"
#include "OptionParser.h"
#include "renumberResidues.h"
#include "MslOut.h"

// STL Includes
#include<iostream>
using namespace std;

using namespace MSL;


// MslOut 
static MslOut MSLOUT("renumberResidues");


int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	System pdb;
	pdb.readPdb(opt.pdb);

	if (opt.preserveFileOrder){
	  cout << "PRESERVING FILE ORDER\n";
	  PDBReader pdb_in;
	  pdb_in.open(opt.pdb);
	  pdb_in.read();
	  pdb_in.close();

	  AtomPointerVector ats = pdb_in.getAtomPointers();
	  int residue = opt.startRes-1;
	  string lastPosId = "";
	  for (uint i = 0; i < ats.size();i++){
	    string curPosId = MslTools::getPositionId(ats(i).getChainId(),ats(i).getResidueNumber(),ats(i).getResidueIcode());

	    if (lastPosId != curPosId){
	      residue++;
	      lastPosId = curPosId;
	    }

	    cout << "Setting "<<ats(i).toString()<< " to "<<residue<<endl;
	    ats(i).setResidueNumber(residue);
	  }
	  PDBWriter pdb_out;
	  pdb_out.open(MslTools::stringf("%s_renum.pdb",MslTools::getFileName(opt.pdb).c_str()));
	  pdb_out.write(ats);
	  pdb_out.close();
	  exit(0);
	  
	} else if (opt.ref == ""){
	  int residue = opt.startRes;
	  for (uint p =0; p < pdb.positionSize();p++){
	    pdb.getPosition(p).setResidueNumber(residue);
	    residue++;
	  }
	} else {

	  System ref;
	  ref.readPdb(opt.ref);


	  if (opt.useDistance){
	    System newSys;
	    
	    vector<int> noMatch;

	    for (uint p = 0; p < pdb.positionSize();p++){
	      string matchId = "";
	      bool matchFound = false;

	      for (uint r = 0; r < ref.positionSize();r++){
		double dist = ref.getPosition(r).getAtom("CA").distance(pdb.getPosition(p).getAtom("CA"));
		if (ref.getPosition(r).atomExists("CA") && 
		    pdb.getPosition(p).atomExists("CA") &&
		    dist < 1.0){

		  fprintf(stdout,"Setting %8s to %1s,%4d\n",pdb.getPosition(p).getPositionId().c_str(),ref.getPosition(r).getChainId().c_str(),ref.getPosition(r).getResidueNumber());
		  pdb.getPosition(p).setResidueNumber(ref.getPosition(r).getResidueNumber());
		  pdb.getPosition(p).setResidueIcode(ref.getPosition(r).getResidueIcode());		  
		  pdb.getPosition(p).setChainId(ref.getPosition(r).getChainId());
		  newSys.addAtoms(pdb.getPosition(p).getAtomPointers());

		  matchFound = true;
		  break;
		} // END IF REF = PDB POSITION 


	      } // END FOR REF.POS

	      if (!matchFound){
		noMatch.push_back(p);
	      }

	    }


	    // Assume that a neighboring position is defined. This will work for a loop with missing residues where the ref PDB has approriate numbering...(i.e. skipping loop residue numbers).
	    for (uint i = 0; i < noMatch.size();i++){
	      Position &pos = pdb.getPosition(noMatch[i]);

	     

	      //string newChain = pdb.getPosition(noMatch[i]-1).getChainId();
	      //int newResNum = pdb.getPosition(noMatch[i]-1).getResidueNumber()+1;
	      int closestId = 0;
	      double closestDist = MslTools::doubleMax;
	      for (uint j = 0; j < newSys.positionSize();j++){
		double dist = newSys.getPosition(j).getAtom("CA").distance(pos.getAtom("CA"));
		if (newSys.getPosition(j).atomExists("CA") && 
		    pos.atomExists("CA") &&
		    dist < closestDist){

		  // If we find a distance that as just as close, use the first one we found rather than this one..
		  if (abs(dist - closestDist) > 0.5){
		    closestDist = dist;
		    closestId = j;
		  }

		}
	      }
	      
	      string newChain = newSys.getPosition(closestId).getChainId();
	      int newResNum = newSys.getPosition(closestId).getResidueNumber()+1;

	      fprintf(stdout,"Missing ref position: %8s setting to %1s,%4d\n",pos.toString().c_str(),newChain.c_str(),newResNum);
	      pos.setResidueNumber(newResNum);
	      pos.setResidueIcode("");
	      pos.setChainId(newChain);
	      newSys.addAtoms(pos.getAtomPointers());
	    }
	    pdb = newSys;

	  } else {

	    if (ref.positionSize() != pdb.positionSize()){	
	      cerr << "Ref position size is: "<<ref.positionSize()<<" pdb position size is: "<<pdb.positionSize()<<endl;
	      exit(111);
	    }
	    for (uint p = 0; p < ref.positionSize();p++){
	      pdb.getPosition(p).setResidueNumber(ref.getPosition(p).getResidueNumber());
	      pdb.getPosition(p).setResidueIcode(ref.getPosition(p).getResidueIcode());
	    }
	  }
	}

	
	pdb.writePdb(MslTools::stringf("%s_renum.pdb",MslTools::getFileName(opt.pdb).c_str()));
	
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);


	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "renumberResidues --pdb foo.pdb --startRes NUM \n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}


	opt.startRes = OP.getInt("startRes");
	if (OP.fail()){
	  opt.startRes = 1;
	}

	opt.ref = OP.getString("ref");
	if (OP.fail()){
	  opt.ref = "";
	}

	opt.preserveFileOrder = OP.getBool("preserveFileOrder");
	if (OP.fail()){
	  opt.preserveFileOrder = false;
	}

	opt.useDistance = OP.getBool("useDistance");
	if (OP.fail()){
	  opt.useDistance = false;
	}
	return opt;
}



