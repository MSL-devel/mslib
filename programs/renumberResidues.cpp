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
#include "MslTools.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "FastaReader.h"
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


	//AtomSelection asel(pdb.getAtomPointers());
	//System asys;
	//asys.addAtoms(asel.select("chain A"));
	//pdb = asys;

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
	    ats(i).setResidueIcode("");
	    //ats(i).setChainId("C");
	  }
	  PDBWriter pdb_out;
	  pdb_out.open(opt.outPdb);
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

	  //	  if (opt.useSeqPatterns){
	    
	    // Find sequence patterns begin,end ; begin, end
	    // Renumber from ref upto first sequence pattern, insertion codes until next sequence pattern ; repeat till all sequence patterns are done.



	    //	  }


	  if (opt.useDistance){
	    System newSys;
	    
	    vector<int> noMatch;
	    int newResidueNumber = opt.startRes;
	    string insertionCodes = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
	    int insertionCodeIndex = 0;
	    bool lastMatchWasReal = false;
	    for (int p = 0; p < pdb.positionSize();p++){
	      string matchId = "";
	      bool matchFound = false;
	      cout << "P: "<<p<<" "<<pdb.getPosition(p).toString()<<endl;
	      for (uint r = 0; r < ref.positionSize();r++){
		//cout << "\tR: "<<r<<" "<<ref.getPosition(r).toString()<<endl;

		if (ref.getPosition(r).atomExists("CA") && 
		    pdb.getPosition(p).atomExists("CA")) {
		  double dist = ref.getPosition(r).getAtom("CA").distance(pdb.getPosition(p).getAtom("CA"));
		  if (dist < 1.2){

		  // Last residue was not a match, consider this to be not a match
		  if (noMatch.size() > 0 && abs(noMatch.back()-p) == 1){
		    if (!lastMatchWasReal){
		      lastMatchWasReal = true;		      
		      matchFound = false;
		      break;
		    }
		  }
		  lastMatchWasReal = false;

		  newResidueNumber = ref.getPosition(r).getResidueNumber();
		  fprintf(stdout,"Setting %8s to %1s,%4d\n",pdb.getPosition(p).getPositionId().c_str(),ref.getPosition(r).getChainId().c_str(),newResidueNumber);
		  pdb.getPosition(p).setResidueNumber(newResidueNumber);
		  pdb.getPosition(p).setResidueIcode(ref.getPosition(r).getResidueIcode());		  
		  pdb.getPosition(p).setChainId(ref.getPosition(r).getChainId());
		  newSys.addAtoms(pdb.getPosition(p).getAtomPointers());

		  matchFound = true;
		  break;
		  } // DIST < 1.0
		
		} // CA exists in pdb and ref positions..
	      } // END FOR REF.POS

	      if (!matchFound){


		// Add using Icodes
		if (opt.useIcodes){

		  // New missing segment
		  if (noMatch.size() == 0 || noMatch.back()+1 != p){
		    if (noMatch.size() != 0){
		      cout << "noMactch.back() == "<<noMatch.back()<<" p = "<<p<<endl;
		    }
		    newResidueNumber++;
		    insertionCodeIndex = 0;
		  } else {
		    insertionCodeIndex++;
		  }


		  if (p != 0){
		    if (insertionCodeIndex > insertionCodes.size()){
		      cerr << "ERROR too many insertions. died at position: "<<p<<" "<<pdb.getPosition(p).toString()<<endl;
		      exit(234);
		    }
		  string iCode = insertionCodes.substr(insertionCodeIndex,1);
		  fprintf(stdout,"Setting %8s to %1s,%4d,%s\n",pdb.getPosition(p).getPositionId().c_str(),pdb.getPosition(p-1).getChainId().c_str(),newResidueNumber,iCode.c_str());
		  pdb.getPosition(p).setResidueNumber(newResidueNumber);
		  pdb.getPosition(p).setResidueIcode(iCode);		  
		  pdb.getPosition(p).setChainId(pdb.getPosition(p-1).getChainId());
		  newSys.addAtoms(pdb.getPosition(p).getAtomPointers());

		  }		  
		}  

		noMatch.push_back(p);

	      }

	    }
	    
	    if (!opt.useIcodes){

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
	    } // if !opt.useIcodes


	    pdb = newSys;

	  } else {
	    
	    if (opt.fasta != ""){
	      
	        FastaReader fin(opt.fasta);
		fin.open();
		if (!fin.read()){
		  cerr << "ERROR reading "<<opt.fasta<<endl;
		  exit(143);
		}
		fin.close();

		string refSeq = fin.getSequence(MslTools::getFileName(opt.ref));
		if (refSeq == ""){
		  cerr << "ERROR couldn't find "<<MslTools::getFileName(opt.ref)<<" as a name in : "<<opt.fasta<<endl;
		  exit(2342);
		}

		string pdbSeq = fin.getSequence(MslTools::getFileName(opt.pdb));
		if (pdbSeq == ""){
		  cerr << "ERROR couldn't find "<<MslTools::getFileName(opt.pdb)<<" as a name in : "<<opt.fasta<<endl;
		  exit(2342);
		}

		if (refSeq.size() != pdbSeq.size()){
		  cerr << "ERROR refSeq and pdbSeq do not have the same length: "<<refSeq.size()<<" , "<<pdbSeq.size()<<endl;
		  exit(2342);
		}

		int nonDashRefCounter = 0;
		int nonDashPdbCounter = 0;
		int newResNum   = ref.getPosition(0).getResidueNumber();
		string newChain = ref.getPosition(0).getChainId();
		string newIcode = ref.getPosition(0).getResidueIcode();
		int sequentialInsertion = 0;
		for (uint i = 0; i < pdbSeq.size();i++){

		  // Skip over '-' in pdbSeq
		  if (pdbSeq[i] == '-') {

		    // Increment nonDash counter
		    if (refSeq[i] != '-'){
		      nonDashRefCounter++;
		    }

		    continue;
		  }

		  // Need to de-couple positions and chains... due to setChainId in Position object.
		  pdb.getPosition(nonDashPdbCounter).setParentChain(NULL);

		  if (refSeq[i] == '-'){
		    sequentialInsertion++;

		    if (sequentialInsertion == 26){
		      newIcode = "0";
		    }

		    fprintf(stdout, "Renumber %10s ", pdb.getPosition(nonDashPdbCounter).getPositionId().c_str());
		    pdb.getPosition(nonDashPdbCounter).setResidueNumber(newResNum);
		    pdb.getPosition(nonDashPdbCounter).setResidueIcode(newIcode);
		    pdb.getPosition(nonDashPdbCounter).setChainId(newChain); // be careful this may set for whole system/chain. need to check this for multi-chain files

		    fprintf(stdout,"to %10s [ %10s ] INSERTION\n",pdb.getPosition(nonDashPdbCounter).getPositionId().c_str(),ref.getPosition(nonDashRefCounter).getPositionId().c_str());
		    
		    newIcode[0]++;
		  } else {


		    fprintf(stdout, "Renumber %10s ", pdb.getPosition(nonDashPdbCounter).getPositionId().c_str());
		    string chain  = ref.getPosition(nonDashRefCounter).getChainId();
		    int resNum    = ref.getPosition(nonDashRefCounter).getResidueNumber();
		    string icode  = ref.getPosition(nonDashRefCounter).getResidueIcode();

		    pdb.getPosition(nonDashPdbCounter).setResidueNumber(resNum);
		    pdb.getPosition(nonDashPdbCounter).setResidueIcode(icode);
		    pdb.getPosition(nonDashPdbCounter).setChainId(chain); // be careful this may set for whole system/chain. need to check this for multi-chain files
		    fprintf(stdout, "to %10s [ %10s ]\n",pdb.getPosition(nonDashPdbCounter).getPositionId().c_str(),ref.getPosition(nonDashRefCounter).getPositionId().c_str());
		    
		    newIcode = "A";
		    if (icode != ""){
		      newIcode = icode;
		      newIcode[0]++;
		    }

		    sequentialInsertion = 0;
		    newResNum = resNum;
		    newChain = chain;
		    nonDashRefCounter++;
		  }

		  nonDashPdbCounter++;

		} // END FOR pdbSeq		

	      
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
	}

	
	pdb.writePdb(opt.outPdb);
	
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

	opt.useIcodes = OP.getBool("useIcodes");
	if (OP.fail()){
	  opt.useIcodes = false;
	}

	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  opt.outPdb = MslTools::stringf("%s_renum.pdb",MslTools::getFileName(opt.pdb).c_str());
	}

	opt.fasta = OP.getString("fasta");
	if (OP.fail()){
	  opt.fasta = "";
	}

	return opt;
}



