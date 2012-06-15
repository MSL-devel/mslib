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
#include "Transforms.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "searchFragmentDatabase.h"
#include "RegEx.h"
#include "MslOut.h"
#include "SasaCalculator.h"

// STL Includes
#include<iostream>
using namespace std;

using namespace MSL;


// MslOut 
static MslOut MSLOUT("searchFragmentDatabase");


int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	Timer tm;
	double start = tm.getWallTime();

	System ref;
	ref.readPdb(opt.pdb);
	
	AtomSelection sel(ref.getAtomPointers());
	//AtomPointerVector refCa = sel.select("resn ALA+CYS+ASP+GLU+PHE+GLY+HIS+ILE+LYS+LEU+MET+ASN+PRO+GLN+ARG+SER+THR+VAL+TRP+TYR and name CA");
	AtomPointerVector refCa = sel.select("name CA");

	// Read a list of PDBs into a single atom vector.
	MSLOUT.stream()<<"Load FragDB: "<<opt.database<<endl;
	AtomPointerVector fragDB;
	fragDB.load_checkpoint(opt.database);

	MSLOUT.stream()<<fragDB.size()<<" residues loaded"<<endl;
	Transforms t;
	int match = 1;
	for (uint i = 0 ; i < fragDB.size()-refCa.size()-1;i++){

	  Atom &at1 = fragDB(i);
	  Atom &at2 = fragDB(i+refCa.size()-1);


	  // Different PDB files
	  if (at1.getSegID() != at2.getSegID()) {
	    //MSLOUT.stream()<< "Working on "<<at1.getSegID()<<endl;
	    continue;
	  }

	  // Different chains
	  if (at1.getChainId() != at2.getChainId()){
	    continue;
	  }
	  //MSLOUT.stream()<< "Working on "<<at1.toString()<<" *** "<<at2.toString()<<endl;

	  if (opt.regex != ""){

	    bool foundRegEx = false;
	    stringstream ss;
	    for (uint j = i;j <= i+refCa.size()-1;j++){
	      ss << MslTools::getOneLetterCode(fragDB(j).getResidueName());
	    }
	    //MSLOUT.stream() << "Seq: "<<ss.str()<<" regex: "<<opt.regex<<endl;

	    vector<string> results; // don't store any results, but required arguement
	    if (MslTools::regex(ss.str(),opt.regex,results)){
		foundRegEx = true;
	    }


	    if (!foundRegEx){
	      continue;
	    } 

	    MSLOUT.stream()<< "FOUND REGEX "<<at1.getSegID()<<" "<<at1.getChainId()<<" "<<at1.getResidueNumber()<<" : "<<ss.str()<<endl;
	  }



	  AtomPointerVector test;
	  for (uint j = 0; j < refCa.size();j++){
	    test.push_back(fragDB[i+j]);
	  }
	  test.saveCoor("pre");

	  if (!t.rmsdAlignment(test,refCa)){
	    test.applySavedCoor("pre");
	    continue;
	  }
	  
	  double rmsd = test.rmsd(refCa);	
	  if (rmsd < opt.rmsd){

	    char matchChar[1000];
	    
	    sprintf(matchChar, "MATCH %4s, %1s - %4d - %3s, chain %1s and resi %4d-%-4d, %8.3f",
			   at1.getSegID().c_str(),
			   at1.getChainId().c_str(),at1.getResidueNumber(),at1.getResidueName().c_str(), 
			   at1.getChainId().c_str(),at1.getResidueNumber(),test.back()->getResidueNumber(),
			   rmsd);
	    
	    string matchString = matchChar;


	    if (opt.pdbDir == ""){
	      char tmp[100];
	      sprintf(tmp,"/tmp/hitCa-%06d.pdb",match);
	      stringstream ss;
	      ss << tmp;

	      PDBWriter pout;
	      pout.open(ss.str());
	      pout.write(test);
	      pout.close();

	    } else {

	      stringstream ss;
	      ss << opt.pdbDir <<"/"<<at1.getSegID()<<".pdb";

	      System sys;
	      sys.readPdb(ss.str());
	      for (uint ats = 0; ats < sys.getAtomPointers().size();ats++){
		sys.getAtom(ats).setSegID("");
	      }
	      AtomSelection sel2(sys.getAtomPointers());

	      ss.str("");
	      char tmp[100];
	      sprintf(tmp,"chain %1s and resi %d-%-d and name CA",at1.getChainId().c_str(),at1.getResidueNumber(),test.back()->getResidueNumber());
	      ss << tmp;
	      AtomPointerVector caAts = sel2.select(ss.str());


	      if (!t.rmsdAlignment(caAts,test,sys.getAtomPointers())){
		MSLOUT.stream() << matchString << endl;
		MSLOUT.stream() << "Problem aligning all atoms using the C-alpha trace"<<endl;
		MSLOUT.stream() << "\tTrying to align with: "<<ss.str()<<endl;
		MSLOUT.stream() << "\tSystem: "<<sys.getSizes()<<endl;
		MSLOUT.stream() << "\tSelected: "<<caAts.size()<<" atoms using '"<<ss.str()<<"'"<<endl;
		continue;
	      } 
	    
	      ss.str("");
	      char tmp2[100];
	      sprintf(tmp2,"hitAll-%06d-%4s.pdb",match,at1.getSegID().c_str());
	      ss << tmp2;


	      PDBWriter pout;
	      pout.open(ss.str());
	      pout.write(sys.getAtomPointers());
	      pout.close();

	      matchString += " "+ss.str();

	      if (opt.printSasa){

		
		// Get Sasa from fragment by itself and in context of its PDB
		ss.str("");
		char tmp[100];
		sprintf(tmp,"frag, chain %1s and resi %d-%-d",at1.getChainId().c_str(),at1.getResidueNumber(),test.back()->getResidueNumber());
		ss << tmp;
		AtomPointerVector fragAllAts = sel2.select(ss.str());

		SasaCalculator sasa(fragAllAts);
		sasa.calcSasa();
		double fragSasa = sasa.getTotalSasa();


		AtomPointerVector fragEnvAllAts = sel2.select("all WITHIN 5 OF frag");
		SasaCalculator sasa2(fragEnvAllAts);
		sasa2.calcSasa();
		double envSasa = sasa2.getTotalSasa();

		double tertSasa = envSasa - fragSasa;

		char sasaStr[30];
		sprintf(sasaStr," %8.3f",tertSasa);
		
		matchString += (string)sasaStr;
		
	      }
	      
	      // Filter for matches that have specific residues buried
//	      if (opt.tertSASA != ""){
//		ss.str("");
//		char tmp[100];
//		sprintf(tmp,"chain %1s and resi %d-%-d",at1.getChainId().c_str(),at1.getResidueNumber(),test.back()->getResidueNumber());
//		ss << tmp;
//		
//		AtomPointerVector allAts = sel2.select(ss.str());
//	      }


	    } // I PDBDIR NOT DEFINED


	    MSLOUT.stream() << matchString<<endl;
	    match++;

	    if (match > opt.maxMatches) {
	      	cout <<"Done due to maxMatches. took: "<<(tm.getWallTime() - start)<<" seconds."<<endl<<endl;
		exit(0);
	    }

	  } // IF RMSD < TOLERANCE

	  test.applySavedCoor("pre");
	  
	  
	}

	cout <<"Done. took: "<<(tm.getWallTime() - start)<<" seconds."<<endl<<endl;
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
		cout << "searchFragmentDatabase --database FRAG_DB \n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.database = OP.getString("database");
	if (OP.fail()){
		cerr << "ERROR 1111 database not specified.\n";
		exit(1111);
	}
	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  cout << "WARNING using RMSD of 2.0"<<endl;
	  opt.rmsd = 2.0;
	}
	opt.pdbDir = OP.getString("pdbDir");
	if (OP.fail()){
	        cout << "WARNING: pdbDir not set, this means only Ca pdbs will print out"<<endl;
		opt.pdbDir = "";
	}

	opt.regex = OP.getString("regex");
	if (OP.fail()){
	  opt.regex ="";
	}

	opt.printSasa = OP.getBool("printSasa");
	if (OP.fail()){
	  opt.printSasa = false;
	}
	opt.maxMatches = OP.getInt("maxMatches");
	if (OP.fail()){
	  opt.maxMatches = 1000;
	  cerr << "WARNING maxMatches set to "<<opt.maxMatches<<endl;
	}
	return opt;
}



