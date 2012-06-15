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
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "AtomContainer.h"
#include "FuseChains.h"
#include "insertLoopIntoTemplate.h"

#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace MSL;
using namespace std;
using namespace MslTools;

int main(int argc, char *argv[]) {

        // Parse command line options into Options structure
	Options opt = setupOptions(argc, argv);

	// Template is the structure to insert into
	System templateInput;
	templateInput.readPdb(opt.templatePDB);
	AtomPointerVector ats = templateInput.getAtomPointers();
	if (opt.templateChain != ""){
	  if (templateInput.chainExists(opt.templateChain)){
	    ats = templateInput.getChain(opt.templateChain).getAtomPointers();
	  } else {
	    cerr << "ERROR 333 templatePDB does not have chain ("<<opt.templateChain<<")"<<endl;
	    exit(3333);
	  }
	}

	System templatePDB;
	templatePDB.addAtoms(ats);
	if (templatePDB.chainSize() > 1 ){
	  cerr << "ERROR 3333 currently insertLoopIntoTemplate assumes templatePDB has a single chain, I found "<<templatePDB.chainSize()<<" chains."<<endl;
	  exit(3333);
	}

	cout << "Template chain "<<templatePDB.getChain(0).getChainId()<<" has "<<templatePDB.getChain(0).positionSize()<< " residues."<<endl;
 
	// Fragment is the structure to take peice of structure from
	System fragmentPDB;
	fragmentPDB.readPdb(opt.fragmentPDB);

	// Find the shortest chain to use as a fragment
	int shortestChain = 0;
	for (uint c = 1; c < fragmentPDB.chainSize();c++){
	  if (opt.fragmentChain == fragmentPDB.getChain(c).getChainId()){
	    shortestChain = c;
	    break;
	  } 

	  if (fragmentPDB.getChain(c).positionSize() < fragmentPDB.getChain(shortestChain).positionSize()){
	    shortestChain = c;
	  }
	}

	Chain &fragChain = fragmentPDB.getChain(shortestChain);
	cout << "Fragment Chain is "<<fragChain.getChainId()<<" from "<<opt.fragmentPDB<< " "<<fragChain.positionSize()<< " residues"<<endl;

	AtomContainer fusedProtein;
	FuseChains fuse;
	fusedProtein.addAtoms(fuse.fuseInsert(templatePDB.getChain(0), fragChain,opt.templateStem1,opt.templateStem2,opt.includeTemplateStems));
	

	if (opt.checkCaCaDistances){

	  bool chainOk = true;
	  int chainBreak1 = 0;
	  int chainBreak2 = 0;
	  for (uint i = 0; i < fusedProtein.size()-1;i++){
	    if (fusedProtein[i].getName() != "CA") continue;

	    int nextCa = 0;
	    for (uint j = i+1; j < fusedProtein.size()-1;j++){
	      if (fusedProtein[j].getName() == "CA") { nextCa = j; break; }
	    }

	    if (fusedProtein[i].distance(fusedProtein[nextCa]) > 4.0){
		chainOk = false;
		chainBreak1 = i;
		chainBreak2 = nextCa;
		break;
	      }
	  }

	  if (!chainOk){
	    fprintf(stdout, "ERROR fusion of %s with %s has chain break at %6s -> %6s, no file output\n",opt.templatePDB.c_str(),opt.fragmentPDB.c_str(),fusedProtein[chainBreak1].getPositionId().c_str(),fusedProtein[chainBreak2].getPositionId().c_str());
	    exit(2543);
	  }

	}

	// Add additional chains from the fragmentPDB
	int chainIdIndex = 0;
	for (uint i = 0; i < fragmentPDB.chainSize();i++){
	    if (i == shortestChain) continue;

	    string chains = "BCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	    string newChainId = chains.substr(chainIdIndex,1);	
	    while (fragmentPDB.chainExists(newChainId)){

	        chainIdIndex++;
	        newChainId = chains.substr(chainIdIndex,1);	

		if (chainIdIndex > chains.length()){
		  cerr << "ERROR 5555 insertLoopIntoTemplate couldn't find a free chain id to use for fragmentPDB chain "<< fragmentPDB.getChain(i).getChainId()<<endl;
		  exit(5555);
		}
	    } 

	    fprintf(stdout, "Adding chain %1s from fragmentPDB as chain %s\n",fragmentPDB.getChain(i).getChainId().c_str(),newChainId.c_str());

	    // Add each atom and change its chain id.
	    for (uint a = 0; a < fragmentPDB.getChain(i).getAtomPointers().size();a++){
		      fusedProtein.addAtom(fragmentPDB.getChain(i).getAtom(a));
		      fusedProtein.getAtom(fusedProtein.size()-1).setChainId(newChainId);
	    }

	    // Increment the chain index
	    chainIdIndex++;
	}

	int numClashes = 0;
	if (opt.clashCheck){
	  for (uint i = 0; i < fusedProtein.size();i++){
	    if (fusedProtein[i].getName() != "CA") continue;

	    for (uint j = i+1;j < fusedProtein.size();j++){
	      if (fusedProtein[j].getName() != "CA") continue;

	      if (fusedProtein[i].distance(fusedProtein[j]) < 2.5){
		numClashes++;
	      }

	    }
	  }

	}

	char fname_base[200];
	sprintf(fname_base, "%s_%s",MslTools::getFileName(opt.templatePDB).c_str(),MslTools::getFileName(opt.fragmentPDB).c_str());
	int numClashTolerance = opt.numClashes;
	if (numClashes <= opt.numClashes) {
	  char fname[200];
	  sprintf(fname, "%s.pdb",fname_base);
	  cout << "Fused Protein: "<< fname << " number of clashes: "<<numClashes<<endl;
	  fusedProtein.writePdb(fname);
	} else {
	  fprintf(stdout, " ERROR fusion of %s with %s has %d clashes which is more than the tolerance of %d clashes, no file output\n",opt.templatePDB.c_str(),opt.fragmentPDB.c_str(),numClashes,numClashTolerance);
	}

	if (opt.outputRosettaFiles != ""){
	  cout << "Output rosetta files: "<<opt.outputRosettaFiles<<endl;
	  System fusedSystem;
	  fusedSystem.addAtoms(fusedProtein.getAtomPointers());

	  vector<string> insertedAndNeighboringPositions;
	  vector<string> insertedPositions = fuse.getInsertedPositions();
	  cout << "Number of inserted positions: "<<insertedPositions.size()<<endl;
	  for (uint i = 0; i < insertedPositions.size();i++){
	    Position &pos = fusedSystem.getPosition(insertedPositions[i]);
	    insertedAndNeighboringPositions.push_back(insertedPositions[i]);

	    // Get Neighbors
	    vector<int> neighbors = pos.getCurrentIdentity().findNeighbors(4.0);
	    for (uint j =0;j < neighbors.size();j++){
	      insertedAndNeighboringPositions.push_back(fusedSystem.getPosition(neighbors[j]).getPositionId());
	    }
	  }
	  
	  cout << "Number of inserted positions + neighbors: "<<insertedAndNeighboringPositions.size()<<endl;
	  // Write out Resfile.ALA and Resfile.Design
	  fstream resfileA;
	  resfileA.open(MslTools::stringf("%s_resfile.ALA",opt.outputRosettaFiles.c_str()).c_str(),std::ios::out);
	  resfileA << "NATRO"<<endl<<"start"<<endl;
	  fstream resfileD;
	  resfileD.open(MslTools::stringf("%s_resfile.DES",opt.outputRosettaFiles.c_str()).c_str(),std::ios::out);
	  resfileD << "NATRO"<<endl<<"start"<<endl;
	  for (uint i = 0; i < insertedAndNeighboringPositions.size();i++){
	    Position &pos = fusedSystem.getPosition(insertedAndNeighboringPositions[i]);
	    string lineALA = MslTools::stringf("%d%1s %1s PIKAA A",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    if (pos.getResidueName() == "GLY"){
	      lineALA = MslTools::stringf("%d%1s %1s PIKAA G",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    }
	    if (pos.getResidueName() == "PRO"){
	      lineALA = MslTools::stringf("%d%1s %1s PIKAA P",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    }

	    resfileA << lineALA << endl;

	    string lineDES = MslTools::stringf("%d%1s %1s ALLAA",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    if (pos.getResidueName() == "GLY"){
	      lineDES = MslTools::stringf("%d%1s %1s PIKAA G",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    }
	    if (pos.getResidueName() == "PRO"){
	      lineDES = MslTools::stringf("%d%1s %1s PIKAA P",pos.getResidueNumber(),pos.getResidueIcode().c_str(),pos.getChainId().c_str());
	    }
	    resfileD << lineDES <<endl;
	  }
	  resfileA.close();
	  resfileD.close();

	  // Write out insertedLoops.txt
	  fstream posFile;
	  posFile.open(MslTools::stringf("%s_insertedPositions.txt",opt.outputRosettaFiles.c_str()).c_str(), std::ios::out);
	  for (uint i = 0; i < insertedPositions.size();i++){
	    Position &pos = fusedSystem.getPosition(insertedPositions[i]);
	    posFile << MslTools::stringf("%s%d%s ",pos.getChainId().c_str(), pos.getResidueNumber(),pos.getResidueIcode().c_str());
	  }
	  posFile << endl;
	  posFile.close();
	}
	
								     

}
	

Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;

	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage: insertSelectionIntoTemplate " << endl;
		cout << endl;
		cout << "\n";
		cout << "template PDB\n";
		cout << "templateStem1 positionId\n";
		cout << "templateStem2 positionId\n";
		cout << "fragment PDB\n";
		cout << endl;
		exit(0);
	}

	opt.templatePDB = OP.getString("template");
	if (OP.fail()){
		cerr << "ERROR 1111 no template specified."<<endl;
		exit(1111);
	}
	opt.templateChain = OP.getString("templateChain");
	if (OP.fail()){
	  opt.templateChain = "";
	}

	opt.fragmentPDB = OP.getString("fragment");
	if (OP.fail()){
		cerr << "ERROR 1111 no fragment specified."<<endl;
		exit(1111);
	}
	opt.fragmentChain = OP.getString("fragmentChain");
	if (OP.fail()){
	  opt.fragmentChain = "NO_CHAIN_INPUT_USE_SHORTEST_CHAIN_IN_FRAGMENT_PDB";
	}

	opt.templateStem1= OP.getString("templateStem1");
	if (OP.fail()){
	  opt.templateStem1 = "";
	}
	opt.templateStem2 = OP.getString("templateStem2");
	if (OP.fail()){
	  opt.templateStem2 = "";
	}
	opt.includeTemplateStems = OP.getBool("includeTemplateStems");

	opt.clashCheck = OP.getBool("clashCheck");
	opt.numClashes = OP.getInt("numClashes");	
	opt.checkCaCaDistances = OP.getBool("checkCaCaDistances");


	opt.outputRosettaFiles = OP.getString("outputRosettaFiles");
	if (OP.fail()){
	  opt.outputRosettaFiles = "";
	}
	
	return opt;
}
