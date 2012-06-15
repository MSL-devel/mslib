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

#include <string>
#include <algorithm>

#include "getSurroundingResidues.h"
#include "OptionParser.h"
#include "System.h"
#include "ResidueSelection.h"
#include "MslTools.h"
#include "Transforms.h"
#include "AtomContainer.h"

using namespace std;
using namespace MSL;


int main(int argc, char *argv[]){


	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Read-in list of PDBS
	vector<string> listPDBs;
	if (!MslTools::readTextFile(listPDBs,opt.pdblist)){
	  cerr << "ERROR 2222 in reading text file "<<opt.pdblist<<endl;
	  exit(2222);
	}

	AtomContainer refResidueAtoms;
	bool alignResidue = false;
	map<string,int> residueCount;
	map<string,PDBWriter *> pdbWriters;

	char resSel[100];
	sprintf(resSel,"resn %s", opt.residue.c_str());

	for (uint i = 0; i < listPDBs.size();i++){
	  fprintf(stdout,"Working on %s\n",MslTools::getFileName(listPDBs[i]).c_str());
	  System sys;
	  sys.readPdb(listPDBs[i]);


	  ResidueSelection sel(sys);
	  vector<Residue *> focusOnResidues = sel.select(resSel);

	  fprintf(stdout,"\tNumber of residues(%3s) = %10d out of %10d total residues\n",opt.residue.c_str(),(unsigned int)focusOnResidues.size(),sys.positionSize());

	  // For each residue
	  for (uint j = 0; j < focusOnResidues.size();j++){
	    Residue *res = focusOnResidues[j];

	    // For each atom type
	    vector<int> allResidueIndices;
	    for (uint a = 0; a < opt.searchCenterAtoms.size();a++){
	      vector<int> resIndices = res->findNeighbors(opt.distance,opt.searchCenterAtoms[a]);

	      allResidueIndices.insert(allResidueIndices.begin(),resIndices.begin(),resIndices.end());
	    }



	    // Store the end of unique residue indices
	    sort(allResidueIndices.begin(),allResidueIndices.end());
	    vector<int>::iterator new_end = unique(allResidueIndices.begin(),allResidueIndices.end());


	    // Skip rest if no surrounding residues were found
	    if (allResidueIndices.size() == 0) continue;

	    // Transform atom vector by currentResidue -> firstResidue  (using the opt.alignByAtoms)
	    if (alignResidue){

	      
	            // Create a moveable atom vetor
	            AtomPointerVector ats;
		    vector<int>::iterator resIt = allResidueIndices.begin();
		    int numRes = 0;
		    for (;resIt != new_end;resIt++){ 
		      Residue &r = sys.getResidue(*resIt);
		      AtomPointerVector &rAts = r.getAtomPointers();
		      ats.insert(ats.begin(),rAts.begin(),rAts.end());
		      numRes++;
		    }

		    fprintf(stdout, "\t Number of surrouding atoms %5d from %5d residues around %1s,%5d,%s\n", (unsigned int)ats.size(),numRes,res->getChainId().c_str(),res->getResidueNumber(),res->getResidueIcode().c_str());

		    // Add atoms of this residue as well
		    ats.insert(ats.begin(),res->getAtomPointers().begin(),res->getAtomPointers().end());

			  

	            AtomContainer alignAtsToRef;
		    for (uint a = 0; a < opt.alignByAtoms.size();a++){
		      if (res->atomExists(opt.alignByAtoms[a])){
			alignAtsToRef.addAtom(res->getAtom(opt.alignByAtoms[a]));
		      } else {
			cerr << "ERROR 3333 in finding alignByAtoms, atom("<<opt.alignByAtoms[a]<<" in structure "<<MslTools::getFileName(listPDBs[i])<<endl;
			exit(3333);
		      }
		    }

		    
		    // Do alignment
		    Transforms t;
		    t.rmsdAlignment(alignAtsToRef.getAtomPointers(),refResidueAtoms.getAtomPointers(), ats);

		    ats.clear();
		    // Now our surrounding residues are aligned into a reference frame.

	    } else {

    	           // Store residue as reference
	    	    alignResidue = true;
		    
		    for (uint a = 0; a < opt.alignByAtoms.size();a++){
		      if (res->atomExists(opt.alignByAtoms[a])){
			refResidueAtoms.addAtom(res->getAtom(opt.alignByAtoms[a]));
		      } else {
			cerr << "ERROR 3334 in finding alignByAtoms, atom("<<opt.alignByAtoms[a]<<" in structure "<<MslTools::getFileName(listPDBs[i])<<endl;
			exit(3334);
		      }
		    }

		    fprintf(stdout,"\t REF RESIDUE with %2d atoms\n",refResidueAtoms.size());
	    }


	    // For each residue write out a separate PDB file
	    vector<int>::iterator resIt = allResidueIndices.begin();
	    for (;resIt != new_end;resIt++){ 
	      Residue &r = sys.getResidue(*resIt);


	      if (opt.nmrpdb){

		char outName[100];
                if (r.getResidueName().length() == 2) {
 		   sprintf(outName,"%3s_%2s.pdb", res->getResidueName().c_str(), r.getResidueName().c_str());
                }
                else if (r.getResidueName().length() == 1) {
 		   sprintf(outName,"%3s_%1s.pdb", res->getResidueName().c_str(), r.getResidueName().c_str());
                }
                else {
 		   sprintf(outName,"%3s_%3s.pdb", res->getResidueName().c_str(), r.getResidueName().c_str());
                }

		map<string,PDBWriter *>::iterator pdbWritersIt;
		pdbWritersIt = pdbWriters.find(outName);
		if (pdbWritersIt != pdbWriters.end()){
		  AtomPointerVector bothAts;
		  bothAts.insert(bothAts.begin(),r.getAtomPointers().begin(),r.getAtomPointers().end());
		  bothAts.insert(bothAts.begin(),res->getAtomPointers().begin(),res->getAtomPointers().end());

		  pdbWritersIt->second->write(bothAts,true,false,true);
		  bothAts.clear();

		} else {
		  pdbWriters[outName] = new PDBWriter();
		  pdbWriters[outName]->open((string)outName);
		  AtomPointerVector bothAts;
		  bothAts.insert(bothAts.begin(),r.getAtomPointers().begin(),r.getAtomPointers().end());
		  bothAts.insert(bothAts.begin(),res->getAtomPointers().begin(),res->getAtomPointers().end());

		  pdbWriters[outName]->write(bothAts,true,false,true);
		  bothAts.clear();
		}


	      } else {
		char outName[100];
		sprintf(outName,"%3s_%3s_%08d.pdb", res->getResidueName().c_str(), r.getResidueName().c_str(), residueCount[r.getResidueName()]++);
		      
		PDBWriter pout;
		pout.open(outName);		      
		pout.write(res->getAtomPointers());
		pout.write(r.getAtomPointers());
		pout.close();
	      }

	    }

	  }


	  
	  
	}
	
	
	// Close all NMR pdbwriters  SHOULD DELETE them too!
	if (opt.nmrpdb){
	  

	  map<string,PDBWriter *>::iterator it;
	  for (it = pdbWriters.begin(); it != pdbWriters.end();it++){
	    it->second->close();
	  }



	}
	
}


Options setupOptions(int theArgc, char * theArgv[]){

	// Create the options
	Options opt;
	

	// Parse the options
	OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.readArgv(theArgc, theArgv);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: getSurroundingResidues conf" << endl;
		cout << endl;
		cout << "pdblist LIST\n";
		cout << "residue RES_TYPE\n";
		cout << "searchCenterAtoms \"CA CB CG1\"\n";
		cout << "alignByAtoms \"N CA C\"\n";
		cout << endl;
		exit(0);
	}

	opt.pdblist = OP.getString("pdblist");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdblist specified."<<endl;
		exit(1111);
	}

	opt.residue = OP.getString("residue");
	if (OP.fail()){
		cerr << "ERROR 1111 no residue specified."<<endl;
		exit(1111);
	}


        opt.searchCenterAtoms = OP.getStringVectorJoinAll("searchCenterAtoms");
	if (OP.fail()){
	      cerr << "ERROR 1111 no searchCenterAtoms specified."<<endl;
	      exit(1111);
	}

        opt.alignByAtoms       = OP.getStringVectorJoinAll("alignByAtoms");
	if (OP.fail()){

	  // Default to backbone atoms of alpha-aminoacid peptide
	  cerr << "WARNING alignByAtoms not specified so defaulting to N,CA,C"<<endl;
	  opt.alignByAtoms.push_back("N");
	  opt.alignByAtoms.push_back("CA");
	  opt.alignByAtoms.push_back("C");
	}


	opt.distance = OP.getDouble("distance");
	if (OP.fail()){
	  opt.distance = 8;
	  cerr << "WARNING distance not specified using " << opt.distance <<endl;
	}

	opt.nmrpdb = OP.getBool("nmrpdb");
	if (OP.fail()){
	  opt.nmrpdb = false;
	}
	   
	return opt;
}
