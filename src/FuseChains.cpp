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

#include "FuseChains.h"
#include "AtomSelection.h"

using namespace MSL;



FuseChains::FuseChains(){
}

FuseChains::FuseChains(const FuseChains &_fuse){
  
}

FuseChains::~FuseChains(){
}




AtomPointerVector & FuseChains::fuseInsert(Chain &_template, Chain &_insert, string templateStem1posId, string templateStem2posId){

	// Clear previous fused chains
	fusedChains.removeAllAtoms();

     
     	// Discover stems by distance  [ STEMS = INSERTION POINTS ]
	int stem1index = -1;
	int stem2index = -1;
	if (templateStem1posId == "" || templateStem2posId == ""){

	  // Get the C-alpha's in order on the shortest chain (otherwise you can have HOH residues causing problems here)
	  AtomSelection sel(_insert.getAtomPointers());
	  AtomPointerVector frag = sel.select("name CA");

	  // Search for a template residue that is on top of the first AminoAcid in the fragment chain
	  for (uint j = 0; j < _template.size();j++){
	    Residue &template1 = _template.getPosition(j).getCurrentIdentity();
	    if (template1.atomExists("CA") && frag[0]->distance(template1.getLastFoundAtom()) < 0.3){
	            stem1index = j;
		    break;
	      }
	  }

	  // Search for a template residue that is on top of the last AminoAcid in the fragment chain
	  for (uint j = 0; j < _template.size();j++){
	    Residue &template1 = _template.getPosition(j).getCurrentIdentity();
	    if (template1.atomExists("CA") && frag[frag.size()-1]->distance(template1.getLastFoundAtom()) < 0.3){
	            stem2index = j;
		    break;
	      }
	  }
	  frag.clear();
	  
	} else {

	  // If stems are specified on command-line then just use those positionId's.
	  stem1index = _template.getPosition(templateStem1posId).getIndexInChain();
	  stem2index = _template.getPosition(templateStem2posId).getIndexInChain();
	}

	// ERROR Checking that stems have been Identified
	if (stem1index == -1 || stem2index == -1){
	  cerr << "ERROR 3333 Stem residues in template chain  are undefined\n";
	  exit(3333);
	}

	// Now create the stem positions
	Position &stem1 = _template.getPosition(stem1index);
	Position &stem2 = _template.getPosition(stem2index);

	// Let the user know what the stems are..
	fprintf(stdout, "Stem residues: %s , %s\n",stem1.toString().c_str(),stem2.toString().c_str());

	
	cout << "Add template chain upto stem1"<<endl;
	for (uint i = 0; i < stem1.getIndexInChain();i++){
	  fusedChains.addAtoms(_template.getPosition(i).getAtomPointers());
	}


	int newResidueNumber = stem1.getResidueNumber();
	bool stem1added = false;
	for (uint i = 0; i < _insert.size();i++){
	  Residue &fragRes = _insert.getPosition(i).getCurrentIdentity();
	  if (!fragRes.atomExists("CA")) continue;
	  if (fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem1.getCurrentIdentity()("CA")) < 0.3) {
	       fprintf(stdout, "Add insert residue(%s) as template stem1 residue\n",fragRes.toString().c_str());

	       AtomContainer newStemRes;
	       newStemRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newStemRes.size();a++){
		 newStemRes[a].setChainId(stem1.getChainId());
		 newStemRes[a].setResidueNumber(newResidueNumber);
		 newStemRes[a].setResidueIcode(stem1.getResidueIcode());
	       }
	       
	       fusedChains.addAtoms(newStemRes.getAtomPointers());

	       stem1added = true;
	       
	       newResidueNumber++;
	       continue;
	  }



	  if (stem1added && fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem2.getCurrentIdentity()("CA")) < 0.3) {
	       fprintf(stdout, "Add insert residue(%s) as template stem2 residue\n",fragRes.toString().c_str());

	       AtomContainer newStemRes;
	       newStemRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newStemRes.size();a++){
		 newStemRes[a].setChainId(stem2.getChainId());
		 newStemRes[a].setResidueNumber(newResidueNumber);
	       }
	       
	       fusedChains.addAtoms(newStemRes.getAtomPointers());
	       newResidueNumber++;
	       // We should be done by now..
	       break;
	  }
	  
	  if (stem1added) {

	      
	       AtomContainer newFragRes;
	       newFragRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newFragRes.size();a++){
		 newFragRes[a].setChainId(stem1.getChainId());
		 newFragRes[a].setResidueNumber(newResidueNumber);
	       }

	       fprintf(stdout, "Add insert residue(%s) as newly inserted residue(%s,%4d)\n",fragRes.toString().c_str(),stem1.getChainId().c_str(),newResidueNumber);	       

	       fusedChains.addAtoms(newFragRes.getAtomPointers());
	       newResidueNumber++;
	  }

	}


	// Add the rest of the template protein
	cout << "Add rest of template chain"<<endl;
	for (uint i = stem2.getIndexInChain()+1; i < _template.size();i++){
	  AtomContainer tmp;
	  tmp.addAtoms(_template.getPosition(i).getAtomPointers());
	  for (uint t = 0; t< tmp.size();t++){
	    tmp[t].setResidueNumber(newResidueNumber);
	  }
	  fusedChains.addAtoms(tmp.getAtomPointers());
	  newResidueNumber++;
	}


	return fusedChains.getAtomPointers();

}
