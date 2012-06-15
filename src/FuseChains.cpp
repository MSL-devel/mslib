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

#include "FuseChains.h"
#include "AtomSelection.h"
#include "System.h"

using namespace MSL;

#include "MslOut.h"
static MslOut MSLOUT("FuseChains");

FuseChains::FuseChains(){
}

FuseChains::FuseChains(const FuseChains &_fuse){
}

FuseChains::~FuseChains(){
}




AtomPointerVector & FuseChains::fuseInsert(Chain &_template, Chain &_insert, string templateStem1posId, string templateStem2posId, bool _includeTemplateStems){


	// Clear previous fused chains
	fusedChains.removeAllAtoms();
	insertedResidues.clear();
     
     	// Discover stems by distance  [ STEMS = INSERTION POINTS ]
	// Also, make a local insert chain that contains only residues including and between stems..
	System localInsert;

	int stem1index = -1;
	int stem2index = -1;
	if (templateStem1posId == "" || templateStem2posId == ""){

	  int insert1index = -1;
	  int insert2index = -1;
	  // Search for a template residue that is on top of the first AminoAcid in the fragment chain
	  for (uint j = 0; j < _template.positionSize();j++){
	    Residue &template1 = _template.getPosition(j).getCurrentIdentity();

	    for (uint k = 0; k < _insert.positionSize();k++){

	      Residue &insert1 = _insert.getPosition(k).getCurrentIdentity();
	      if (template1.atomExists("CA") && insert1.atomExists("CA") && insert1.getLastFoundAtom().distance(template1.getLastFoundAtom()) < 0.5){
	            stem1index = j;
		    insert1index = k;
		    MSLOUT.stream() << "Template residue "<<template1.getIdentityId()<<" and Insert residue "<<insert1.getIdentityId()<<" are close"<<endl;
		    break;
	      }
	    }
	    if (insert1index != -1) break;
	  }

	  // Search for a template residue that is on top of the last AminoAcid in the fragment chain
	  for (uint j = 0; j < _template.positionSize();j++){
	    Residue &template1 = _template.getPosition(j).getCurrentIdentity();

	    for (uint k = insert1index+2; k < _insert.positionSize();k++){
	      Residue &insert1 = _insert.getPosition(k).getCurrentIdentity();
	      if (template1.atomExists("CA") && insert1.atomExists("CA") && insert1.getLastFoundAtom().distance(template1.getLastFoundAtom()) < 0.5){
	            stem2index = j;
		    insert2index = k;
		    MSLOUT.stream() << "Template residue "<<template1.getIdentityId()<<" and Insert residue "<<insert1.getIdentityId()<<" are close"<<endl;
		    break;
	      }
	    }
	    if (insert2index != -1) break;
	  }

	  for (uint k = 0; k < _insert.positionSize();k++){
	    if (k >= insert1index && k <= insert2index){
	      localInsert.addAtoms(_insert.getPosition(k).getAtomPointers());
	    }
	  }


	} else {

	  // If stems are specified on command-line then just use those positionId's.
	  stem1index = _template.getPosition(templateStem1posId).getIndexInChain();
	  stem2index = _template.getPosition(templateStem2posId).getIndexInChain();

	  localInsert.addAtoms(_insert.getAtomPointers());
	}

	// ERROR Checking that stems have been Identified
	if (stem1index == -1 || stem2index == -1){
	  cerr << "ERROR 3333 Stem residues in template chain are undefined "<<stem1index<<" and "<<stem2index<<endl;
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


	// Add Inserted residues
	int newResidueNumber = stem1.getResidueNumber();
	bool stem1added = false;
	for (uint i = 0; i < localInsert.positionSize();i++){
	  Residue &fragRes = localInsert.getPosition(i).getCurrentIdentity();
	  if (!fragRes.atomExists("CA")) continue;
	  if (fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem1.getCurrentIdentity()("CA")) < 0.5) {

	       if (_includeTemplateStems){
		 fprintf(stdout, "Add use template stem1 residue as template stem1 residue\n");
		 fusedChains.addAtoms(stem1.getAtomPointers());
		 insertedResidues.push_back(stem1.getPositionId());
	       } else {

		 fprintf(stdout, "Add insert residue(%s) as template stem1 residue\n",fragRes.toString().c_str());

		 AtomContainer newStemRes;
		 newStemRes.addAtoms(fragRes.getAtomPointers());
		 for (uint a = 0; a < newStemRes.size();a++){
		   newStemRes[a].setChainId(stem1.getChainId());
		   newStemRes[a].setResidueNumber(newResidueNumber);
		   newStemRes[a].setResidueIcode(stem1.getResidueIcode());
		 }
	       
		 fusedChains.addAtoms(newStemRes.getAtomPointers());

		 insertedResidues.push_back(newStemRes[0].getPositionId());
		 
	       }

	       stem1added = true;
	       
	       newResidueNumber++;
	       continue;
	  }



	  if (stem1added && fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem2.getCurrentIdentity()("CA")) < 0.5) {

	       if (_includeTemplateStems){
		 fprintf(stdout, "Add use template stem2 residue as template stem2 residue\n");
		 for (uint a = 0; a < stem2.atomSize();a++){
		   stem2.getAtom(a).setResidueNumber(newResidueNumber);		   
		 }
		 fusedChains.addAtoms(stem2.getAtomPointers());
		 insertedResidues.push_back(stem2.getPositionId());
	       } else {
		 fprintf(stdout, "Add insert residue(%s) as template stem2 residue\n",fragRes.toString().c_str());

		 AtomContainer newStemRes;
		 newStemRes.addAtoms(fragRes.getAtomPointers());
		 for (uint a = 0; a < newStemRes.size();a++){
		   newStemRes[a].setChainId(stem2.getChainId());
		   newStemRes[a].setResidueNumber(newResidueNumber);
		 }
	       
		 fusedChains.addAtoms(newStemRes.getAtomPointers());
		 insertedResidues.push_back(newStemRes[0].getPositionId());
	       }

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
	       insertedResidues.push_back(newFragRes[0].getPositionId());
	       newResidueNumber++;
	  }

	}

	
	

	// Add the rest of the template protein
	cout << "Add rest of template chain"<<endl;
	for (uint i = stem2.getIndexInChain()+1; i < _template.positionSize();i++){
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
