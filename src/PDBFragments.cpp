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


#include "PDBFragments.h"
#include "BBQTable.h"
#include "Transforms.h"
#include "AtomSelection.h"
#include "MslOut.h"

// BOOST Includes
#include <boost/regex.hpp>


static MslOut MSLOUT("PDBFragments");
using namespace MSL;
using namespace std;


/*
Input:    
    _sys                    :   MSL System
    _stemResidues           :   id of "stem" residues inside chain   (current needs to be size of 4, 2 on each side of the locally sampling fragment)
    _numResiduesInFragment :*  number of residues wanted between "stems" (default is calculated from firstOfLastTwoStems.getResidueNumber() - lastOfFirstTwoStems.getResidueNumber() )

Output:
    int number of fragments found.
 */
int PDBFragments::searchForMatchingFragments(System &_sys, vector<string> &_stemResidues,int _numResiduesInFragment, string _regex){

	if (fragDB.size() <= 4){
		cerr << "ERROR 1232 PDBFragments::searchForMatchingFragments(), fragment database is too small, maybe you forgot to load it using loadFragmentDatabase() ?"<<endl;
		exit(1232);
	}
	int numFrags = 0;
	int illegalQuads = 0;

	// Remove last set of results
	lastResults.clear();

	// Clear our set of matched sequences
	matchedSequences.clear();


	// If stem residues is 0, the use all the atoms of _ats to search
	if (_stemResidues.size() == 2){
	} else {

		// Get two stem vectors
		if (_stemResidues.size() % 2 != 0){
			cerr << "ERROR 1967 PDBFragments::searchForMatchingFragments() , number of stem residues is not even\n";
			exit(1967);
		}


		// Create 2 lists (AtomPointerVectors) of the cterm and nterm stem.
		int numStemResidues =  _stemResidues.size() / 2;
		AtomPointerVector stem1;
		AtomPointerVector stem2;
		for (uint i = 0; i < _stemResidues.size();i++){
			if (i < numStemResidues){
				stem1.push_back(&_sys.getPosition(_stemResidues[i]).getAtom("CA"));
			} else {
				stem2.push_back(&_sys.getPosition(_stemResidues[i]).getAtom("CA"));
			}

		}

		AtomPointerVector stemBBats;
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("N"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("CA"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("C"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("N"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("CA"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("C"));

		//MSLOUT.stream() << "STEM sizes: "<<stem1.size()<<","<<stem2.size()<<endl;

		// Check for small stems
		if (stem1.size() <= 1 || stem2.size() <= 1){
			cerr<< "PDBFragments::searchForMatchingFragments() . Stems are too small."<<stem1.size()<<","<<stem2.size()<<endl;
		}

		// Use natural sequence length if num residues in fragment not set
		if (_numResiduesInFragment == -1){
			_numResiduesInFragment = stem2(0).getResidueNumber() - stem1(stem1.size()-1).getResidueNumber() - 1;
		}

		MSLOUT.stream() << "Number of residues between stems: "<<_numResiduesInFragment<<endl;

		// Store stem-to-stem distance-squared vector (for filtering candidates below)
		AtomPointerVector stems = stem1 + stem2;
		vector<double> stemDistanceSq;
		for (uint c = 0; c < stem1.size();c++){
			  fprintf(stdout,"C: %1s %4d %3s\n",
			  stem1[c]->getChainId().c_str(),
			  stem1[c]->getResidueNumber(),
			  stem1[c]->getResidueName().c_str());
			for (uint n = 0; n < stem2.size();n++){

				double distSq = stem1(c).distance2(stem2(n));			
				  fprintf(stdout,"\tN: %1s %4d %3s = %8.3f\n",
				  stem2[n]->getChainId().c_str(),
				  stem2[n]->getResidueNumber(),
				  stem2[n]->getResidueName().c_str(),distSq);

				stemDistanceSq.push_back(distSq);
			}
		}


		// BBQTable is used for adding backbone-atoms to a c-alpha trace.
		BBQTable bbqT;
		//bbqT.setDebugFlag(true);
		if (fragType == caOnly){

			// BBQ Table for adding backbone atoms
			bbqT.openReader(bbqTable);
		}

		MSLOUT.stream() << "FragDB.size(): "<<fragDB.size();
		// Now loop over all fragments in database checking for ones of the correct size
		double tol = 16; // Tolerance of distance to be deviant from stems in Angstroms^2
		int matchIndex = 0; // Index for keeping track of matches
		for (uint i = 0 ; i < fragDB.size()-(_numResiduesInFragment+stem1.size()+stem2.size());i++){

			// Get proposed ctermStem
			AtomPointerVector ctermStem;
			for (uint n = 0; n < stem1.size();n++){
				ctermStem.push_back(fragDB[i+n]);
			}

			// Get proposed ntermStem
			AtomPointerVector ntermStem;
			for (uint n = 0; n < stem2.size();n++){
				ntermStem.push_back(fragDB[i+stem1.size()+_numResiduesInFragment+n]);
			}

			// Check for same PDB
			if (ctermStem(0).getSegID() != ntermStem(ntermStem.size()-1).getSegID()){
			        //fprintf(stdout,"New PDB: %4s to %4s\n",ctermStem[0]->getSegID().c_str(),ntermStem(ntermStem.size()-1).getSegID().c_str());

				// Jump ahead to move past all pdb1 v pdb2 tests
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}


			// Check for same Chain
			if (ctermStem(0).getChainId() != ntermStem(ntermStem.size()-1).getChainId()){
				/*
				fprintf(stdout,"New Chain: %4s %1s to %4s %1s\n",
					ctermStem[0]->getSegID().c_str(),
					ctermStem[0]->getChainId().c_str(),					
					ntermStem(ntermStem.size()-1).getSegID().c_str(),
					ntermStem(ntermStem.size()-1).getChainId().c_str());
				*/
				// Jump ahead to move past all chain1 v chain2 tests
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}


			// Check for a gap (PDB Structures have gaps).
			if ( abs(ctermStem[ctermStem.size()-1]->getResidueNumber() - ntermStem[0]->getResidueNumber()) != _numResiduesInFragment+1){
				/*
				fprintf(stdout," Gap found at residue %4s %1s %4d %3s and %4s %1s %4d %3s\n",
					ctermStem[ctermStem.size()-1]->getSegID().c_str(),
					ctermStem[ctermStem.size()-1]->getChainId().c_str(),
					ctermStem[ctermStem.size()-1]->getResidueNumber(),
					ctermStem[ctermStem.size()-1]->getResidueName().c_str(),
					ntermStem[0]->getSegID().c_str(),
					ntermStem[0]->getChainId().c_str(),
					ntermStem[0]->getResidueNumber(),
					ntermStem[0]->getResidueName().c_str());
				*/
				// Jump ahead to move past gap
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}
			

			// Distance Filter... are the proposed stems close enough ?
			bool passDistanceFilter = true;
			int index = 0;

			// Check against compiled list of distance-squared's
			for (uint c = 0; c < ctermStem.size();c++){
				/*
				  fprintf(stdout," Checking residue %1s %4d %3s\n",
				  ctermStem[c]->getChainId().c_str(),
				  ctermStem[c]->getResidueNumber(),
				  ctermStem[c]->getResidueName().c_str());
				*/
				for (uint n = 0; n < ntermStem.size();n++){

					double distSq = ctermStem[c]->distance2(*ntermStem[n]);

					/*
					  fprintf(stdout," \t%1s %4d %3s = %8.3f vs %8.3f\n",
					  ntermStem[n]->getChainId().c_str(),
					  ntermStem[n]->getResidueNumber(),
					  ntermStem[n]->getResidueName().c_str(),distSq,stemDistanceSq[index]);
					*/
					if (abs(stemDistanceSq[index++] - distSq) > tol){
						passDistanceFilter = false;
						break;
					}
				}

				if (!passDistanceFilter){
					break;
				}
			}



			// Continue if the distance filter not passed
			if (!passDistanceFilter) {
				continue;
			}
			MSLOUT.stream() << "PASSED DISTANCE FILTER!"<<endl;

			// align and print winning fragment

			// Make a copy of frag stem, so not to effect original atoms
			AtomPointerVector fragStem;
			for (uint ct = 0; ct < ctermStem.size();ct++){
				fragStem.push_back(new Atom(ctermStem(ct)));
			}
			for (uint nt = 0; nt < ntermStem.size();nt++){
				fragStem.push_back(new Atom(ntermStem(nt)));
			}

			// Align the fragStem to stem to see if it matches well enough...
			Transforms tm;
			fragStem.saveCoor("pre");
			tm.rmsdAlignment(fragStem,stems);
			
			double rmsd = fragStem.rmsd(stems);
			MSLOUT.stream() << "RMSD: "<<rmsd<<endl;
			// Continue if the RMSD filter not passed
			if (rmsd > 0.5){
				continue;
			}




			matchIndex++;
			fragStem.applySavedCoor("pre");

			//lastResults->writePdb("/tmp/preAdd.pdb");

			// Create new atoms for middle atoms..
			AtomPointerVector tmp;
			index = 1;
			string matchSeq = "";
			for (uint f = 0; f < ctermStem.size();f++){
			        //fragStem(f).setResidueNumber(index++);
				//fragStem(f).setResidueName("FRG");
				//fragStem(f).setChainId("A");
				//fragStem(f).setSegID("");
				tmp.push_back(&fragStem(f));
				//cout << "C-ADD: "<<fragStem(f)<<" "<<fragStem(f).getSegID()<<endl;
				matchSeq += MslTools::getOneLetterCode(fragStem(f).getResidueName());
			}


			for (int f = 0; f < _numResiduesInFragment;f++){
			        Atom *a =new Atom(fragDB(i+ctermStem.size()+f));
				//a->setCoor(fragDB(i+ctermStem.size()+f).getCoor());
				//a->setResidueNumber(index++);
				//a->setResidueName("FRG");
				//a->setChainId("A");
				//a->setSegID("");
				tmp.push_back(a);
				//cout << "F-ADD: "<<*a<<" "<<a->getSegID()<<endl;

				matchSeq += MslTools::getOneLetterCode(fragDB(i+ctermStem.size()+f).getResidueName());
			} 		


			for (uint f = 0; f < ntermStem.size();f++){
				//fragStem(f+ctermStem.size()).setResidueNumber(index++);
				//fragStem(f+ctermStem.size()).setResidueName("FRG");
				//fragStem(f+ctermStem.size()).setChainId("A");
				//fragStem(f+ctermStem.size()).setSegID("");
				tmp.push_back(&fragStem(f+ctermStem.size()));
				matchSeq += MslTools::getOneLetterCode(fragStem(f+ctermStem.size()).getResidueName());
				//cout << "N-ADD: "<<fragStem(f+ctermStem.size())<<" "<<fragStem(f+ctermStem.size()).getSegID()<<endl;
			}
			string key = MslTools::stringf("%06d-%s-%1s_%04d%s-%1s_%04d%s",
						       matchIndex,
						       fragDB(i+ctermStem.size()).getSegID().c_str(),
						       fragDB(i+ctermStem.size()).getChainId().c_str(),
						       fragDB(i+ctermStem.size()).getResidueNumber(),
						       fragDB(i+ctermStem.size()).getResidueIcode().c_str(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getChainId().c_str(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getResidueNumber(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getResidueIcode().c_str());

			matchedSequences[key] = matchSeq;

			if (_regex != ""){


			  if (!boost::regex_search(matchSeq.c_str(),boost::regex(_regex))){
			    continue;
			  } else {
			    MSLOUT.stream() << "RegEx Matched."<<endl;
			  }

			}
			

			Chain tmpChain;
			tmpChain.addAtoms(tmp);			

			
			tm.rmsdAlignment(fragStem,stems,tmpChain.getAtomPointers());


			fprintf(stdout,"(%4s and chain %1s and resi %3d-%3d)  %8.3f ",
				ctermStem[0]->getSegID().c_str(),
				ctermStem[0]->getChainId().c_str(),
				ctermStem[0]->getResidueNumber(),
				ntermStem[ntermStem.size()-1]->getResidueNumber(),
			        rmsd);

			
			bool successful = true;
			if (fragType == caOnly){
			  if (pdbDir != ""){

			              Atom &at1 = tmpChain.getAtom(0);
			              Atom &at2 = tmpChain.getAtom(tmpChain.atomSize()-1);

				      string allAtomFileName = MslTools::stringf("%s/%s.pdb",pdbDir.c_str(),at1.getSegID().c_str());

				      System allAtomSys;
				      allAtomSys.readPdb(allAtomFileName);
				      for (uint ats = 0; ats < allAtomSys.getAtomPointers().size();ats++){
					allAtomSys.getAtom(ats).setSegID("");
				      }

				      AtomSelection sel2(allAtomSys.getAtomPointers());

			    	      stringstream ss;
				      char tmpstr[100];
				      sprintf(tmpstr,"chain %1s and resi %d-%-d and name CA",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
				      ss << tmpstr;
				      AtomPointerVector caAts = sel2.select(ss.str());


				      if (!tm.rmsdAlignment(caAts,tmpChain.getAtomPointers(),allAtomSys.getAtomPointers())){
					MSLOUT.stream() << "Problem aligning all atoms using the C-alpha trace"<<endl;
					MSLOUT.stream() << "PDB: "<<allAtomFileName<<endl;
					MSLOUT.stream() << "\tTrying to align with: "<<ss.str()<<endl;
					MSLOUT.stream() << "\tSystem: "<<allAtomSys.getSizes();
					MSLOUT.stream() << "\tSelected: "<<caAts.size()<<" atoms using '"<<ss.str()<<"'"<<endl;
					MSLOUT.stream() << "\tReference: "<<tmpChain.getAtomPointers().size()<<" atoms"<<endl;
					MSLOUT.stream() << tmpChain.getAtomPointers();
					continue;
				      } 

				      // Get BB atoms from stem-equivalent residues
				      AtomPointerVector allAts_stemBBats;
				      Position &stem1eq  = allAtomSys.getPosition(MslTools::getPositionId(at1.getChainId(),at1.getResidueNumber(),at1.getResidueIcode()));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("C"));

				      Position &stem2eq  = allAtomSys.getPosition(MslTools::getPositionId(at2.getChainId(),at2.getResidueNumber()+stem2.size()-1,at2.getResidueIcode()));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("C"));
				      	
				      double bbRMSD = allAts_stemBBats.rmsd(stemBBats);

				      fprintf(stdout, "%8.3f",bbRMSD);
				      if (bbRMSD > 1.5){
					continue;
				      }

				      ss.str("");
				      char tmpstr2[100];
				      sprintf(tmpstr2,"chain %1s and resi %d-%-d",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
				      ss << tmpstr2;
				      AtomPointerVector allAts = sel2.select(ss.str());
				      lastResults.push_back(new AtomContainer());
				      lastResults.back()->addAtoms(allAts);




			  } else if (bbqTable != ""){
			  
				illegalQuads += bbqT.fillInMissingBBAtoms(tmpChain);

				if (illegalQuads > 0){
				  successful = false;
			  	  exit(0);
				}

				lastResults.push_back(new AtomContainer(tmpChain.getAtomPointers()));
			  }
				
			} else {
			  MSLOUT.stream() <<" ALL ATOMS "<<endl;
			}
			fprintf(stdout,"\n");

			fragStem.applySavedCoor("pre");

			if (successful){
				numFrags++;
			}

			/*					
			stringstream ss;
			PDBWriter pout;
			pout.open(ss);
			pout.write(sys.getAtomPointers());
			pout.close();

			matchingFragments.push(pair<double,string>(rmsd,ss.str()));
			*/


		}

		//fprintf(stdout, "Best RMSD: %8.3f\n",matchingFragments.top().first);


		fprintf(stdout,"Number of succesful fragments: %10d, illegal quads: %10d\n",numFrags,illegalQuads);

	}
	
	
	return numFrags;
}


