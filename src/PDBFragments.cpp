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

#include "PDBFragments.h"
#include "BBQTable.h"
#include "Transforms.h"

/*
Input:    
    _ch                    :   MSL Chain
    _stemResidueIndices    :   indices of "stem" residues inside chain   (current needs to be size of 4, 2 on each side of the locally sampling fragment)
    _numResiduesInFragment :*  number of residues wanted between "stems" (default is calculated from firstOfLastTwoStems.getResidueNumber() - lastOfFirstTwoStems.getResidueNumber() )

Output:
    int number of fragments found.
    lastResults (private System*) contains all matching fragments as "altConformations" in each atom.
 */
int PDBFragments::searchForMatchingFragments(Chain &_ch, vector<int> &_stemResidueIndices,int _numResiduesInFragment){

	int numFrags = 0;
	int illegalQuads = 0;

	// Reset the lastResults
	if (lastResults != NULL){
		delete(lastResults);
	}
	lastResults = new System();

	// If stem residues is 0, the use all the atoms of _ats to search
	if (_stemResidueIndices.size() == 2){
	} else {

		// Get two stem vectors
		if (_stemResidueIndices.size() % 2 != 0){
			cerr << "ERROR 1967 PDBFragments::searchForMatchingFragments() , number of stem residues is not even\n";
			exit(1967);
		}


		// Create 2 lists (AtomVectors) of the cterm and nterm stem.
		int numStemResidues =  _stemResidueIndices.size() / 2;
		AtomVector stem1;
		AtomVector stem2;
		for (uint i = 0; i < _stemResidueIndices.size();i++){
			if (i < numStemResidues){
				stem1.push_back(&_ch(_stemResidueIndices[i])("CA"));
			} else {
				stem2.push_back(&_ch(_stemResidueIndices[i])("CA"));
			}

		}

		//cout << "STEM sizes: "<<stem1.size()<<","<<stem2.size()<<endl;

		// Check for small stems
		if (stem1.size() <= 1 || stem2.size() <= 1){
			cerr<< "PDBFragments::searchForMatchingFragments() . Stems are too small."<<stem1.size()<<","<<stem2.size()<<endl;
		}

		// Use natural sequence length if num residues in fragment not set
		if (_numResiduesInFragment == -1){
			_numResiduesInFragment = stem2(0).getResidueNumber() - stem1(stem1.size()-1).getResidueNumber() - 1;
		}


		// Store stem-to-stem distance-squared vector (for filtering candidates below)
		AtomVector stems = stem1 + stem2;
		vector<double> stemDistanceSq;
		for (uint c = 0; c < stem1.size();c++){
			/*
			  fprintf(stdout,"C: %1s %4d %3s\n",
			  stem1[c]->getChainId().c_str(),
			  stem1[c]->getResidueNumber(),
			  stem1[c]->getResidueName().c_str());
			*/
			for (uint n = 0; n < stem2.size();n++){

				double distSq = stem1(c).distance2(stem2(n));			
				/*
				  fprintf(stdout,"\tN: %1s %4d %3s = %8.3f\n",
				  stem2[n]->getChainId().c_str(),
				  stem2[n]->getResidueNumber(),
				  stem2[n]->getResidueName().c_str(),distSq);
				*/

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

		// Now loop over all fragments in database checking for ones of the correct size
		double tol = 4; // Tolerance of distance to be deviant from stems in Angstroms
		for (uint i = 0 ; i < fragDB.size()-(_numResiduesInFragment+stem1.size()+stem2.size());i++){

			// Get proposed ctermStem
			AtomVector ctermStem;
			for (uint n = 0; n < stem1.size();n++){
				ctermStem.push_back(fragDB[i+n]);
			}

			// Get proposed ntermStem
			AtomVector ntermStem;
			for (uint n = 0; n < stem2.size();n++){
				ntermStem.push_back(fragDB[i+stem1.size()+_numResiduesInFragment+n]);
			}

			// Check for same PDB
			if (ctermStem(0).getSegID() != ntermStem(ntermStem.size()-1).getSegID()){
				fprintf(stdout,"New PDB: %4s to %4s\n",ctermStem[0]->getSegID().c_str(),ntermStem(ntermStem.size()-1).getSegID().c_str());

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


			// align and print winning fragment

			// Make a copy of frag stem, so not to effect original atoms
			AtomVector fragStem;
			for (uint ct = 0; ct < ctermStem.size();ct++){
				fragStem.push_back(new Atom(ctermStem(ct)));
			}
			for (uint nt = 0; nt < ntermStem.size();nt++){
				fragStem.push_back(new Atom(ntermStem(nt)));
			}

			// Align the fragStem to stem to see if it matches well enough...
			Transforms tm;
			fragStem.saveCoor("pre");
			tm.align(fragStem,stems);
			
			double rmsd = fragStem.rmsd(stems);

			// Continue if the RMSD filter not passed
			if (rmsd > 0.1){
				continue;
			}
			fragStem.applySavedCoor("pre");

			//lastResults->writePdb("/tmp/preAdd.pdb");

			// Create new atoms for middle atoms..
			AtomVector tmp;
			index = 1;
			for (uint f = 0; f < ctermStem.size();f++){
				fragStem(f).setResidueNumber(index++);
				fragStem(f).setResidueName("FRG");
				fragStem(f).setChainId("A");
				tmp.push_back(&fragStem(f));

			}
		
			for (int f = 0; f < _numResiduesInFragment;f++){
				Atom *a =new Atom("CA","C");
				a->setCoor(fragDB(i+ctermStem.size()+f).getCoor());
				a->setResidueNumber(index++);
				a->setResidueName("FRG");
				a->setChainId("A");
				tmp.push_back(a);
			} 		

			for (uint f = 0; f < ntermStem.size();f++){
				fragStem(f+ctermStem.size()).setResidueNumber(index++);
				fragStem(f+ctermStem.size()).setResidueName("FRG");
				fragStem(f+ctermStem.size()).setChainId("A");
				tmp.push_back(&fragStem(f+ctermStem.size()));
			}

			Chain tmpChain;
			tmpChain.addAtoms(tmp);			

			
			tm.align(fragStem,stems,tmpChain.getAtoms());


			fprintf(stdout,"(%4s and chain %1s and resi %3d-%3d)  %8.3f\n",
				ctermStem[0]->getSegID().c_str(),
				ctermStem[0]->getChainId().c_str(),
				ctermStem[0]->getResidueNumber(),
				ntermStem[0]->getResidueNumber(),
			        rmsd);

			
			bool successful = true;
			if (fragType == caOnly){
				illegalQuads += bbqT.fillInMissingBBAtoms(tmpChain);

				if (illegalQuads > 0){
					successful = false;
					exit(0);
				}
				
			}

			if (successful){
				numFrags++;

				lastResults->addAtoms(tmpChain.getAtoms());
			}

			/*					
			stringstream ss;
			PDBWriter pout;
			pout.open(ss);
			pout.write(sys.getAtoms());
			pout.close();

			matchingFragments.push(pair<double,string>(rmsd,ss.str()));
			*/


		}

		//fprintf(stdout, "Best RMSD: %8.3f\n",matchingFragments.top().first);


		fprintf(stdout,"Number of succesful fragments: %10d, illegal quads: %10d\n",numFrags,illegalQuads);

	}
	
	
	return numFrags;
}
