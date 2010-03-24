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
#include <iostream>
#include <cstdlib>

#include "System.h"
#include "AtomSelection.h"
#include "ResidueSelection.h"
#include "MslTools.h"


using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to select atoms and residues in MSL
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_selecting_atoms_and_residues <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to select atoms in the System (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	string file = "example0003.pdb";
	file = (string)argv[1] + "/" + file;
	cout << "Create an System and read the atoms from " << file << endl;

	// Read the PDB into a new System called "sys"
	System sys;
	if (!sys.readPdb(file)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	/************************************************************************
	 *   AtomSelection 
	 *
         *    Selection By Name:
	 *        Can select based on selection flags that are inside each Atom. 
	 *      By default every atom has an "all" selection flag set.
	 *      
	 *        You can create selection flags for atoms by:
	 *         1. Using Atom object:             a.setSelectionFlag("foobar")
	 *         2. Using AtomSelection object:    as.select("foobar, chain A")
	 *
	 *               
	 *    Selection By Property:          
	 *          Chain A:     as.select("chain A")
	 ************************************************************************/
	
	AtomSelection sel(sys.getAtomPointers());


	// By default every atom has been labeled as 'all'
	AtomPointerVector allAts = sel.select("all");
	cout << endl<<"Get selection named 'all' using the command: " << endl;
	cout << "   sel.select(\"all\")" << endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	
	for (AtomPointerVector::iterator it = sys.getAtomPointers().begin();it != sys.getAtomPointers().end(); it++){
		cout << "\tAtom " << *(*it) << " ; flag for all: "<<(*it)->getSelectionFlag("all")<<endl;
	}

	// Select chain A atoms and store as "chA" for easy retrieval
	cout << endl<<"Create a selection for chain A atoms 'chA'"<<endl;
	cout << "   sel.select(\"chA,chain A\")" << endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run

	AtomPointerVector chainAAts = sel.select("chA, chain A");
	for (AtomPointerVector::iterator it = sys.getAtomPointers().begin();it != sys.getAtomPointers().end(); it++){
		cout << "\tAtom " << *(*it) << " ; flag for chA: "<<(*it)->getSelectionFlag("chA")<<endl;
	}

	// Retreive selection "chA"
	if (sel.selectionExists("chA")){
	  AtomPointerVector get_chA = sel.getSelection("chA");
	  cout << endl<<"Retreive again the selected atoms by calling the selection with: " << endl;
	  cout << "   AtomPointerVector get_chA = sel.select(\"chA\");" << endl;
	  cout <<"Retreived Atom Count: "<<get_chA.size()<<" original selection Atom Count: "<<chainAAts.size()<<endl;
	}
	
	
	// Expanding the selection logic
	cout << "Select backbone atoms using selection logic"<<endl;
	cout << "    sel.select(\"bb, (((name CA OR name N) OR name C) OR name O)\");"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run

	AtomPointerVector backbone    = sel.select("bb, (((name CA OR name N) OR name C) OR name O)");
	for (uint i = 0; i < backbone.size();i++){
	  cout << *backbone[i]<<endl;
	}
	cout << endl;
	cout << "Select backbone atoms using short-hand selection logic"<<endl;
	cout << "    sel.select(\"bb2, name CA+N+C+O\");"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	AtomPointerVector bbShortHand = sel.select("bb2, name CA+N+C+O");
	for (uint i = 0; i < bbShortHand.size();i++){
	  cout << *bbShortHand[i]<<endl;
	}

	// Selections using names from previous selections
	cout << "Create a selection using the name of a previous selection.."<<endl;
	cout << "    AtomPointerVector test = sel.select(\"new, bb AND name O\");"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	AtomPointerVector test = sel.select("new, bb AND name O");
	for (uint i = 0; i < test.size();i++){
	  cout << *test[i]<<endl;
	}

	// Selections using NOT operator
	cout << "Create a selection using the NOT operator.."<<endl;
	cout <<"      AtomPointerVector notBB = sel.select(\"notBB, not name CA+C+O+N\")"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	AtomPointerVector notBB = sel.select("notBB, not name CA+C+O+N");
	for (uint i = 0; i < notBB.size();i++){
	  cout << *notBB[i]<<endl;
	}

	// Selections using range "-" operator.
	cout << "Create a selection using the '-' operator.."<<endl;
	cout <<"      AtomPointerVector res26 = sel.select(\"res26, resi 2-6\")"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	AtomPointerVector res26 = sel.select("res26, resi 2-6");
	for (uint i = 0; i < res26.size();i++){
	  cout << *res26[i]<<endl;
	}

	// Selection using distances to a sub-selection
	cout << "Create an atom selection named 'oxygensHbondSer', based on a distance from another selection 'name OG'\n";
	cout << "    AtomPointerVector oxygensHbondSer = sel.select(\"oxygensHbondSer, name O WITHIN 3.5 OF name OG\");"<<endl;
	cout << " Hit any character and enter to execute test"<<endl;
	cin >> file; // dummy line to wait for input to run
	AtomPointerVector oxygensHbondSer = sel.select("oxygensHbondSer, name O WITHIN 3.5 OF name OG");
	for (uint i = 0; i < oxygensHbondSer.size();i++){
	  cout << *oxygensHbondSer[i]<<endl;
	}


	return 0;
}
