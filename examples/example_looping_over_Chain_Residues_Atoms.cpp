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
#include <iostream>
#include <cstdlib>

#include "System.h"
#include "MslTools.h"


using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to loop over chains, residues and
 *  atoms in the System
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_looping_over_Chain_Residues_Atoms <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to loop over chains, residues, atoms in the System (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string file = "example0002.pdb";
	file = (string)argv[1] + "/" + file;
	cout << "Create an System and read the atoms from " << file << endl;

	// read the PDB into a new System "container"
	System container;
	if (!container.readPdb(file)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	// the sequence can be conveniently be printed using the << operator
	cout << "Print the sequence of the molecule" << endl;
	cout << container << endl;
	cout << endl;
	cout << "=============================" << endl;

	
	/************************************************************************
	 *  The () operator (with unsigned int) gives the element down the hierarchy
	 *  in the following order 
	 * 
	 *  System
	 *      |--> Chain: A, B, C...
	 *              |--> Position: 1, 2, 3...
	 *                         |--> Identity: ALA, GLY, LEU...
	 *                                    |--> Atom: CA, CB, CG...
	 *
	 *  The support for multiple identities (with only one being active at any time)
	 *  is mainly for protein design.  A PDB file normally has only one identity per
	 *  position.  The second example bypasses the identities completely.
	 ************************************************************************/
	Chain & chn = container(0); // get the first chain in the System (chain "A") by reference
	Position & pos = chn(0);    // get the first position in the Chain (1) by reference
	Residue & res = pos(0);     // get the first identity in the Position (ASP) by reference
	Atom & atm = res(0);        // get the first atom in the Residue (N) by reference

	cout << "First chain of System: " << chn << endl;
	cout << "First position of the Chain: " << pos << endl;
	cout << "First identity of the Position: " << res << endl;
	cout << "First atom of the Residue: " << atm << endl;
	cout << endl;
	cout << "=============================" << endl;

	/************************************************************************
	 *  This code goes down the hierarchy in the system and prints the elements
	 ************************************************************************/
	cout << "Print the hierarchy System/Chain/Position/Identity/Atom" << endl;
	cout << endl;
	cout << "The System has " << container.chainSize() << " chains" << endl;
	// loop over each chain and get the number of positions
	for (unsigned int i=0; i<container.chainSize(); i++) {
		cout << "   Chain " << i << " (id: " << container(i).getChainId() << ") has " << container(i).positionSize() << " positions" << endl;

		// loop over each position and get the number of identities
		for (unsigned int j=0; j<container(i).positionSize(); j++) {
			cout << "      Position " << j << " (id: " << container(i)(j).getPositionId() << ") has " << container(i)(j).identitySize() << " identity" << endl;

			// loop over each identity and print the name
			for (unsigned int k=0; k<container(i)(j).identitySize(); k++) {
				cout << "         Identity " << k << " (id: " << container(i)(j)(k).getIdentityId() << ") is a " << container(i)(j)(k).getResidueName() << endl;

				// loop over each atom and print the name
				for (unsigned int l=0; l<container(i)(j)(k).size(); l++) {
					cout << "            Atom " << l << " (id: " << container(i)(j)(k)(l).getAtomOfIdentityId() << ") is a " << container(i)(j)(k)(l).getName() << endl;
				}
			}
		}
	}

	cout << endl;
	cout << "=============================" << endl;

	/************************************************************************
	 *  Using the [] operator at all levels (System, Chain, Position, Residue) gives
	 *  an Atom
	 * 
	 ************************************************************************/
	Atom & aSys = container[0]; // get the first atom (N) of tye system by reference
	Atom & aChn = chn[0];       // get the first atom (N) of the chain by reference
	Atom & aPos = pos[0];       // get the first atom (N) of the position by reference
	Atom & aRes = res[0];       // get the first atom (N) of the identity by reference

	cout << "First atom of System: " << aSys << endl;
	cout << "First atom of the Chain: " << aChn << endl;
	cout << "First atom of the Position: " << aPos << endl;
	cout << "First atom of the Residue: " << aRes << endl;
	cout << endl;
	cout << "=============================" << endl;

	/************************************************************************
	 *  While looping throught the structure can skip the identity level by using 
	 *  the [] at the Position (note: instead of .size(), which gives the number of
	 * identities, we use atomSize() to get the number of atoms in the active identity.
	 ************************************************************************/
	cout << "Skip the identity level (a regular PDB has 1 identity in each position and print the hierarchy System/Chain/Position/Atom" << endl;
	cout << endl;
	cout << "The System has " << container.chainSize() << " chains" << endl;
	// loop over each chain and get the number of positions
	for (unsigned int i=0; i<container.chainSize(); i++) {
		cout << "   Chain " << i << " (id: " << container(i).getChainId() << ") has " << container(i).positionSize() << " positions" << endl;

		// loop over each position and get the number of atoms
		for (unsigned int j=0; j<container(i).positionSize(); j++) {
			cout << "      Position " << j << " (id: " << container(i)(j).getPositionId() << ") is a " << container(i)(j).getResidueName() << " and has " << container(i)(j).atomSize() << " atoms" << endl;

			// loop over each atom and print the name, note using atomSize instead of size
			for (unsigned int k=0; k<container(i)(j).atomSize(); k++) {
				// note using the [] operator
				cout << "            Atom " << k << " (id: " << container(i)(j)[k].getAtomId() << ") is a " << container(i)(j)[k].getName() << endl;
			}
		}
	}

	return 0;
}
