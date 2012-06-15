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

#include "AtomContainer.h"
#include "AtomSelection.h"
#include "MslTools.h"


using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to select atoms in MSL
 *******************************************************************/

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_selecting_atoms <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to select atoms (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	string file = "example0003.pdb";
	file = (string)argv[1] + "/" + file;
	cout << "Create an AtomContainer and read the atoms from " << file << endl;

	// Read the PDB into a new AtomContainer called "atoms"
	AtomContainer atoms;
	if (!atoms.readPdb(file)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << file << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
		cout << "OK" << endl;
	}
	cout << endl;

	/************************************************************************
	 *   The AtomSelection is the object that operate selection over atoms using
	 *   logical statements (i.e. "NAME CA OR RESI 37")
	 *
	 *   The AtomSelection takes the atom pointers as a constructor argument
	 ************************************************************************/
	AtomSelection sel(atoms.getAtomPointers());

	/************************************************************************
	 *  The syntax is similar to PyMOL selection syntax.  Currently the following
	 *  keywords are supported
	 *    NAME      atom name (CA)
	 *    RESI      residue number (37, or even 37A with an insertion code)
	 *    RESN      residue name (LEU)
	 *    CHAIN     chain id (A)
	 *    HASCOOR   atom with defined coordinates
	 *    ALL       all atoms
	 *
	 *  The keywords can be mixed with the AND, OR and XOR operators for
	 *  complex logic (using parenthesis for priority)
	 *
	 *     NAME CA AND RESI 37
	 *     NAME CA AND (CHAIN B OR RESN LEU)
	 *  
	 *  Currently it is also supported a WITHIN X OF statement
	 *     NAME CA WITHIN 5 OF CHAIN B
	 *  will select all the CA that are within 5A from any atom on chain B,
	 *  (including all CAs on chain B)
	 ************************************************************************/

	cout << endl;
	cout << "=====================================" << endl;

	// Select all CA atoms and store as "allCAs" for easy retrieval
	cout << "Select by atom name:"<<endl;
	cout << "   Selection name: \"bb\"" << endl;
	cout << "   Logic         : \"name CA" << endl;
	sel.select("allCAs, name CA");

	// check the number of atoms selected
	cout << "   " << sel.size("allCAs") << " have been selected out of " << sel.size("all") << endl;

	// get the pointers of the atoms selectes and print them
	AtomPointerVector allCAs = sel.getSelection("allCAs");
	cout << allCAs;

	cout << endl;
	cout << "=====================================" << endl;

	// Expanding the selection logic
	cout << "Select all backbone atoms using selection logic:" << endl;
	cout << "   Selection name: \"bb\"" << endl;
	cout << "   Logic         : \"name CA OR name N OR name C OR name O\""<<endl;

	// the atom pointers can be obtained directly when the selection is created
	AtomPointerVector bb = sel.select("bb, name CA OR name N OR name C OR name O");
	cout << "   " << sel.size("bb") << " have been selected out of " << sel.size("all") << endl;
	cout << bb;

	cout << endl;
	cout << "=====================================" << endl;

	// using + to separate options
	cout << "The same selection can be written using + to separate multiple options:"<< endl;
	cout << "   Selection name: \"bb2\"" << endl;
	cout << "   Logic         : \"name CA+N+C+O\"" << endl; 
	AtomPointerVector bb2 = sel.select("bb2, name CA+N+C+O");
	cout << "   " << sel.size("bb2") << " have been selected out of " << sel.size("all") << endl;

	cout << endl;
	cout << "=====================================" << endl;

	// Expanding the selection logic with more keywords
	cout << "Complex selection: the side-chain carboxilic oxygen atoms on chain B:" << endl;
	cout << "   Selection name: carbox" << endl;
	cout << "   Logic         : \"((resn ASP AND name OD1+OD2) OR (resn GLU AND name OE1+OE2)) AND chain B\""<<endl;
	AtomPointerVector carbox = sel.select("carbox, ((resn ASP AND name OD1+OD2) OR (resn GLU AND name OE1+OE2)) AND chain B");
	cout << "   " << sel.size("carbox") << " have been selected out of " << sel.size("all") << endl;
	cout << carbox;

	cout << endl;
	cout << "=====================================" << endl;

	// Using previous selections
	cout << "Complex selection can use previously defined selections.  Let's select the backbone of LEU B 9 from the \"bb\" selection."<<endl;
	cout << "   Selection name: B9bb" << endl;
	cout << "   Logic         : \"bb AND resi 9 AND chain B\""<<endl;
	AtomPointerVector B9bb = sel.select("B9bb, bb AND resi 9 AND chain B");
	cout << "   " << sel.size("B9bb") << " have been selected out of " << sel.size("all") << endl;
	cout << B9bb;

	cout << endl;
	cout << "=====================================" << endl;

	// Using NOT operator
	cout << "The selection logic can use the NOT operator: let's select all CAs from residues excluding the hydrophobic residues"<<endl;
	cout << "   Selection name: CAnoHB" << endl;
	cout << "   Logic         : \"name CA and NOT resn ALA+VAL+LEU+ILE+MET+PHE+TRP+TYR\""<<endl;
	AtomPointerVector CAnoHB = sel.select("CAnoHB, name CA and NOT resn ALA+VAL+LEU+ILE+MET+PHE+TRP+TYR");
	cout << "   " << sel.size("CAnoHB") << " have been selected out of " << sel.size("all") << endl;
	cout << CAnoHB;

	cout << endl;
	cout << "=====================================" << endl;

	// Selecting ranges of residues, i.e. 2-4
	cout << "Residues can be selected by number using ranges. Let's select all CAs of residue 2-8 on chain B"<<endl;
	cout << "   Selection name: B2_8CA" << endl;
	cout << "   Logic         : \"name CA and resi 2-8 AND chain B\""<<endl;
	AtomPointerVector B2_8CA = sel.select("B2_8CA, name CA and resi 2-8 AND chain B");
	cout << "   " << sel.size("B2_8CA") << " have been selected out of " << sel.size("all") << endl;
	cout << B2_8CA;


	cout << endl;
	cout << "=====================================" << endl;

	// Selection using WITHIN X OF
	cout << "Atoms near other atoms can be selected using WITHIN X OF.  Let's select all CA within 7A from the CB or ALA B 10\n";
	cout << "   Selection name: wi7fromB10CB" << endl;
	cout << "   Logic         : \"name CA WITHIN 7 OF (chain B AND resi 10 AND name CB)\""<<endl;
	AtomPointerVector wi7fromB10CB = sel.select("wi7fromB10CB, name CA WITHIN 7 OF (chain B AND resi 10 AND name CB)");
	cout << "   " << sel.size("wi7fromB10CB") << " have been selected out of " << sel.size("all") << endl;
	cout << wi7fromB10CB;


	return 0;
}
