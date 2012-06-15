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
#include "RegEx.h"

using namespace std;
using namespace MSL;
int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_regular_expressions <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to use regular expressions in MSL (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string file = "example0004.pdb";
	file = (string)argv[1] + "/" + file;
	cout << "Create an AtomContainer and read the atoms from " << file << endl;

	System sys;
	if (!sys.readPdb(file)) {
	  // reading failed, error handling code here
	  cerr << "ERROR could not read in "<<file<<endl;
	  exit(0);
	}
 

	// Check to make sure chain A exits in sys
	if (!sys.chainExists("A")){
	  // error code here.
	  cerr << "ERROR chain A does not exist in file "<<file<<endl;
	  exit(0);
	}
 
	// Get a Chain object
	Chain &ch = sys.getChain("A");
 
	// Regular Expression Object
	RegEx re;
 
	// Find 3 Prolines surrounded by two Glycines on one side and three Glycines on the other
	string regex = "V{2}IL";
 
	// Now do a sequence search...
	vector<pair<int,int> > matchingResidueIndices = re.getResidueRanges(ch,regex);
 
 
	// Loop over each match.
	for (uint m = 0; m < matchingResidueIndices.size();m++){
 
	  // Loop over each residue for this match
	  int match = 1;
	  for (uint r = matchingResidueIndices[m].first; r <= matchingResidueIndices[m].second;r++){
 
	    // Get the residue
	    Residue &res = ch.getResidue(r);
 
	    // .. do something cool with matched residues ...
	    cout << "MATCH("<<match<<"):  RESIDUE: "<<res.toString()<<endl;
	  }
	}

}


