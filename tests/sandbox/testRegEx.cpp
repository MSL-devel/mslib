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

#include "testData.h"
#include "RegEx.h"
#include "System.h"

using namespace MSL;
using namespace std;



int main(){



	writePdbFile();
	System sys;
	sys.readPdb("/tmp/xtalLattice.pdb");
	             
	cout << "\n****\nSearch chain A of /tmp/xtalLattice for A..K\n"<<endl;
	if (sys.chainExists("A")){
		RegEx re;
		vector<pair<int,int> > matches = re.getResidueRanges(sys("A"),"A..K");

		cout <<endl;
		for (uint i = 0; i < matches.size();i++){

			for (uint j = matches[i].first; j <= matches[i].second;j++){
				Residue &r = sys("A").getResidue(j);

				cout << "R: "<<r.toString()<<endl;
			}
			cout << " -- "<<endl;
		}

	}




	// Generic string matching
	string regExpression = "^CLUSTAL W (\\S+)\\s+([\\S\\s]+)$";
        string lineToMatch="CLUSTAL W 2.1 multiple sequence alignment";
        string lineNotToMatch="LUTAL W 2.1 multiple sequence alignment";

	vector<string> results;

	// Line which should match RegEx
	if (MslTools::regex(lineToMatch,regExpression,results)){

	  cout << "PASS test 1\n";

	  // Print out the matched '()' substrings
	  for (uint i = 0; i < results.size();i++){
	    cout << "Match["<<i<<"]: "<<results[i]<<endl;
	  }

	} else {
	  cerr << "ERROR did not match line properly.\n";
	}


	// Line which should NOT match RegEx
	if (MslTools::regex(lineNotToMatch,regExpression,results)){
	  cerr << "ERROR matched line when it was not suppose to\n";
	} else {
	  cout << "PASS test 2 , did not match a string that it should not match\n";
	}


}
