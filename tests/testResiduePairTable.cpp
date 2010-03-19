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


#include "ResiduePairTable.h"
#include "ResiduePairTableReader.h"
#include "PDBReader.h"
#include "System.h"
#include "testData.h"


using namespace std;

using namespace MSL;



int main(){

	cout << "Read /export/home/dwkulp/projects/PeptideReceptors/skolnick.mulipot.txt"<<endl;
	ResiduePairTableReader rptr("/export/home/dwkulp/projects/PeptideReceptors/skolnick.mulipot.txt");
	
	rptr.open();
	rptr.read();
	rptr.close();

	cout << "Read a string pdb 'fourHelixBundle'"<<endl;
	PDBReader pdbin(fourHelixBundle);
	pdbin.read();
	pdbin.close();

	cout << "Get Statistics for each residue"<<endl;
	ResiduePairTable &rpt = rptr.getResiduePairTable();

	System sys(pdbin.getAtomPointers());
	cout << "Number of residues: "<<sys.getChain("A").size()<<endl;
	for (uint i = 0; i < sys.getChain("A").size();i++){
		for (uint j = 0; j < sys.getChain("A").size();j++){
			if (i == j) continue;


			Residue & r1 = sys.getChain("A").getResidue(i);
			Residue & r2 = sys.getChain("A").getResidue(j);

			double val = rpt.getValue(r1.getResidueName(), r2.getResidueName());
			fprintf(stdout, "%-15s %1s %3s %3d - %1s %3s %3d = %8.3f\n", "4HelixBundle" , "A", r1.getResidueName().c_str(),r1.getResidueNumber(), "A", r2.getResidueName().c_str(),r2.getResidueNumber(),val);
		}

		
		
	}



}
