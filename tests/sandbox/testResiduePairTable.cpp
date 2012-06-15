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
	cout << "Number of residues: "<<sys.getChain("A").positionSize()<<endl;
	for (uint i = 0; i < sys.getChain("A").positionSize();i++){
		for (uint j = 0; j < sys.getChain("A").positionSize();j++){
			if (i == j) continue;


			Residue & r1 = sys.getChain("A").getResidue(i);
			Residue & r2 = sys.getChain("A").getResidue(j);

			double val = rpt.getValue(r1.getResidueName(), r2.getResidueName());
			fprintf(stdout, "%-15s %1s %3s %3d - %1s %3s %3d = %8.3f\n", "4HelixBundle" , "A", r1.getResidueName().c_str(),r1.getResidueNumber(), "A", r2.getResidueName().c_str(),r2.getResidueNumber(),val);
		}

		
		
	}



}
