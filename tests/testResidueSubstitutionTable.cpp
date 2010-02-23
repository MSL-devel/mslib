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

#include "ResidueSubstitutionTable.h"
#include "ResidueSubstitutionTableReader.h"
#include "PDBReader.h"
#include "System.h"
#include "testData.h"


using namespace std;

using namespace MSL;



int main(int argc, char **argv){
        if (argc < 3) {
            cout << "Usage: testResidueSubstitutionTable <AASubstitutionTable> <pdb>" << endl;
            exit(0);
        }

	cout << "Read " << string(argv[1]) << endl;
	ResidueSubstitutionTableReader rstr(argv[1]);
	
	rstr.open();
	rstr.read();
	rstr.close();

	cout << "Read pdb " << string(argv[2]) << endl;
	PDBReader pdbin(argv[2]);
        pdbin.open();
	pdbin.read();
	pdbin.close();

	cout << "See if the pdb has any non-natural amino acids."<<endl;
	ResidueSubstitutionTable &rst = rstr.getResidueSubstitutionTable();

	System sys(pdbin.getAtoms());
	cout << "Number of residues: "<<sys.getChain("A").size()<<endl;

        PDBWriter pdbout("test.pdb");
        pdbout.open();

	for (uint i = 0; i < sys.residueSize();i++){
            Residue &currResidue = sys.getResidue(i);
            Residue newResidue = rst.replaceResidue(currResidue);
            pdbout.write(newResidue.getAtoms(), false);

            if(rst.isResidueInSubstitutionTable(currResidue))
                cout << "Residue " << currResidue.getResidueName() << " found in list." << endl;
            else
                cout << "Residue " << currResidue.getResidueName() << " NOT found in list." << endl;
	}

        pdbout.close();


}
