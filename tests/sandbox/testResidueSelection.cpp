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
#include <string>

#include "ResidueSelection.h"
#include "testData.h"
using namespace MSL;
using namespace std;


int main(){

         // Write out test PDB files
         writePdbFile();
	 System sys;
	 sys.readPdb("/tmp/pdbDimer.pdb");
        

	// Atom Selection Object
	ResidueSelection sel(sys);
	sel.setDebugFlag(true);

	// By default residues of 'all' built in, so we should be able to select it.
	vector<Residue *> allRes = sel.select("all");

	// Print out atoms that were selected
	cout << endl<<"Create selection named 'all' using the command: " << endl;
	cout << "   sel.select(\"all\")" << endl;
	for (uint i = 0; i < allRes.size();i++){
	  cout << "\tResidue " << allRes[i]->toString() << " ; flag: "<<allRes[i]->getSelectionFlag("all")<<endl;
	}

	

	vector<Residue *> hisRes = sel.select("his, resn HIS");
	cout << endl<<"Create selection using the command: " << endl;
	cout << "   sel.select(\"resn HIS\")" << endl;
	for (uint i = 0; i < allRes.size();i++){
	  cout << "\tResidue " << allRes[i]->toString() << " ; flag: "<<allRes[i]->getSelectionFlag("his")<<endl;
	}



	vector<Residue *> cbRes = sel.select("cb, name cb");
	cout << endl<<"Create selection using the command: " << endl;
	cout << "   sel.select(\"name CB\")" << endl;
	for (uint i = 0; i < allRes.size();i++){
	  cout << "\tResidue " << allRes[i]->toString() << " ; flag: "<<allRes[i]->getSelectionFlag("cb")<<endl;
	}


	vector<Residue *> cbcgRes = sel.select("cbcg, name CB and name CG");
	cout << endl<<"Create selection using the command: " << endl;
	cout << "   sel.select(\"name CB and name CG\")" << endl;
	for (uint i = 0; i < allRes.size();i++){
	  cout << "\tResidue " << allRes[i]->toString() << " ; flag: "<<allRes[i]->getSelectionFlag("cbcg")<<endl;
	}




	

	
	
}


