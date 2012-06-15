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

#include "MslTools.h"

using namespace std;

using namespace MSL;


#include "VectorPair.h"
#include "SasaCalculator.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"

int main(int argc, char *argv[]) {
	if (argc < 1){
		cerr << "USAGE: testVectorPair PDB\n";
		exit(0);
	}
	string file = (string)argv[1];
	
	System sys;
	sys.readPdb(file);

	for (uint i = 0; i < sys.chainSize();i++){
		Chain &ch = sys.getChain(i);

		for (uint r1 = 0; r1 < ch.positionSize();r1++){
			Residue &res1 = ch.getResidue(r1);

			if (res1.getResidueName() == "GLY") continue;

			if (!(res1.atomExists("CA") && res1.atomExists("CB"))) continue;
			for (uint r2 = r1+1;r2 < ch.positionSize();r2++){
				Residue &res2 = ch.getResidue(r2);

				if (res2.getResidueName() == "GLY") continue;
				if (!(res2.atomExists("CA") && res2.atomExists("CB"))) continue;

				VectorPair vp(res1("CA").getCoor(), res1("CB").getCoor(),res2("CA").getCoor(),res2("CB").getCoor(),res1.getPositionId(),res2.getPositionId());

				vp.calcAll();
				
				// Get SASA?
				AtomPointerVector atoms;
				atoms = res1.getAtomPointers() + res2.getAtomPointers();
				SasaCalculator sas(atoms);
				sas.calcSasa();
				double sasa = sas.getTotalSasa();
				
				atoms.clear();
				atoms = res1.getAtomPointers();
				SasaCalculator sas1(atoms);
				sas1.calcSasa();
				sasa -= sas1.getTotalSasa();
				
				atoms.clear();
				atoms = res2.getAtomPointers();
				SasaCalculator sas2(atoms);
				sas2.calcSasa();
				sasa -= sas2.getTotalSasa();
				
				if (abs(sasa) < 0.0001){
					sasa = 0.0;
				} 

				fprintf(stdout, "%s %s %1s %1s %-8.3f %-5d %3s %3s\n",file.c_str(),vp.toString().c_str(),res1.getChainId().c_str(),res2.getChainId().c_str(),sasa,res2.getResidueNumber() - res1.getResidueNumber(),res1.getResidueName().c_str(),res2.getResidueName().c_str());
			}
		}
	}
			       
	
	
}
