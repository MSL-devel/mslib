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


#include "HelixFusion.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "System.h"
#include "testData.h"

#include <ostream>


using namespace std;

using namespace MSL;


int main(){

	cout << "Read a string pdb 'fourHelixBundle'"<<endl;
	PDBReader pdbin;
	pdbin.read(fourHelixBundle);
	pdbin.close();

	System sys;
	sys.addAtoms(pdbin.getAtomPointers());
	sys.writePdb("/tmp/preFusion.pdb");
	HelixFusion hf;
	cout << "Chain Sizes: "<<sys("A").positionSize()<<" "<<sys("B").positionSize()<<endl;
	hf.setChains(sys("A"),sys("B"));

	hf.fusionByAtomicAlignment(0.3);




	for (uint f = 0; f< hf.getNumberFusions();f++){
		Chain &newChain = hf.getFusedChain(f);
		PDBWriter pdbout;
		char name[80];
		sprintf(name,"/tmp/fusedChain-%04d.pdb",f);
		pdbout.open((string)name);
		pdbout.write(newChain.getAtomPointers());
		pdbout.close();
	}



	hf.fusionByHelicalFrames();
}
