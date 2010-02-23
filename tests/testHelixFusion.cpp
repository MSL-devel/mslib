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
	sys.addAtoms(pdbin.getAtoms());
	sys.writePdb("/tmp/preFusion.pdb");
	HelixFusion hf;
	cout << "Chain Sizes: "<<sys("A").size()<<" "<<sys("B").size()<<endl;
	hf.setChains(sys("A"),sys("B"));

	hf.fusionByAtomicAlignment(0.3);




	for (uint f = 0; f< hf.getNumberFusions();f++){
		Chain &newChain = hf.getFusedChain(f);
		PDBWriter pdbout;
		char name[80];
		sprintf(name,"/tmp/fusedChain-%04d.pdb",f);
		pdbout.open((string)name);
		pdbout.write(newChain.getAtoms());
		pdbout.close();
	}



	hf.fusionByHelicalFrames();
}
