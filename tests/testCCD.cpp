/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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
#include <vector>
#include <map>

#include "PDBReader.h"
#include "PDBWriter.h"
#include "CartesianGeometry.h"
#include "AtomSelection.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "CCD.h"

#include "testData.h"

using namespace std;

using namespace MSL;




int main() {

	// Write a test pdb file '/tmp/pdbDimer.pdb'
	writePdbFile();

	PDBReader rin;
	rin.open("/tmp/pdbDimer.pdb");
	rin.read();
	rin.close();

	AtomSelection sel(rin.getAtomPointers());
	//AtomPointerVector &loop = sel.select("(resi 25-29 and name CA)");
	AtomPointerVector &loop = sel.select("(resi 22-32 and name CA)");

	
	PDBWriter wout;
	wout.open("/tmp/loop.pdb");
	wout.write(loop);
	wout.close();


	Atom *fixedEnd = new Atom(loop(loop.size()-1));

	for (uint i = 0; i < loop.size()-1;i++){
		fprintf(stdout,"%4d %4d %8.3f\n", i,i+1,loop(i).distance(loop(i+1)));
	}

	// Now move the loop
	Transforms t;
	RandomNumberGenerator rng;
	rng.setTimeBasedSeed();
	for (uint i = 1; i < loop.size()-2;i++){
		
		CartesianPoint axis = loop(i+1).getCoor()-loop(i).getCoor();
		axis = axis.getUnit();
		axis += loop(i).getCoor();


		double angle = rng.getRandomIntLimit(15);
		for (uint j=i+2; j < loop.size();j++){

			t.rotate(loop(j),angle, axis, loop(i).getCoor());
		}

		

	}


	wout.open("/tmp/loop.random.pdb");
	wout.write(loop);
	wout.close();



	CCD abinitio;
	abinitio.closeFragment(loop,*fixedEnd);
	
	for (uint i = 0; i < loop.size()-1;i++){
		fprintf(stdout,"%4d %4d %8.3f\n", i,i+1,loop(i).distance(loop(i+1)));
	}
	wout.open("/tmp/loop.closed.pdb");
	wout.write(loop);
	wout.close();
}
