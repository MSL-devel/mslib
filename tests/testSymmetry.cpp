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

#include "Symmetry.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "CartesianGeometry.h"
#include "Matrix.h"
#include "testData.h"

using namespace std;

int main(){

	// Coiling an ideal, Z-aligned helix
	PDBReader pin;
	pin.read(idealHelix);

	AtomPointerVector ideal;
	ideal = pin.getAtoms();
	pin.close();

	PDBReader pin2;
	pin2.read(idealHelix);
	AtomPointerVector ideal2;
	ideal2 = pin2.getAtoms();
	pin2.close();

		  
	CartesianPoint x(6.5,0,0);
	ideal.translate(x);


	Matrix zRot = CartesianGeometry::instance()->getZRotationMatrix(45);
	CartesianPoint x2(10,0,0);
	ideal2.translate(x2);
	ideal2.rotate(zRot);


	PDBWriter pout;
	Symmetry s;
	for (uint i = 1; i < 12; i++){


		s.applyCN(ideal,i);

		stringstream fname;
		fname << "C"<<i<<".pdb";
		pout.open(fname.str());
		pout.write(s.getAtoms());
		pout.close();


		if (i % 2 == 1){
			s.applyDN(ideal2,i);

			fname.str("");
			fname << "D"<<i<<".pdb";
			pout.open(fname.str());
			pout.write(s.getAtoms());
			pout.close();
		}
	}




}
