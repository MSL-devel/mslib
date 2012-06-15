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

#include "Symmetry.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "CartesianGeometry.h"
#include "Matrix.h"
#include "Transforms.h"
#include "testData.h"

using namespace std;

using namespace MSL;


int main(){

	// Coiling an ideal, Z-aligned helix
	PDBReader pin;
	pin.read(idealHelix);

	AtomPointerVector ideal;
	ideal = pin.getAtomPointers();
	pin.close();

	PDBReader pin2;
	pin2.read(idealHelix);
	AtomPointerVector ideal2;
	ideal2 = pin2.getAtomPointers();
	pin2.close();

	Transforms tr;
		  
	CartesianPoint x(6.5,0,0);
//	ideal.translate(x);
	tr.translate(ideal, x);


	Matrix zRot = CartesianGeometry::getZRotationMatrix(45);
	CartesianPoint x2(10,0,0);
//	ideal2.translate(x2);
	tr.translate(ideal2, x2);
//	ideal2.rotate(zRot);
	tr.rotate(ideal2, zRot);


	PDBWriter pout;
	Symmetry s;
	for (uint i = 1; i < 12; i++){


		s.applyCN(ideal,i);

		stringstream fname;
		fname << "C"<<i<<".pdb";
		pout.open(fname.str());
		pout.write(s.getAtomPointers());
		pout.close();


		if (i % 2 == 1){
			s.applyDN(ideal2,i);

			fname.str("");
			fname << "D"<<i<<".pdb";
			pout.open(fname.str());
			pout.write(s.getAtomPointers());
			pout.close();
		}
	}




}
