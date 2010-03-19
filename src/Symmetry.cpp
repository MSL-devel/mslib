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

#include "Symmetry.h"
#include "CartesianGeometry.h"

using namespace MSL;
using namespace std;


Symmetry::Symmetry(){
}

Symmetry::~Symmetry(){
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
}




// Generic C_N symmetry written from C2 template by David Slochower
void Symmetry::applyCN(AtomPointerVector &_ats, int _N){
	
	// Find out how many matrices to make:
	double angle = 360.0/_N;
	int k = 1;
	atoms.clear();
        for (uint i = 0; i < _ats.size(); i++) {
	    atoms.push_back(new Atom(_ats(i)));
	}
	for (double j = angle; j < 360; j += angle) {

	  Matrix zMat = CartesianGeometry::instance()->getZRotationMatrix(j);
	  AtomPointerVector zRot;
	  for (uint i = 0; i < _ats.size(); i++) {
	    zRot.push_back(new Atom(_ats(i)));
	  }
	  //zRot.rotate(zMat);
	  Transforms tr;
	  tr.rotate(zRot, zMat);
	  char chainLetter = int('A') + (k);

	  for (uint i = 0; i < _ats.size(); i++){
	    string chainID;
	    chainID += chainLetter;
	    zRot(i).setChainId(chainID);
	    atoms.push_back(zRot[i]);
	  }
	  k++;
	  zRot.clear();
	}
}

void Symmetry::applyD2(AtomPointerVector &_ats){

	/*
	  Formal description of D2 please..
	 */

	// Rotate 180 around each axis
	Matrix xMat = CartesianGeometry::instance()->getXRotationMatrix(180);
	Matrix yMat = CartesianGeometry::instance()->getYRotationMatrix(180);
	Matrix zMat = CartesianGeometry::instance()->getZRotationMatrix(180);

	atoms.clear();
	AtomPointerVector xRot,yRot,zRot;
	for (uint i =0; i < _ats.size();i++){
		atoms.push_back(new Atom(_ats(i)));
		xRot.push_back(new Atom(_ats(i)));
		yRot.push_back(new Atom(_ats(i)));
		zRot.push_back(new Atom(_ats(i)));
	}

	Transforms tr;
//	xRot.rotate(xMat);
//	yRot.rotate(yMat);
//	zRot.rotate(zMat);
	tr.rotate(xRot, xMat);
	tr.rotate(yRot, yMat);
	tr.rotate(zRot, zMat);


	for (uint i =0; i < xRot.size();i++){
		xRot(i).setChainId("B");
		atoms.push_back(xRot[i]);
	}
	for (uint i =0; i < yRot.size();i++){
		yRot(i).setChainId("C");
		atoms.push_back(yRot[i]);
	}
	for (uint i =0; i < zRot.size();i++){
		zRot(i).setChainId("D");
		atoms.push_back(zRot[i]);
	}
	
	xRot.clear();	
	yRot.clear();	
	zRot.clear();	
	
}

void Symmetry::applyDN(AtomPointerVector &_ats, int _N){

	/*
	  Formal description of D2 please..
	 */

	
	applyCN(_ats,_N);

	int atsize = atoms.size();
	string alphabet = "ABCDEFGHIJKLMNOPQRSTUVQXYZ";
	for (uint i = 0 ; i < atsize;i++){
		atoms.push_back(new Atom(atoms(i)));

		atoms.back()->setCoor(-atoms.back()->getX(),-atoms.back()->getY(),-atoms.back()->getZ());

		int index = alphabet.find(atoms.back()->getChainId());
		atoms.back()->setChainId(alphabet.substr(index+_N,1));
	}

	
}

// void Symmetry::applyC2(AtomPointerVector &_ats){
// 	/*
// 	  Formal description of C2 please..
// 	 */

// 	// Single Axis of rotation, 360/2 = 180 degrees.
// 	Matrix zMat = CartesianGeometry::instance()->getZRotationMatrix(180);

// 	atoms.clear();
// 	AtomPointerVector zRot;
// 	for (uint i =0; i < _ats.size();i++){	
// 		atoms.push_back(new Atom(_ats(i)));
// 		zRot.push_back(new Atom(_ats(i)));
// 	}

// 	zRot.rotate(zMat);
// 	for (uint i =0; i < zRot.size();i++){

// 		zRot(i).setChainId("B");
// 		atoms.push_back(zRot[i]);
// 	}

// 	zRot.clear();
// }

// void Symmetry::applyC3(AtomPointerVector &_ats){
// 	/*
// 	  Formal description of C3 please..
// 	 */

// 	// Single Axis of rotation, 360/3 = 120 degrees.
// 	Matrix zMat1 = CartesianGeometry::instance()->getZRotationMatrix(120);
// 	Matrix zMat2 = CartesianGeometry::instance()->getZRotationMatrix(240);

// 	atoms.clear();
// 	AtomPointerVector zRot1;
// 	AtomPointerVector zRot2;
// 	for (uint i =0; i < _ats.size();i++){	
// 		atoms.push_back(new Atom(_ats(i)));
// 		zRot1.push_back(new Atom(_ats(i)));
// 		zRot2.push_back(new Atom(_ats(i)));
// 	}

// 	zRot1.rotate(zMat1);
// 	for (uint i =0; i < _ats.size();i++){

// 		zRot1(i).setChainId("B");
// 		atoms.push_back(zRot1[i]);

// 	}

// 	zRot2.rotate(zMat2);
// 	for (uint i =0; i < _ats.size();i++){

// 		zRot2(i).setChainId("C");
// 		atoms.push_back(zRot2[i]);	
// 	}

// 	zRot1.clear();
// 	zRot2.clear();
// }


// void Symmetry::applyCNanti(AtomPointerVector &_ats, int _N){


// 	// Find out how many matrices to make:
// 	double angle = 360.0/_N;




// 	atoms.clear();
//         for (uint i = 0; i < _ats.size(); i++) {
// 	    atoms.push_back(new Atom(_ats(i)));
// 	}

// 	int k = 1;
// 	for (double j = angle; j < 360; j += angle) {


//   	  // Get Zrot Matrix
// 	  Matrix zMat = CartesianGeometry::instance()->getZRotationMatrix(j);

// 	  // Apply Zrot Matrix
// 	  AtomPointerVector zRot;
// 	  for (uint i = 0; i < _ats.size(); i++) {
// 	    zRot.push_back(new Atom(_ats(i)));
// 	  }
// 	  zRot.rotate(zMat);

// 	  // Anti-parallel when approriate... (every other one)
// 	  if (k % 2 == 1){
// 		  for (uint i = 0; i < _ats.size(); i++) {
// 			  zRot(i).setCoor(zRot(i).getX(), zRot(i).getY(), -zRot(i).getZ());
// 		  }
// 	  }

// 	  char chainLetter = int('A') + (k++);

// 	  for (uint i = 0; i < _ats.size(); i++){
// 	    string chainID;
// 	    chainID += chainLetter;
// 	    zRot(i).setChainId(chainID);
// 	    atoms.push_back(zRot[i]);
// 	  }

// 	  // Clear tmp zRot atom vector
// 	  zRot.clear();
// 	}

	

// }
