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
#include "CartesianPoint.h"
#include "Matrix.h"
#include "CartesianGeometry.h"

using namespace MSL;
using namespace std;


Symmetry::Symmetry(){
}

Symmetry::~Symmetry(){
	deletePointers();
}

Symmetry::Symmetry(const Symmetry & _symmetry) {
	copy(_symmetry);
}

void Symmetry::operator=(const Symmetry & _symmetry) {
	copy(_symmetry);
}

void Symmetry::deletePointers() {
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
}

void Symmetry::copy(const Symmetry & _symmetry) {
	deletePointers();
	for (AtomPointerVector::const_iterator k=_symmetry.atoms.begin(); k!=_symmetry.atoms.end(); k++) {
		atoms.push_back(new Atom(**k));
	}
}

// Generic C_N symmetry written from C2 template by David Slochower
void Symmetry::applyCN(AtomPointerVector &_ats, int _N, const CartesianPoint & _primaryAxis){
	
	// clear any previous result
	deletePointers();

	// Find out how many matrices to make:
	//double angle = 360.0/_N;
	double angle = 0.0;
	int cycle = 1;
	atoms.clear();
	for (uint i = 0; i < _ats.size(); i++) {
		atoms.push_back(new Atom(_ats(i)));
	}
	//for (double j = angle; j < 360; j += angle) {
	for (unsigned int i = 1; i < _N; i++) {
		angle += 360.0/_N;

		Matrix rotMat = CartesianGeometry::getRotationMatrix(angle, _primaryAxis);
		AtomPointerVector axisRot;
		for (uint i = 0; i < _ats.size(); i++) {
			axisRot.push_back(new Atom(_ats(i)));
		}
		//zRot.rotate(zMat);
		Transforms tr;
		tr.rotate(axisRot, rotMat);
		char chainLetter = int('A') + (cycle);

		for (uint i = 0; i < _ats.size(); i++){
			string chainID;
			chainID += chainLetter;
			axisRot[i]->setChainId(chainID);
			atoms.push_back(axisRot[i]);
		}
		cycle++;
		axisRot.clear();
	}
}

void Symmetry::applyDN(AtomPointerVector &_ats, int _N){
	// Programs by default may assume that the primary axis is the z-axis, and the template atoms
        // are on the x-axis.  A proper secondary axis would then be at the angle of 180 degrees/n.

	CartesianPoint defaultPrimaryAxis(0,0,1);

        double secondaryAngle = M_PI/_N;
	CartesianPoint defaultSecondaryAxis(cos(secondaryAngle),sin(secondaryAngle),0);

	applyDN(_ats, _N, defaultPrimaryAxis, defaultSecondaryAxis);
}

void Symmetry::applyDN(AtomPointerVector &_ats, int _N, const CartesianPoint & _primaryAxis, const CartesianPoint & _secondaryAxis){

	/*
	  Formal description of DN please..
          Use this for antiparallel bundles/coils, where 2*N is the total number of helices/objects.
          First a parallel half of C_N symmetry around the primary axis, then it is rotated 180 degrees
          around the secondaryAxis.   
	 */

	// clear any previous result
	deletePointers();

	double rotationAngle = 180.;
	Matrix rotMat = CartesianGeometry::getRotationMatrix(rotationAngle, _secondaryAxis);

	applyCN(_ats,_N, _primaryAxis);

	string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

	AtomPointerVector axisRot;
	for (uint i = 0; i < atoms.size(); i++) {
		axisRot.push_back(new Atom(atoms(i)));
	}
	Transforms tr;
	tr.rotate(axisRot, rotMat);

	for (uint i = 0; i < axisRot.size(); i++){
		int index = alphabet.find(axisRot[i]->getChainId());
		axisRot[i]->setChainId(alphabet.substr(index+_N,1));
		atoms.push_back(axisRot[i]);
	}
	axisRot.clear();
}

