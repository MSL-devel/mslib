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
#include "CCD.h"
#include "AtomContainer.h"

using namespace MSL;
using namespace std;


CCD::CCD(){
	useBBQ = false;
}

CCD::CCD(string _BBQTableForBackboneAtoms){
	bbqT.openReader(_BBQTableForBackboneAtoms);
	useBBQ = true;
}

CCD::~CCD(){

}

void CCD::localSample(AtomPointerVector &_av,int numFragments, int maxDegree){

	// Take one atom off each end for BBQ purposes
	Atom *forBBQ_N = new Atom(_av(_av.size()-1));
	Atom *forBBQ_C = new Atom(_av(0));

	_av.erase(_av.begin());
	_av.erase(_av.end()-1);

	AtomPointerVector forBBQ;
	forBBQ.push_back(forBBQ_N);
	forBBQ.push_back(forBBQ_C);



	// Store output PDB string
	stringstream ss;

	RandomNumberGenerator rng;
	rng.setTimeBasedSeed();
	PDBWriter pout;

	AtomContainer allAtoms;

	_av.saveCoor("pre");
	Atom *fixedAtom = NULL;
	bool reversed = false;
	for (uint d = 0; d < numFragments; d++){
		reversed = false;

		// Figure out which way to order the loop
		if (rng.getRandomDouble() > 0.5){
			std::reverse(_av.begin(), _av.end());
			reversed = true;
		}

		// Fixed atom is last Atom.
		fixedAtom = new Atom(_av(_av.size()-1));

		// Now move the _av
		Transforms t;
		for (uint i = 1; i < _av.size()-2;i++){
		
			CartesianPoint axis = _av(i+1).getCoor();

		        double angle = rng.getRandomInt(maxDegree);
			for (uint j=i+2; j < _av.size();j++){
				t.rotate(_av(j),angle, axis, _av(i).getCoor());
			}

		
		}

		closeFragment(_av,*fixedAtom);

		if (reversed){
			std::reverse(_av.begin(), _av.end());
		}

		Chain c;

		if (useBBQ){
			cout << "BBQ TIME"<<endl;
			c.addAtoms(forBBQ);
			c.addAtoms(_av);
			bbqT.fillInMissingBBAtoms(c);
		} else {
			c.addAtoms(_av);
		}
		

		allAtoms.addAtoms(forBBQ);
		allAtoms.addAtoms(c.getAtomPointers());


		ss << "MODEL "<<endl;
		stringstream tmp;
		pout.open(tmp);
		pout.write(c.getAtomPointers());
		pout.close();

		ss << tmp.str()<< "ENDMDL\n";


		_av.applySavedCoor("pre");
		delete(fixedAtom);
		fixedAtom = NULL;
	}

	closedSystem.addAtoms(allAtoms.getAtomPointers());
	closedSystem_NMRString = ss.str();
}

void CCD::closeFragment(AtomPointerVector &_av, Atom &_fixedEnd){
	

	int numIterations = 0;
	bool converged = false;
	while (true) {
		converged = false;

		for (uint i = _av.size()-3; i > 1;i--){
			//cout << "Pivot index: "<<i<<" "<<_av(i).getResidueNumber()<<endl;

			// Get angle of rotation
			double angle = getMinimumAngle(_av, i,_fixedEnd);
			if (angle == MslTools::doubleMax) continue;
			//cout << "ANGLE: "<<angle*180/M_PI<<endl;

			// Rotate the fragment from pivotIndex+1 to end
			rotateFragment(_av, i, angle);
		
			double dist  = _fixedEnd.distance(_av(_av.size()-1));
			//cout << "\n***Dist: "<<dist<<endl<<endl;;
			if (dist < 0.02){
				converged = true;
				break;
			}
		}
		if (converged) break;

		for (uint i = 0; i < _av.size()-2;i++){
			//cout << "Pivot index: "<<i<<" "<<_av(i).getResidueNumber()<<endl;

			// Get angle of rotation
			double angle = getMinimumAngle(_av, i,_fixedEnd);
			if (angle == MslTools::doubleMax) continue;
			//cout << "ANGLE: "<<angle*180/M_PI<<endl;

			// Rotate the fragment from pivotIndex+1 to end
			rotateFragment(_av, i, angle);
		
			double dist  = _fixedEnd.distance(_av(_av.size()-1));
			//cout << "\n***Dist: "<<dist<<endl<<endl;;
			if (dist < 0.02){
				converged = true;
				break;
			}
		}

		numIterations++;
		if (converged) break;
		if (numIterations > 1000000) break;
	}

	if (converged) {
		fprintf(stdout, "Fragment closure in %6d\n", numIterations);
	}
}




double CCD::getMinimumAngle(AtomPointerVector &_av, int _indexOfPivot, Atom &_fixedEnd){
	
	// Get Pid and Pih
	CartesianPoint pid  = (_fixedEnd.getCoor()            - _av(_indexOfPivot).getCoor());
	CartesianPoint pih  = (_av(_av.size()-1).getCoor()    - _av(_indexOfPivot).getCoor());
	
	CartesianPoint pivotBondVector =  (_av(_indexOfPivot+1).getCoor() - _av(_indexOfPivot).getCoor());
	pivotBondVector = pivotBondVector.getUnit();

	// Compute Ks
	double k1 = (pid * pivotBondVector);
	k1       *= (pih * pivotBondVector);
	double k2 = (pid * pih);
	CartesianPoint cross = pivotBondVector.cross(pih);
	double k3 = (pid * cross);

//	cout << "Ks: "<<k1<<" "<<k2<<" "<<k3<<endl;
//	cout << "Pid: "<<pid.toString()<<endl;
//	cout << "Pih: "<<pih.toString()<<endl;
//	cout << "z: "<<pivotBondVector.toString()<<endl;

	// Compute 1st derivative at 0.
	double psi = atan(k3 / (k2 - k1));

	// Compute second derivative
	double secondDerviative = (k1-k2)*cos(psi) - k3* sin(psi);
	//cout << "Second derviative: "<<secondDerviative<<endl;
	if (secondDerviative >= 0) {
		cout << "Failed second derivate test"<<endl;
		return MslTools::doubleMax;
	}


	double distToMaxOfPsi   = k1 * (1 - cos(psi))      + k2 * cos(psi)      + k3 * sin(psi);
	double distToMaxOfPsiPi = k1 * (1 - cos(psi+M_PI)) + k2 * cos(psi+M_PI) + k3 * sin(psi+M_PI);
	

	if (distToMaxOfPsi > distToMaxOfPsiPi){
		return psi;
	}

	return psi+M_PI;
	
}


void CCD::rotateFragment(AtomPointerVector &_av, int _indexOfPivot, double _angleOfRotation){


	// Axis of rotation
	CartesianPoint axisOfRotation = _av(_indexOfPivot+1).getCoor();

	// Transform each atom downstream.
	Transforms t;
	for (uint i = _indexOfPivot+2; i < _av.size();i++){
		t.rotate(_av(i),_angleOfRotation, axisOfRotation,_av(_indexOfPivot).getCoor());
	}
	
}
