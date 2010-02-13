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

#include "AtomPointerVector.h"

#include "CartesianGeometry.h"

AtomPointerVector::AtomPointerVector(){
	name = "";
	archiveType = "binary";
}

AtomPointerVector::AtomPointerVector(unsigned int _size, Atom * _pointer) {
	if (_size > 0) {
		this->insert(this->begin(), _size, _pointer);
	}
}

AtomPointerVector::AtomPointerVector(const AtomPointerVector & _atoms){
	name = _atoms.name;
	archiveType = _atoms.archiveType;
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		this->push_back(*k);
	}
}

void AtomPointerVector::operator=(const AtomPointerVector & _atoms){
	name = _atoms.name;
	archiveType = _atoms.archiveType;
	this->clear();
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		this->push_back(*k);
	}
}

AtomPointerVector AtomPointerVector::operator+(const AtomPointerVector & _atoms) {
	AtomPointerVector out = *this;
	out += _atoms;
	return out;
}

void AtomPointerVector::operator+=(const AtomPointerVector & _atoms) {
	// adds two atom vectors (but does not add common elements)
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		bool found = false;
		for (AtomPointerVector::iterator l=begin(); l!=end(); l++) {
			if (*k == *l) {
				found = true;
				break;
			}
		}
		if (!found) {
			this->push_back(*k);
		}
	}
}

AtomPointerVector AtomPointerVector::operator-(const AtomPointerVector & _atoms) {
	AtomPointerVector out = *this;
	out -= _atoms;
	return out;
}

void AtomPointerVector::operator-=(const AtomPointerVector & _atoms) {
	// removes all elements belonging to _atoms that are present in this
	for (AtomPointerVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		for (AtomPointerVector::iterator l=begin(); l!=end(); l++) {
			if (*k == *l) {
				erase(l);
				l--;
			}
		}
	}
}

AtomPointerVector::~AtomPointerVector(){
}


void AtomPointerVector::updateGeometricCenter(unsigned int _updateStamp){

	/**************************************************
	 *   An optional integer argument allow not to
	 *   recompute the center if it was recently done,
	 *   as determined by the identity of the stamp
	 **************************************************/

	if (updateStamp != 0 && updateStamp == _updateStamp) {
		// it was already updated
		return;
	}
	updateStamp = _updateStamp;
	geometricCenter.setCoor(0.0, 0.0, 0.0);
	
	for (uint i =0; i < (*this).size();i++){
		geometricCenter += (*this)[i]->getCoor();
	}
	geometricCenter /= (double)(*this).size();

}



double AtomPointerVector::rmsd(const AtomPointerVector &_av) const {

	// Definitions
	double sumR = 0.0;
	//int counter;

	// Initializations
	//counter =0; 
	//sumR = 0.0;

	if (size() != _av.size()) {
		cerr << "ERROR 38919: different number of atoms in rmsd calculation (" << size() << " != " << _av.size() << " in double AtomPointerVector::rmsd(const AtomPointerVector &_av) const" << endl;
		exit(38919);
		//return 999;
	}

	if (size() == 0) {
		return 0.0;
	}

	for (unsigned int i = 0; i < size(); i++){
		sumR += pow((*this)[i]->distance(*_av[i]), 2.0);
	}
	// Loop..
//	for (uint i = 0; i < (*this).size(); i++){
//	
//		// Calculations
//		double dist = sqrt( 
//		(((*_av[i])[0] - (*(*this)[i])[0])  * ((*_av[i])[0] - (*(*this)[i])[0])) + 
//		(((*_av[i])[1] - (*(*this)[i])[1])  * ((*_av[i])[1] - (*(*this)[i])[1])) + 
//		(((*_av[i])[2] - (*(*this)[i])[2])  * ((*_av[i])[2] - (*(*this)[i])[2])));  
//	
//		sumR += dist*dist;
//		counter++;
//	}

//<<<<<<< AtomPointerVector.cpp
	return sqrt(sumR/size());
//=======
//	 // Initializations
//	 counter =0; sumR = 0.0;
//
//	 if ((*this).size() != _av.size()) {
//	   return 999;
//	 }
//
//	 // Loop..
//	 for (uint i = 0; i < (*this).size(); i++){
//
//	   // Calculations
//	   double dist = sqrt( 
//			      (((*_av[i])[0] - (*(*this)[i])[0])  * ((*_av[i])[0] - (*(*this)[i])[0])) + 
//			      (((*_av[i])[1] - (*(*this)[i])[1])  * ((*_av[i])[1] - (*(*this)[i])[1])) + 
//			      (((*_av[i])[2] - (*(*this)[i])[2])  * ((*_av[i])[2] - (*(*this)[i])[2])));  
//
//	   sumR += dist*dist;
//	   counter++;
//	 }
//
//	 return sqrt(sumR/(counter));
//>>>>>>> 1.10
}

void AtomPointerVector::translate(double _x, double _y, double _z){
	translate(CartesianPoint(_x, _y, _z));
}

void AtomPointerVector::translate(const CartesianPoint &_vec){
	for (uint i = 0 ; i < (*this).size();i++){
		if ((*this)[i]->hasCoor()) {
			(*this)[i]->setCoor((*this)(i).getCoor() + _vec); 
		}
	}
}         


void AtomPointerVector::rotate(const Matrix &_rotMat){

	for (uint i = 0 ; i < (*this).size();i++){

		//(*this)(i).setCoor(CartesianGeometry::instance()->matrixTransposeTimesCartesianPoint((*this)(i).getCoor(),_rotMat)); 
		if ((*this)[i]->hasCoor()) {
			(*this)[i]->setCoor(CartesianGeometry::instance()->matrixTimesCartesianPoint((*this)(i).getCoor(),_rotMat)); 
		}
	}

}
void AtomPointerVector::saveCoor(string _coordName){
	for (uint i = 0; i < (*this).size();i++){
		(*this)(i).saveCoor(_coordName);
	}
}

bool AtomPointerVector::applySavedCoor(string _coordName){
	bool result = 1;
	for (uint i = 0; i < (*this).size();i++){
		result &= (*this)(i).applySavedCoor(_coordName);
	}

	return result;
}

void AtomPointerVector::clearSavedCoor(){
	for (uint i = 0; i < (*this).size();i++){
		(*this)(i).clearSavedCoor();
	}
}


void AtomPointerVector::deletePointers(){
	for (AtomPointerVector::iterator k=begin(); k!=end(); k++) {
		delete *k;
	}
	(*this).clear();
}
