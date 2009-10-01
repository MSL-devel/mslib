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

#include "Transforms.h"

Transforms::Transforms() {
	frame["O"] = CartesianPoint(0.0, 0.0, 0.0);
	frame["X"] = CartesianPoint(1.0, 0.0, 0.0);
	frame["Y"] = CartesianPoint(0.0, 1.0, 0.0);

	lastRotMatrix = Matrix(3, 3, 0.0);
	lastTranslation = CartesianPoint(0,0,0);
	saveHistory_flag = false;
}

Transforms::Transforms(const Transforms & _transform) {
	saveHistory_flag = _transform.saveHistory_flag;
	frame = _transform.frame;
	lastRotMatrix = _transform.lastRotMatrix;
	lastTranslation = _transform.lastTranslation;
}

Transforms::~Transforms() {
}

void Transforms::translate(Atom & _atom, const CartesianPoint & _p) {
	_atom.getCoor() += _p;
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second += _p;
		}
	}
	lastTranslation = _p;
}

void Transforms::Xrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getXRotationMatrix(_degrees);
	_atom.getCoor() *= lastRotMatrix;
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Yrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getYRotationMatrix(_degrees);
	_atom.getCoor() *= lastRotMatrix;
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Zrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getZRotationMatrix(_degrees);
	_atom.getCoor() *= lastRotMatrix;
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::rotate(Atom & _atom, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter) {
	Matrix m = CartesianGeometry::instance()->getRotationMatrix(_degrees, _axisFromRotCenter - _rotCenter);
	rotate(_atom, m, _rotCenter);
}

void Transforms::rotate(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter) {
	_atom.getCoor() -= _rotCenter;
	_atom.getCoor() *= _rotMatrix;
	_atom.getCoor() += _rotCenter;
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= _rotCenter;
			k->second *= _rotMatrix;
			k->second += _rotCenter;
		}
	}
	lastRotMatrix = _rotMatrix;
}

void Transforms::align(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {
	if (align(_atom.getCoor(), _target, _rotCenter)) {
		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _rotCenter;
				k->second *= lastRotMatrix;
				k->second += _rotCenter;
			}
		}
	}
}

void Transforms::orient(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	if (orient(_atom.getCoor(), _target, _axis1, _axis2)) {
		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _axis1;
				k->second *= lastRotMatrix;
				k->second += _axis1;
			}
		}
	}
}

void Transforms::translate(AtomVector & _atoms, CartesianPoint _p) {
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() += _p;
	} 
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second += _p;
		}
	}
	lastTranslation = _p;
	/*
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
			*(*m) += _p;
		}
	} 
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second += _p;
		}
	}
	*/
}

void Transforms::Xrotate(AtomVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getXRotationMatrix(_degrees);
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() *= lastRotMatrix;
	}
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Yrotate(AtomVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getYRotationMatrix(_degrees);
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() *= lastRotMatrix;
	}
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Zrotate(AtomVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::instance()->getZRotationMatrix(_degrees);
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() *= lastRotMatrix;
	}
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::rotate(AtomVector & _atoms, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter) {
	Matrix m = CartesianGeometry::instance()->getRotationMatrix(_degrees, _axisFromRotCenter - _rotCenter);
	rotate(_atoms, m, _rotCenter);
}

void Transforms::rotate(AtomVector & _atoms, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter) {
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() -= _rotCenter;
		(*k)->getCoor() *= _rotMatrix;
		(*k)->getCoor() += _rotCenter;
	}
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= _rotCenter;
			k->second *= _rotMatrix;
			k->second += _rotCenter;
		}
	}
	lastRotMatrix = _rotMatrix;
	/*
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
			*(*m) -= _rotCenter;
			*(*m) *= _rotMatrix;
			*(*m) += _rotCenter;
		}
	}
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= _rotCenter;
			k->second *= _rotMatrix;
			k->second += _rotCenter;
		}
	}
	*/
}


void Transforms::align(AtomVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {
	CartesianPoint refCopy(_reference);
	if (align(refCopy, _target, _rotCenter)) {
		for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() -= _rotCenter;
			(*k)->getCoor() *= lastRotMatrix;
			(*k)->getCoor() += _rotCenter;
		}

		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _rotCenter;
				k->second *= lastRotMatrix;
				k->second += _rotCenter;
			}
		}
	}
}

bool Transforms::align(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {

	/***********************************************************
	 *  This function aligns the vector to the target vector.
	 *
	 *  NOTE: the definition of rotation of a dihedral angle and
	 *         the rotation matrix appear to be opposite, that is why
	 *         the rotation around axis is the same value 'angle'
	 *         of the dihedral, and not '-angle'
	 ***********************************************************/
	lastRotMatrix = Matrix(3, 3, 0.0);
	lastRotMatrix[0][0] = 1;
	lastRotMatrix[1][1] = 1;
	lastRotMatrix[2][2] = 1;
	
	CartesianPoint newTarget = _target - _rotCenter;
	if (newTarget.length() == 0) {
		return false;
	}

	_object -= _rotCenter;
	if (_object.length() == 0) {
		return false;
	}

	CartesianPoint cx = _object.cross(newTarget);

	/********* SPECIAL CASES *********************/
	if (cx.length() == 0) {
		// they are parallel
		if (_object * newTarget < 0.0) {
			// they are anti-parallel
			lastRotMatrix[0][0] = -1;
			lastRotMatrix[1][1] = -1;
			lastRotMatrix[2][2] = -1;
		} else {
			// they are parallel, nothing to do
			_object += _rotCenter;
			return false;
		}
	} else {
		/********* NORMAL CASES **********************
		 * Use the cross product to define the axis or rotation
		 *********************************************/
		double degrees = CartesianGeometry::instance()->angle(_object, newTarget);
		lastRotMatrix = CartesianGeometry::instance()->getRotationMatrix(degrees, cx);
	}
	_object *= lastRotMatrix;
	_object += _rotCenter;
	return true;


}

void Transforms::orient(AtomVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	CartesianPoint refCopy(_reference);
	if (orient(refCopy, _target, _axis1, _axis2)) {
		for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() -= _axis1;
			(*k)->getCoor() *= lastRotMatrix;
			(*k)->getCoor() += _axis1;
		}

		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _axis1;
				k->second *= lastRotMatrix;
				k->second += _axis1;
			}
		}
	}
}
bool Transforms::orient(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {

	/****************************************
	 *  This function axially rotates the point around
	 *  a line ('center' + n * 'axis') to orient 
	 *  it to the same direction from the line with respect 
	 *  to the reference point 'target' so that the
	 *  dihedral is = 0
	 * 
	 *
	 *   _object                _target
	 *          \              /
	 *          _axis1---_axis2
	 *
	 *  NOTE2: the definition of rotation of a dihedral angle and
	 *         the rotation matrix appear to be opposite, that is why
	 *         the rotation around axis is the same value 'angle'
	 *         of the dihedral, and not '-angle'
	 ****************************************/

	if (_axis2 == _axis1 || _axis2 == _target || _axis1 == _target) {
		lastRotMatrix = Matrix(3, 3, 0.0);
		lastRotMatrix[0][0] = 1;
		lastRotMatrix[1][1] = 1;
		lastRotMatrix[2][2] = 1;
		return false;
	}
	double degrees = CartesianGeometry::instance()->dihedral(_object, _axis1, _axis2, _target);
	lastRotMatrix = CartesianGeometry::instance()->getRotationMatrix(degrees, _axis2 - _axis1);
	_object -= _axis1;
	_object *= lastRotMatrix;
	_object += _axis1;

	return true;

}

/*
double Transforms::align(vector<Residue *> &_align, vector<Residue *> &_ref){
	vector<Residue *> tmp;
	return align(_align,_ref,tmp);
}
double Transforms::align(vector<Residue *> &_align, vector<Residue *> &_ref, vector<Residue *> &_move){

	AtomVector all;
	if (_move.size() != 0){
		for (uint i = 0; i < _move.size();i++){
			all += _move[i]->getAtoms();
		}
	}

	return align(_align,_ref,all);

}
double Transforms::align(vector<Residue *> &_align, vector<Residue *> &_ref, AtomVector &_move){

	if (_align.size() != _ref.size()){
		return 999;
	}

	bool emptyMove = false;
	if (_move.size() == 0){
		emptyMove = true;
	}
	AtomVector align;
	AtomVector ref;
	for (uint i = 0; i < _align.size();i++){
		align.push_back(new Atom(_align[i]->getAtom("N")));
		align.push_back(new Atom(_align[i]->getAtom("CA")));
		align.push_back(new Atom(_align[i]->getAtom("C")));

 		ref.push_back(new Atom(_ref[i]->getAtom("N")));
 		ref.push_back(new Atom(_ref[i]->getAtom("CA")));
 		ref.push_back(new Atom(_ref[i]->getAtom("C")));

		if (emptyMove){
			_move += _align[i]->getAtoms();
		}
	}

//	cout << "Number of atoms: "<<align.size()<<" move "<<_move.size()<<" atoms."<<endl;
	if (align.size() != ref.size()){
		align.deletePointers();
		ref.deletePointers();
		return 999;
	}
// 	cout << "Moving:"<<endl;
// 	for (uint i = 0; i < _move.size();i++){
// 		cout << _move(i)<<endl;
// 	}
	if (!this->align(align,ref,_move)) {
		align.deletePointers();
		ref.deletePointers();
		return 999;
	}

// 	cout << "Moved:"<<endl;
// 	for (uint i = 0; i < _move.size();i++){
// 		cout << _move(i)<<endl;
// 	}

	AtomVector resultAlign;
	AtomVector resultRef;
	for (uint i = 0; i < _align.size();i++){
		resultAlign.push_back(&_align[i]->getAtom("N"));
		resultAlign.push_back(&_align[i]->getAtom("CA"));
		resultAlign.push_back(&_align[i]->getAtom("C"));

 		resultRef.push_back(&_ref[i]->getAtom("N"));
 		resultRef.push_back(&_ref[i]->getAtom("CA"));
 		resultRef.push_back(&_ref[i]->getAtom("C"));
	}

	
	double rmsd =  resultAlign.rmsd(resultRef);

	resultAlign.clear();
	resultRef.clear();
	if (emptyMove){
		_move.clear();
	}

	align.deletePointers();
	ref.deletePointers();
	return rmsd;

}
*/

/*
bool Transforms::align(AtomVector &_align, AtomVector &_ref, AtomVector &_moveable){

	if (_align.size() == 0) {
		cerr << "ERROR Transforms::align() ; align vector is emtpy\n";
		return false;
	}
		
	if (_ref.size() == 0){
		cerr << "ERROR Transforms::align() ; ref vector is emtpy\n";
		return false;
	}
	if (_moveable.size() == 0){
		cerr << "ERROR Transforms::align() ; moveable vector is emtpy\n";
		return false;
	}

	// Create a quaternion that minimizes the RMSD between two sets of points (COM1, COM2 center-of-mass get defined as well)
	if (!q.makeQuaternion(_align,_ref)) return false;

	Matrix rotationMatrix(3,3,0.0);
	q.convertToRotationMatrix(rotationMatrix);

	_ref.updateGeometricCenter();
	_align.updateGeometricCenter();
	CartesianPoint GC1 = _ref.getGeometricCenter();
	CartesianPoint GC2 = _align.getGeometricCenter();
	
	
	GC2 *= -1;
	
	for (uint i = 0; i< _moveable.size();i++){

		// Translate moveable to GC2
		//_moveable[i][0] -= GC2[0];  _moveable[i][1] -= GC2[1];  _moveable[i][2] -= GC2[2];
		translate(_moveable(i), GC2); // GC2->origin

		double tmpx,tmpy,tmpz;
 		tmpx = (_moveable(i)[0] * rotationMatrix[0][0]) + (_moveable(i)[1] * rotationMatrix[1][0]) + (_moveable(i)[2] * rotationMatrix[2][0]);
 		tmpy = (_moveable(i)[0] * rotationMatrix[0][1]) + (_moveable(i)[1] * rotationMatrix[1][1]) + (_moveable(i)[2] * rotationMatrix[2][1]);
 		tmpz = (_moveable(i)[0] * rotationMatrix[0][2]) + (_moveable(i)[1] * rotationMatrix[1][2]) + (_moveable(i)[2] * rotationMatrix[2][2]);
		_moveable(i).setCoor(tmpx,tmpy,tmpz);

		// Translate _moveable to GC1
		//_moveable[i][0] += GC1[0];  _moveable[i][1] += GC1[1];  _moveable[i][2] += GC1[2];
		translate(_moveable(i), GC1); // origin->GC1

	}


	lastRotMatrix = rotationMatrix;
	return true;
}
*/

bool Transforms::align(AtomVector &_align, AtomVector &_ref, AtomVector &_moveable){

	// Create a quaternion that minimizes the RMSD between two sets of points (COM1, COM2 center-of-mass get defined as well)
	if (!q.makeQuaternion(_align,_ref)) return false;

	//Matrix rotationMatrix(3,3,0.0);
	q.convertToRotationMatrix(lastRotMatrix);

	_ref.updateGeometricCenter();
	_align.updateGeometricCenter();
	CartesianPoint GC1 = _ref.getGeometricCenter();
	CartesianPoint GC2 = _align.getGeometricCenter();
	
	
	GC2 *= -1;
	
	for (uint i = 0; i< _moveable.size();i++){

		// Translate moveable to GC2
		//_moveable[i][0] -= GC2[0];  _moveable[i][1] -= GC2[1];  _moveable[i][2] -= GC2[2];
		translate(_moveable(i), GC2); // GC2->origin

		double tmpx,tmpy,tmpz;
 		tmpx = (_moveable(i)[0] * lastRotMatrix[0][0]) + (_moveable(i)[1] * lastRotMatrix[1][0]) + (_moveable(i)[2] * lastRotMatrix[2][0]);
 		tmpy = (_moveable(i)[0] * lastRotMatrix[0][1]) + (_moveable(i)[1] * lastRotMatrix[1][1]) + (_moveable(i)[2] * lastRotMatrix[2][1]);
 		tmpz = (_moveable(i)[0] * lastRotMatrix[0][2]) + (_moveable(i)[1] * lastRotMatrix[1][2]) + (_moveable(i)[2] * lastRotMatrix[2][2]);
		_moveable(i).setCoor(tmpx,tmpy,tmpz);

		// Translate _moveable to GC1
		//_moveable[i][0] += GC1[0];  _moveable[i][1] += GC1[1];  _moveable[i][2] += GC1[2];
		translate(_moveable(i), GC1); // origin->GC1

	}

	lastTranslation = GC1 - GC2;
	return true;
}


/*
Matrix createBasisTransformation(vector<vector<double> > &_basis1, vector<vector<double> > &_basis2){
	Matrix m;
	m.initialize(3,3);

	CartesianPoint i,j,k,ip,jp,kp;
	i.setCoor(_basis1[0]);
	j.setCoor(_basis1[1]);
	k.setCoor(_basis1[2]);

	ip.setCoor(_basis2[0]);
	jp.setCoor(_basis2[1]);
	kp.setCoor(_basis2[2]);
    
	//  Transformation Matrix going from basis1 to basis2
	// Dot products
	m[0][0] = i * ip;
	m[0][1] = i * jp;
	m[0][2] = i * kp;
  
	m[1][0] = j * ip;
	m[1][1] = j * jp;
	m[1][2] = j * kp;
  
	m[2][0] = k * ip;
	m[2][1] = k * jp;
	m[2][2] = k * kp;

	return m;
}
*/

Matrix Transforms::createBasisTransformation(vector<vector<double> > &_basis1, vector<vector<double> > &_basis2){
	//Matrix m;
	//m.initialize(3,3);

	CartesianPoint i,j,k,ip,jp,kp;
	i.setCoor(_basis1[0]);
	j.setCoor(_basis1[1]);
	k.setCoor(_basis1[2]);

	ip.setCoor(_basis2[0]);
	jp.setCoor(_basis2[1]);
	kp.setCoor(_basis2[2]);
    
	//  Transformation Matrix going from basis1 to basis2
	// Dot products
	lastRotMatrix[0][0] = i * ip;
	lastRotMatrix[0][1] = i * jp;
	lastRotMatrix[0][2] = i * kp;
  
	lastRotMatrix[1][0] = j * ip;
	lastRotMatrix[1][1] = j * jp;
	lastRotMatrix[1][2] = j * kp;
  
	lastRotMatrix[2][0] = k * ip;
	lastRotMatrix[2][1] = k * jp;
	lastRotMatrix[2][2] = k * kp;

	return lastRotMatrix;
}

bool Transforms::align(AtomVector &_align, AtomVector &_ref){
	return align(_align,_ref,_align);
}


void Transforms::applyHistory(Atom & _atom) {


	CartesianPoint tX = frame["O"] + CartesianPoint(1.0, 0.0, 0.0);
	align(tX, frame["X"], frame["O"]);
	_atom.getCoor() *= lastRotMatrix;

	CartesianPoint tY = frame["O"] + (CartesianPoint(0.0, 1.0, 0.0) * lastRotMatrix);
	orient(tY, frame["Y"], frame["O"], frame["X"]);
	_atom.getCoor() *= lastRotMatrix;
	_atom.getCoor() += frame["O"];
}

void Transforms::applyHistory(AtomVector & _atoms) {
	CartesianPoint tX = frame["O"] + CartesianPoint(1.0, 0.0, 0.0);
	align(tX, frame["X"], frame["O"]);
	Matrix rot1 = lastRotMatrix;

	CartesianPoint tY = frame["O"] + (CartesianPoint(0.0, 1.0, 0.0) * lastRotMatrix);
	orient(tY, frame["Y"], frame["O"], frame["X"]);
	Matrix rot2 = lastRotMatrix * rot1;
	for (AtomVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() *= rot2;
		(*k)->getCoor() += frame["O"];
	}
}



SphericalPoint Transforms::transform(CartesianPoint &_p1){

	double radius = sqrt(_p1.getX() *_p1.getX() + _p1.getY() *_p1.getY() + _p1.getZ() *_p1.getZ());
	double sigma  = atan2(_p1.getY(), _p1.getX());
	//double theta  = atan( sqrt(_p1.getX()*_p1.getX() + _p1.getY()*_p1.getY()) / _p1.getZ());
	double theta  = acos(_p1.getZ() / radius);


	SphericalPoint sp(radius,sigma,theta);

	return sp;
}

void Transforms::RotatePdbAboutZYX(AtomVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _localZ, CartesianPoint & _localY, CartesianPoint & _localX, double _RotationAlongZ, double _RotationAlongY, double _RotationAlongX) {
	Matrix zrot = CartesianGeometry::instance()->getRotationMatrix(_RotationAlongZ, _localZ);
	Matrix yrot = CartesianGeometry::instance()->getRotationMatrix(_RotationAlongY, _localY);
	Matrix xrot = CartesianGeometry::instance()->getRotationMatrix(_RotationAlongX, _localX);
	Matrix combination = zrot*yrot*xrot;

	rotate(_theAtoms,combination,_center);
	_localZ *= combination;
	_localY *= combination;
	_localX *= combination;
}

void Transforms::TranslateRigidBodyPdbResidue(AtomVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _TranslationVector, double _TranslationAmount) {
	CartesianPoint translation = _TranslationVector*_TranslationAmount;
	translate(_theAtoms,translation);
	_center += translation;
}

