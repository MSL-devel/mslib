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

#include "Transforms.h"

using namespace MSL;
using namespace std;
#include "MslOut.h"
static MslOut MSLOUT("Transforms");

Transforms::Transforms() {

	frame["O"] = CartesianPoint(0.0, 0.0, 0.0);
	frame["X"] = CartesianPoint(1.0, 0.0, 0.0);
	frame["Y"] = CartesianPoint(0.0, 1.0, 0.0);

	lastRotMatrix = Matrix(3, 3, 0.0);
	lastTranslation = CartesianPoint(0,0,0);
	saveHistory_flag = false;
	transformAllCoors_flag = true;
	naturalMovementOnSetDOF_flag = false;
}

Transforms::Transforms(const Transforms & _transform) {
	saveHistory_flag = _transform.saveHistory_flag;
	frame = _transform.frame;
	lastRotMatrix = _transform.lastRotMatrix;
	lastTranslation = _transform.lastTranslation;
	transformAllCoors_flag = _transform.transformAllCoors_flag;
	naturalMovementOnSetDOF_flag = _transform.naturalMovementOnSetDOF_flag;
}

Transforms::~Transforms() {
}

/* -- PRIVATE TRANSFORM FUNCTIONS THAT DO NOT UPDATE THE HISTORY AND LAST TRANSFORM MEMORY -- */
void Transforms::translateAtom(Atom & _atom, const CartesianPoint & _p) {
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) += _p;
		}
	} else {
		_atom.getCoor() += _p;
	}
}

void Transforms::rotateAtom(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter) {
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) -= _rotCenter;
			*(*m) *= _rotMatrix;
			*(*m) += _rotCenter;
		}
	} else {
		_atom.getCoor() -= _rotCenter;
		_atom.getCoor() *= _rotMatrix;
		_atom.getCoor() += _rotCenter;
	}
}

bool Transforms::alignAtom(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {
	if (align(_atom.getCoor(), _target, _rotCenter)) {
		if (transformAllCoors_flag && _atom.getNumberOfAltConformations() > 1) {
			// apply to all alt confs
			unsigned int active = _atom.getActiveConformation();
			vector<CartesianPoint *> & pts = _atom.getAllCoor();
			for (unsigned int i=0; i<pts.size(); i++) {
				// the active was already transformed
				if (i != active) {
					*pts[i] -= _rotCenter;
					*pts[i] *= lastRotMatrix;
					*pts[i] += _rotCenter;
				}
			}
		}
		return true;
	}
	return false;
}

bool Transforms::orientAtom(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	if (orient(_atom.getCoor(), _target, _axis1, _axis2)) {
		if (transformAllCoors_flag && _atom.getNumberOfAltConformations() > 1) {
			// apply to all alt confs
			unsigned int active = _atom.getActiveConformation();
			vector<CartesianPoint *> & pts = _atom.getAllCoor();
			for (unsigned int i=0; i<pts.size(); i++) {
				// the active was already transformed
				if (i != active) {
					*pts[i] -= _axis1;
					*pts[i] *= lastRotMatrix;
					*pts[i] += _axis1;
				}
			}
		}
		return true;
	}
	return false;
}

/* ---------------------- END OF PRIVATE TRANSFORM FUNCTONS ----------------------------- */


void Transforms::translate(Atom & _atom, const CartesianPoint & _p) {
	/*
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) += _p;
		}
	} else {
		_atom.getCoor() += _p;
	}
	*/
	translateAtom(_atom, _p);
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second += _p;
		}
	}
	lastTranslation = _p;
}

void Transforms::Xrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::getXRotationMatrix(_degrees);
	/*
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) *= lastRotMatrix;
		}
	} else {
		_atom.getCoor() *= lastRotMatrix;
	}
	*/
	rotateAtom(_atom, lastRotMatrix);
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Yrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::getYRotationMatrix(_degrees);
	/*
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) *= lastRotMatrix;
		}
	} else {
		_atom.getCoor() *= lastRotMatrix;
	}
	*/
	rotateAtom(_atom, lastRotMatrix);
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Zrotate(Atom & _atom, double _degrees) {
	lastRotMatrix = CartesianGeometry::getZRotationMatrix(_degrees);
	/*
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) *= lastRotMatrix;
		}
	} else {
		_atom.getCoor() *= lastRotMatrix;
	}
	*/
	rotateAtom(_atom, lastRotMatrix);
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::rotate(Atom & _atom, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter) {
	Matrix m = CartesianGeometry::getRotationMatrix(_degrees, _axisFromRotCenter - _rotCenter);
	rotate(_atom, m, _rotCenter);
}

void Transforms::rotate(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter) {
	/*
	if (transformAllCoors_flag) {
		for (vector<CartesianPoint *>::iterator m=_atom.getAllCoor().begin(); m!=_atom.getAllCoor().end(); m++) {
			*(*m) -= _rotCenter;
			*(*m) *= _rotMatrix;
			*(*m) += _rotCenter;
		}
	} else {
		_atom.getCoor() -= _rotCenter;
		_atom.getCoor() *= _rotMatrix;
		_atom.getCoor() += _rotCenter;
	}
	*/
	rotateAtom(_atom, _rotMatrix, _rotCenter);
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= _rotCenter;
			k->second *= _rotMatrix;
			k->second += _rotCenter;
		}
	}
	lastRotMatrix = _rotMatrix;
}

bool Transforms::align(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {
	if (alignAtom(_atom, _target, _rotCenter)) {
		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _rotCenter;
				k->second *= lastRotMatrix;
				k->second += _rotCenter;
			}
		}
		return true;
	}
	return false;
}

bool Transforms::orient(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	if (orientAtom(_atom, _target, _axis1, _axis2)) {
		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _axis1;
				k->second *= lastRotMatrix;
				k->second += _axis1;
			}
		}
		return true;
	}
	return false;
}

void Transforms::translate(AtomPointerVector & _atoms, CartesianPoint _p) {
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		translateAtom(*(*k), _p);
	} 
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second += _p;
		}
	}
	lastTranslation = _p;
	/*
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
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

void Transforms::Xrotate(AtomPointerVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::getXRotationMatrix(_degrees);

	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		rotateAtom(*(*k), lastRotMatrix);
	} 
	/*
	if (transformAllCoors_flag) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
				*(*m) *= lastRotMatrix;
			}
		} 
	} else {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() *= lastRotMatrix;
		}
	}
	*/
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Yrotate(AtomPointerVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::getYRotationMatrix(_degrees);

	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		rotateAtom(*(*k), lastRotMatrix);
	} 
	/*
	if (transformAllCoors_flag) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
				*(*m) *= lastRotMatrix;
			}
		} 
	} else {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() *= lastRotMatrix;
		}
	}
	*/
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::Zrotate(AtomPointerVector & _atoms, double _degrees) {
	lastRotMatrix = CartesianGeometry::getZRotationMatrix(_degrees);

	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		rotateAtom(*(*k), lastRotMatrix);
	} 
	/*
	if (transformAllCoors_flag) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
				*(*m) *= lastRotMatrix;
			}
		} 
	} else {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() *= lastRotMatrix;
		}
	}
	*/
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second *= lastRotMatrix;
		}
	}
}

void Transforms::rotate(AtomPointerVector & _atoms, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter) {
	Matrix m = CartesianGeometry::getRotationMatrix(_degrees, _axisFromRotCenter - _rotCenter);
	rotate(_atoms, m, _rotCenter);
}

void Transforms::rotate(AtomPointerVector & _atoms, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter) {
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		rotateAtom(*(*k), _rotMatrix, _rotCenter);
	} 
	/*
	if (transformAllCoors_flag) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
				*(*m) -= _rotCenter;
				*(*m) *= _rotMatrix;
				*(*m) += _rotCenter;
			}
		} 
	} else {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			(*k)->getCoor() -= _rotCenter;
			(*k)->getCoor() *= _rotMatrix;
			(*k)->getCoor() += _rotCenter;
		}
	}
	*/
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= _rotCenter;
			k->second *= _rotMatrix;
			k->second += _rotCenter;
		}
	}
	lastRotMatrix = _rotMatrix;
	/*
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
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


bool Transforms::align(AtomPointerVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _rotCenter) {
	CartesianPoint refCopy(_reference);
	if (align(refCopy, _target, _rotCenter)) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			rotateAtom(*(*k), lastRotMatrix, _rotCenter);
		}
		/*
		if (transformAllCoors_flag) {
			for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
				for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
					*(*m) -= _rotCenter;
					*(*m) *= lastRotMatrix;
					*(*m) += _rotCenter;
				}
			} 
		} else {
			for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
				(*k)->getCoor() -= _rotCenter;
				(*k)->getCoor() *= lastRotMatrix;
				(*k)->getCoor() += _rotCenter;
			}
		}

		*/
		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _rotCenter;
				k->second *= lastRotMatrix;
				k->second += _rotCenter;
			}
		}
		return true;
	}
	return false;
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
		double degrees = CartesianGeometry::angle(_object, newTarget);
		lastRotMatrix = CartesianGeometry::getRotationMatrix(degrees, cx);
	}
	_object *= lastRotMatrix;
	_object += _rotCenter;
	return true;


}

bool Transforms::orient(AtomPointerVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	CartesianPoint refCopy(_reference);
	if (orient(refCopy, _target, _axis1, _axis2)) {
		for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
			rotateAtom(*(*k), lastRotMatrix, _axis1);
		}
		/*
		if (transformAllCoors_flag) {
			for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
				for (vector<CartesianPoint *>::iterator m=(*k)->getAllCoor().begin(); m!=(*k)->getAllCoor().end(); m++) {
					*(*m) -= _axis1;
					*(*m) *= lastRotMatrix;
					*(*m) += _axis1;
				}
			} 
		} else {
			for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
				(*k)->getCoor() -= _axis1;
				(*k)->getCoor() *= lastRotMatrix;
				(*k)->getCoor() += _axis1;
			}
		}

		if (saveHistory_flag) {
			for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
				k->second -= _axis1;
				k->second *= lastRotMatrix;
				k->second += _axis1;
			}
		}
		*/
		return true;
	}
	return false;
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
	double degrees = CartesianGeometry::dihedral(_object, _axis1, _axis2, _target);
	lastRotMatrix = CartesianGeometry::getRotationMatrix(degrees, _axis2 - _axis1);
	_object -= _axis1;
	_object *= lastRotMatrix;
	_object += _axis1;

	return true;

}

bool Transforms::rmsdAlignment(AtomPointerVector &_align, AtomPointerVector &_ref){
	return rmsdAlignment(_align,_ref,_align);
}

bool Transforms::rmsdAlignment(AtomPointerVector &_align, AtomPointerVector &_ref, AtomPointerVector &_moveable){

	// Create a quaternion that minimizes the RMSD between two sets of points (COM1, COM2 center-of-mass get defined as well)
	if (!q.makeQuaternion(_align,_ref)) return false;

	//Matrix rotationMatrix(3,3,0.0);
	q.convertToRotationMatrix(lastRotMatrix);

	//_ref.updateGeometricCenter();
	//_align.updateGeometricCenter();
	CartesianPoint GC1 = _ref.getGeometricCenter();
	CartesianPoint GC2 = _align.getGeometricCenter();
	CartesianPoint pt = _ref.getGeometricCenter() - _align.getGeometricCenter();
	
//	GC2 *= -1;
	
	Matrix rotMatrix = lastRotMatrix.getTranspose();

	for (AtomPointerVector::iterator k=_moveable.begin(); k!=_moveable.end(); k++) {
//		translateAtom(*(*k), GC2); // GC2->origin
//		rotateAtom(*(*k), rotMatrix);
//		translateAtom(*(*k), GC1); // origin->GC1
		translateAtom(*(*k), pt); // GC2->origin
		rotateAtom(*(*k), rotMatrix, GC1);
	}
	/*
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
	*/
	if (saveHistory_flag) {
		for (map<string, CartesianPoint>::iterator k=frame.begin(); k!=frame.end(); k++) {
			k->second -= GC2;
			k->second *= rotMatrix;
			k->second += GC1;
		}
	}

	lastTranslation = GC1 - GC2;


	// Output when Transforms output is turned on.
	MSLOUT.stream() << "Rotation Matrix: "<<lastRotMatrix.toString()<<endl;
	MSLOUT.stream() << "Translation Vector: "<<lastTranslation.toString()<<endl;

	return true;
}

bool Transforms::revertRmsdAlignment(AtomPointerVector &_align, AtomPointerVector &_ref, AtomPointerVector &_moveable){
	CartesianPoint GC1 = _ref.getGeometricCenter();
	CartesianPoint GC2 = _align.getGeometricCenter();
	CartesianPoint oldGC2 = GC1 - lastTranslation;
	CartesianPoint pt = oldGC2 - GC2;

	for (AtomPointerVector::iterator k = _moveable.begin(); k != _moveable.end(); k++) {
		translateAtom(*(*k), pt);
		rotateAtom(*(*k), lastRotMatrix, oldGC2);
	}

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


void Transforms::applyHistory(Atom & _atom) {

	/*
	CartesianPoint tX = frame["O"] + CartesianPoint(1.0, 0.0, 0.0);
	align(tX, frame["X"], frame["O"]);
	_atom.getCoor() *= lastRotMatrix;

	CartesianPoint tY = frame["O"] + (CartesianPoint(0.0, 1.0, 0.0) * lastRotMatrix);
	orient(tY, frame["Y"], frame["O"], frame["X"]);
	_atom.getCoor() *= lastRotMatrix;
	_atom.getCoor() += frame["O"];
	*/

	CartesianPoint tX = frame["O"] + CartesianPoint(1.0, 0.0, 0.0);
	align(tX, frame["X"], frame["O"]);
	Matrix rot1 = lastRotMatrix;

	CartesianPoint tY = frame["O"] + (CartesianPoint(0.0, 1.0, 0.0) * lastRotMatrix);
	orient(tY, frame["Y"], frame["O"], frame["X"]);
	Matrix rot2 = lastRotMatrix * rot1;

	rotateAtom(_atom, rot2);
	translateAtom(_atom, frame["O"]);

}

void Transforms::applyHistory(AtomPointerVector & _atoms) {
	CartesianPoint tX = frame["O"] + CartesianPoint(1.0, 0.0, 0.0);
	align(tX, frame["X"], frame["O"]);
	Matrix rot1 = lastRotMatrix;

	CartesianPoint tY = frame["O"] + (CartesianPoint(0.0, 1.0, 0.0) * lastRotMatrix);
	orient(tY, frame["Y"], frame["O"], frame["X"]);
	Matrix rot2 = lastRotMatrix * rot1;
	/*
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		(*k)->getCoor() *= rot2;
		(*k)->getCoor() += frame["O"];
	}
	*/
	for (AtomPointerVector::iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		rotateAtom(*(*k), rot2);
		translateAtom(*(*k), frame["O"]);
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

/*
void Transforms::RotatePdbAboutZYX(AtomPointerVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _localZ, CartesianPoint & _localY, CartesianPoint & _localX, double _RotationAlongZ, double _RotationAlongY, double _RotationAlongX) {
	Matrix zrot = CartesianGeometry::getRotationMatrix(_RotationAlongZ, _localZ);
	Matrix yrot = CartesianGeometry::getRotationMatrix(_RotationAlongY, _localY);
	Matrix xrot = CartesianGeometry::getRotationMatrix(_RotationAlongX, _localX);
	Matrix combination = zrot*yrot*xrot;

	rotate(_theAtoms,combination,_center);
	_localZ *= combination;
	_localY *= combination;
	_localX *= combination;
}

void Transforms::TranslateRigidBodyPdbResidue(AtomPointerVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _TranslationVector, double _TranslationAmount) {
	CartesianPoint translation = _TranslationVector*_TranslationAmount;
	translate(_theAtoms,translation);
	_center += translation;
}
*/

bool Transforms::setBondDistance(Atom & _atom1, Atom & _atom2, double _distance) {
	/* move _atom2 in the same direction of the _atom1-_atom2 bond (they do not
	   need to be bonded through) to place them at the right distance, and move 
	   any other atom connected to _atom2 (bonded and those bonded that follow) as 
	   well, if the bonded information is stored in the atoms */
	double d =  _distance - _atom1.distance(_atom2);
	CartesianPoint traslation = (_atom2.getCoor() - _atom1.getCoor()).getUnit() * d;

	// build the list of atoms to be moved
	set<Atom*> excludeList;
	excludeList.insert(&_atom1);
	set<Atom*> moveList = _atom2.findLinkedAtoms(excludeList);
	moveList.insert(&_atom2);
	set<Atom*> secondMoveList; // used only for natural rotations
	if (naturalMovementOnSetDOF_flag) {
		excludeList.clear();
		excludeList.insert(&_atom2);
		secondMoveList = _atom1.findLinkedAtoms(excludeList);
		secondMoveList.insert(&_atom1);
	}

	if (!naturalMovementOnSetDOF_flag) {
		// move the atoms
		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			translateAtom(**k, traslation); 
		}
	} else {
		double f1 = (double)moveList.size();
		double f2 = (double)secondMoveList.size();
		
		// normalize
		double tot = f1 + f2;
		f1 /= -tot;
		f2 /= tot;

		// get the translations
		CartesianPoint t1 = traslation * f2;
		CartesianPoint t2 = traslation * f1;

		cout << traslation << " " << f1 << " " << t1 << " " << f2 << " " << t2 << endl;
		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			translateAtom(**k, t1); 
		}
		for (set<Atom*>::iterator k=secondMoveList.begin(); k!=secondMoveList.end(); k++) {
			translateAtom(**k, t2); 
		}
	}
	return true;
}

bool Transforms::setBondAngle(Atom & _atom1, Atom & _atom2, Atom & _atom3, double _angleDegrees) {
	/* rotate _atom3 in the same direction of the _atom1-_atom2 bond (they do not
	   need to be bonded through) to place the angle at the right value, and move 
	   any other atom connected to _atom3 (bonded and those bonded that follow) as 
	   well, if the bonded information is stored in the atoms */

	// bring the angle into the -180 to 180 range and then make it positive
	while (_angleDegrees < -180.0) {
		_angleDegrees += 360.0;
	}
	while (_angleDegrees > 180.0) {
		_angleDegrees -= 360.0;
	}
	if (_angleDegrees < 0.0) {
		_angleDegrees *= -1.0;
	}

	// bring the measured angle into the -180 to 180 range and then make it positive
	double currentAngle = _atom1.angle(_atom2, _atom3);
	while (currentAngle < -180.0) {
		currentAngle += 360.0;
	}
	while (currentAngle > 180.0) {
		currentAngle -= 360.0;
	}
	if (currentAngle < 0.0) {
		currentAngle *= -1.0;
	}
	double rotation = _angleDegrees - currentAngle;

	// calculate the rotation axis
	CartesianPoint pA1Centered = _atom1.getCoor()- _atom2.getCoor();
	CartesianPoint pA3Centered = _atom3.getCoor()- _atom2.getCoor();
	CartesianPoint rotAxis = pA1Centered.cross(pA3Centered);

	// build the list of atoms to be moved
	set<Atom*> excludeList;
	excludeList.insert(&_atom1);
	excludeList.insert(&_atom2);
	set<Atom*> moveList = _atom3.findLinkedAtoms(excludeList);
	moveList.insert(&_atom3);
	set<Atom*> secondMoveList; // used only for natural rotations
	if (naturalMovementOnSetDOF_flag) {
		excludeList.clear();
		excludeList.insert(&_atom2);
		excludeList.insert(&_atom3);
		secondMoveList = _atom1.findLinkedAtoms(excludeList);
		secondMoveList.insert(&_atom1);
	}

	if (!naturalMovementOnSetDOF_flag) {
		// get the rotation matrix
		Matrix m = CartesianGeometry::getRotationMatrix(rotation, rotAxis);

		// move the atoms
		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			rotateAtom(**k, m, _atom2.getCoor()); 
		}
	} else {
		// calculate the geometric centers of the two sets
		CartesianPoint CoM1;
		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			CoM1 += (*k)->getCoor();
		}
		CoM1 /= (double)moveList.size();
		CartesianPoint CoM2;
		for (set<Atom*>::iterator k=secondMoveList.begin(); k!=secondMoveList.end(); k++) {
			CoM2 += (*k)->getCoor();
		}
		CoM2 /= (double)secondMoveList.size();

		// calculate the distances from the axis of rotation * list size
		double f1 = (double)moveList.size() * CartesianGeometry::distanceFromLine(CoM1, _atom2.getCoor(),  _atom3.getCoor());
		double f2 = (double)secondMoveList.size() * CartesianGeometry::distanceFromLine(CoM2, _atom2.getCoor(),  _atom3.getCoor());
		
		// normalize
		double tot = f1 + f2;
		f1 /= -tot;
		f2 /= tot;

		// get the rotation matrices
		Matrix m1 = CartesianGeometry::getRotationMatrix(rotation * f2, rotAxis);
		Matrix m2 = CartesianGeometry::getRotationMatrix(rotation * f1, rotAxis);

		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			rotateAtom(**k, m1, _atom2.getCoor()); 
		}
		for (set<Atom*>::iterator k=secondMoveList.begin(); k!=secondMoveList.end(); k++) {
			rotateAtom(**k, m2, _atom2.getCoor()); 
		}
	}

	return true;
}

bool Transforms::setDihedral(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _angleDegrees, bool _strict) {
	/* rotate _atom4 in the same direction of the _atom2-_atom3 bond (they do not
	   need to be bonded through) to place the dihedral at the right value, and move 
	   any other atom connected to _atom3 (bonded and those bonded that follow) as 
	   well, if the bonded information is stored in the atoms.  If _strict  = true then
	   the other atoms connected to _atom3 won't move, except for _atom4 */

	while (_angleDegrees < -180.0) {
		_angleDegrees += 360.0;
	}
	while (_angleDegrees > 180.0) {
		_angleDegrees -= 360.0;
	}

	// bring the measured angle into the -180 to 180 range and then make it positive
	double currentAngle = _atom1.dihedral(_atom2, _atom3, _atom4);
	while (currentAngle < -180.0) {
		currentAngle += 360.0;
	}
	while (currentAngle > 180.0) {
		currentAngle -= 360.0;
	}
	double rotation = _angleDegrees - currentAngle;

	// calculate the rotation axis
	CartesianPoint rotAxis = _atom3.getCoor()- _atom2.getCoor();

	// build the list of atoms to be moved
	set<Atom*> excludeList;
	set<Atom*> moveList;
	set<Atom*> secondMoveList; // used only for natural rotations
	if (_strict) {
		// move only _atom4 and connected atoms
		excludeList.insert(&_atom1);
		excludeList.insert(&_atom2);
		excludeList.insert(&_atom3);
		moveList = _atom4.findLinkedAtoms(excludeList);
		moveList.insert(&_atom4);
		if (naturalMovementOnSetDOF_flag) {
			excludeList.clear();
			excludeList.insert(&_atom4);
			secondMoveList = _atom3.findLinkedAtoms(excludeList);
		}
	} else {
		// move all atoms connected to _atom3 (excluding _atom2)
		excludeList.insert(&_atom1);
		excludeList.insert(&_atom2);
		moveList = _atom3.findLinkedAtoms(excludeList);
		if (naturalMovementOnSetDOF_flag) {
			excludeList.clear();
			excludeList.insert(&_atom3);
			excludeList.insert(&_atom4);
			secondMoveList = _atom2.findLinkedAtoms(excludeList);
		}
	}

	// move the atoms
	if (!naturalMovementOnSetDOF_flag) {
		// get the rotation matrix
		Matrix m = CartesianGeometry::getRotationMatrix(rotation, rotAxis);

		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			rotateAtom(**k, m, _atom2.getCoor()); 
		}
	} else {
		// calculate the geometric centers of the two sets
		CartesianPoint CoM1;
		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			CoM1 += (*k)->getCoor();
		}
		CoM1 /= (double)moveList.size();
		CartesianPoint CoM2;
		for (set<Atom*>::iterator k=secondMoveList.begin(); k!=secondMoveList.end(); k++) {
			CoM2 += (*k)->getCoor();
		}
		CoM2 /= (double)secondMoveList.size();

		// calculate the distances from the axis of rotation * list size
		double f1 = (double)moveList.size() * CartesianGeometry::distanceFromLine(CoM1, _atom2.getCoor(),  _atom3.getCoor());
		double f2 = (double)secondMoveList.size() * CartesianGeometry::distanceFromLine(CoM2, _atom2.getCoor(),  _atom3.getCoor());
		
		// normalize
		double tot = f1 + f2;
		f1 /= -tot;
		f2 /= tot;

		// get the rotation matrices
		Matrix m1 = CartesianGeometry::getRotationMatrix(rotation * f2, rotAxis);
		Matrix m2 = CartesianGeometry::getRotationMatrix(rotation * f1, rotAxis);

		for (set<Atom*>::iterator k=moveList.begin(); k!=moveList.end(); k++) {
			rotateAtom(**k, m1, _atom2.getCoor()); 
		}
		for (set<Atom*>::iterator k=secondMoveList.begin(); k!=secondMoveList.end(); k++) {
			rotateAtom(**k, m2, _atom2.getCoor()); 
		}
	}


	return true;
}

bool Transforms::setImproper(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _angleDegrees) {
	return setDihedral(_atom1, _atom2, _atom3, _atom4, _angleDegrees, true);
}

