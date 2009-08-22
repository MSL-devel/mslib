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

#include "CartesianGeometry.h"

CartesianGeometry * CartesianGeometry::instance() {
	static CartesianGeometry inst;
	return &inst;
}

CartesianGeometry::CartesianGeometry() {
}

CartesianGeometry::CartesianGeometry(const CartesianGeometry & theInstance)
{
}

void CartesianGeometry::operator= (const CartesianGeometry & theInstance)
{	
}

CartesianPoint CartesianGeometry::addCartesianPoints(const CartesianPoint & _point, const CartesianPoint & _translation) {
	double new_x = _point.getX() + _translation.getX();
	double new_y = _point.getY() + _translation.getY();
	double new_z = _point.getZ() + _translation.getZ();
	CartesianPoint new_point(new_x, new_y, new_z);
	return new_point;
}

CartesianPoint CartesianGeometry::matrixTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix) {
	double x = _rotationMatrix.getElement(0,0)*_point.getX() + _rotationMatrix.getElement(0,1)*_point.getY() + _rotationMatrix.getElement(0,2)*_point.getZ();
	double y = _rotationMatrix.getElement(1,0)*_point.getX() + _rotationMatrix.getElement(1,1)*_point.getY() + _rotationMatrix.getElement(1,2)*_point.getZ();
	double z = _rotationMatrix.getElement(2,0)*_point.getX() + _rotationMatrix.getElement(2,1)*_point.getY() + _rotationMatrix.getElement(2,2)*_point.getZ();
	return CartesianPoint(x,y,z);
}

CartesianPoint CartesianGeometry::matrixTransposeTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix) {
	double x = _rotationMatrix.getElement(0,0)*_point.getX() + _rotationMatrix.getElement(1,0)*_point.getY() + _rotationMatrix.getElement(2,0)*_point.getZ();
	double y = _rotationMatrix.getElement(0,1)*_point.getX() + _rotationMatrix.getElement(1,1)*_point.getY() + _rotationMatrix.getElement(2,1)*_point.getZ();
	double z = _rotationMatrix.getElement(0,2)*_point.getX() + _rotationMatrix.getElement(1,2)*_point.getY() + _rotationMatrix.getElement(2,2)*_point.getZ();
	return CartesianPoint(x,y,z);
}

double CartesianGeometry::distance(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const
{
	CartesianPoint difference = _firstCartesianPoint - _secondCartesianPoint;
	return sqrt(difference * difference);
}

double CartesianGeometry::distance2(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const
{
	CartesianPoint difference = _firstCartesianPoint - _secondCartesianPoint;
	return difference * difference;
}

double CartesianGeometry::angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const {
	return angleRadians(_firstCartesianPoint, _secondCartesianPoint) * 180.0 / M_PI;
}

double CartesianGeometry::angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const {
	return angleRadians(_firstCartesianPoint, _center, _secondCartesianPoint) * 180.0 / M_PI;
}

double CartesianGeometry::angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const {
	double dotp = _firstCartesianPoint.getUnit() * _secondCartesianPoint.getUnit();
	
	// the following necessary for value very
	// close to 1 but just above
	if (dotp > 1.0) {
		dotp = 1.0;
	}
	else if (dotp < -1.0) {
		dotp = - 1.0;
	}
	return acos(dotp);
}

double CartesianGeometry::cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const {
	double dotp = _firstCartesianPoint.getUnit() * _secondCartesianPoint.getUnit();
	
	// the following necessary for value very
	// close to 1 but just above
	if (dotp > 1.0) {
		dotp = 1.0;
	}
	else if (dotp < -1.0) {
		dotp = - 1.0;
	}
	return dotp;
}

double CartesianGeometry::angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const {
	CartesianPoint diff1 = _firstCartesianPoint - _center;
	CartesianPoint diff2 = _secondCartesianPoint - _center;
	return angleRadians(diff1, diff2);
}
double CartesianGeometry::cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const {
	CartesianPoint diff1 = _firstCartesianPoint - _center;
	CartesianPoint diff2 = _secondCartesianPoint - _center;
	return cosAngle(diff1, diff2);
}
double CartesianGeometry::dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const {
	return dihedralRadians(_p1, _p2, _p3, _p4) * 180.0 / M_PI;
}

double CartesianGeometry::dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const {
	CartesianPoint AB = _p1 - _p2;
	CartesianPoint CB = _p3 - _p2;
	CartesianPoint DC = _p4 - _p3;

	if (AB.length() == 0.0 || CB.length() == 0.0 || DC.length() == 0.0) {
		return 0.;
	}
	else {
		CartesianPoint ABxCB = AB.cross(CB).getUnit();
		CartesianPoint DCxCB = DC.cross(CB).getUnit();

		// the following necessary for value very
		// close to 1 but just above
		double dotp = ABxCB * DCxCB;
		if (dotp > 1.0) {
			dotp = 1.0;
		}
		else if (dotp < -1.0) {
			dotp = -1.0;
		}

		double angle = acos(dotp);
		if (ABxCB * DC > 0) {
			angle *= -1;
		}
		return angle;
	}
}
double CartesianGeometry::cosDihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const {
	CartesianPoint AB = _p1 - _p2;
	CartesianPoint CB = _p3 - _p2;
	CartesianPoint DC = _p4 - _p3;

	if (AB.length() == 0.0 || CB.length() == 0.0 || DC.length() == 0.0) {
		return 0.;
	}

	CartesianPoint ABxCB = AB.cross(CB).getUnit();
	CartesianPoint DCxCB = DC.cross(CB).getUnit();

	// the following necessary for value very
	// close to 1 but just above
	double dotp = ABxCB * DCxCB;
	if (dotp > 1.0) {
		dotp = 1.0;
	}
	else if (dotp < -1.0) {
		dotp = -1.0;
	}

	return dotp;
}

/*
CartesianPoint CartesianGeometry::build(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double _distance, const double _angle, const double _dihedral)
{
	// Need some check that points are there?

	// Move to radians representation

	// radii and unit vector from _distAtom to _angleAtom
	CartesianPoint rab = _distAtom - _angleAtom;
	CartesianPoint uab = rab.getUnit();

	// radii and unit vector from _angleAtom to _dihedralAtom
	CartesianPoint rbc = _angleAtom - _dihedralAtom;
	CartesianPoint ubc = rbc.getUnit();

	double rsin = _dihedral * sin(M_PI - _angle);
	double rcos = _dihedral * cos(M_PI - _angle);
	double rsinsin = rsin * sin(M_PI + _dihedral);
	double rsincos = rsin * cos(M_PI + _dihedral);

	/ ****************************************
	* First set component in direction of the
	* (b->a) bond.
	**************************************** /
	CartesianPoint newCartesianPoint = uab * rcos;

	/ ****************************************
	* Then add component perpendicular to bond,
	* and colinear with opposing dihedral
	* vector (b->c)
	**************************************** /
	CartesianPoint ux = rbc - (uab * (rbc * uab));
	CartesianPoint rx = ux.getUnit() * rsincos;
	newCartesianPoint += rx;

	/ ****************************************
	* Then add component perpendicular to bond,
	* and perpendicular to opposing dihedral
	* vector (b->c)
	**************************************** /
	CartesianPoint uy = uab.cross(ubc);
	CartesianPoint ry = uy.getUnit() * rsinsin;
	newCartesianPoint += ry;

	newCartesianPoint += _distAtom;
	return newCartesianPoint;
}
*/

CartesianPoint CartesianGeometry::build(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral) {
	return buildRadians(_distAtom, _angleAtom, _dihedralAtom, _distance, _angle * M_PI / 180, _dihedral * M_PI / 180);
}

CartesianPoint CartesianGeometry::buildRadians(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral) {
	/***************************************************
	 * This function sets the coordinates of a cartesian
	 * point A:
	 *  - atoms B C D positions 
	 *  - the distance from atom B
	 *  - the angle A-B-C
	 *  - the dihedral A-B-C-D
	 *  
	 *                    A
	 *                     \
	 *                      B--C
	 *                          \
	 *                           D
	 *
	 * Angles are in RADIANS
	 *
	 * No check points are coded here (the distance and the
	 * angle should not be zero) and the atoms should not
	 * be overlapping or B-C-D be a 180 angle
	 *
	 ***************************************************/

	// unit vector from _distAtom to _angleAtom
	CartesianPoint uab = (_distAtom - _angleAtom).getUnit();

	// radii from _angleAtom to _dihedralAtom
	CartesianPoint rbc = _angleAtom - _dihedralAtom;

	double angle2 = M_PI - _angle;
	double dihe2 = M_PI + _dihedral;
	double rsin = _distance * sin(angle2);
	double rcos = _distance * cos(angle2);
	double rsinsin = rsin * sin(dihe2);
	double rsincos = rsin * cos(dihe2);

	/****************************************
	* First set component in direction of the
	* (b->a) bond:
	*  >> (uab * rcos)
	*
	* Then add component perpendicular to bond,
	* and colinear with opposing dihedral
	* vector (b->c)
	*  >>  (rbc - (uab * (rbc * uab))).getUnit() * rsincos
	*
	* Then add component perpendicular to bond,
	* and perpendicular to opposing dihedral
	* vector (b->c)
	*  >> (uab.cross(rbc)).getUnit() * rsinsin
	*
	* Then move it relative to the distance atom
	*  >> _distAtom
	*
	****************************************/
	return (uab * rcos)  +  ((rbc - (uab * (rbc * uab))).getUnit() * rsincos)  +  ((uab.cross(rbc)).getUnit() * rsinsin)  +  _distAtom;

}


void CartesianGeometry::seed(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle) {
	seedRadians(_originCartesianPoint, _distAtom, _angleAtom, _distance12, _distance23, _angle * M_PI / 180);
}

void CartesianGeometry::seedRadians(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle)
{
	if (_distance12 == 0 || _distance23 == 0 || _angle == 0) {
		cerr << "ERROR 5107: cannot seed atoms because of zero values in arguments" <<  endl;
		exit (5107);
	}
	// put _originCartesianPoint at the origin
	_originCartesianPoint.setCoor(0.,0.,0.);
	// put A on the X axis at _distance distance from B
	_distAtom = _distAtom + CartesianPoint(_distance12, 0., 0.);
	// put this on the XY plane at d distance from the origin and 
	// the correct angle from A
	double rsin = _distance23 * sin(M_PI - _angle);
	double rcos = _distance23 * cos(M_PI - _angle);
	// translate *this on the X axis of the coordinates of A
	_angleAtom.setCoor(rcos + _distAtom.getX(), rsin, 0.0);
}

Matrix CartesianGeometry::getRotationMatrix(double degrees, const CartesianPoint & _axis) const {
	double radiants = degrees * M_PI / 180.0;
	double cosRad = cos(radiants);
	double sinRad = sin(radiants);
	double tRad = 1.0 - cosRad;

	// to avoid rounding errors
	if (degrees == 0.0) {
		cosRad = 1.0;
		sinRad = 0.0;
	} else if (degrees == 90.0) {
		cosRad = 0.0;
		sinRad = 1.0;
	} else if (degrees == 180.0) {
		cosRad = -1.0;
		sinRad = 0.0;
	} else if (degrees == 270.0) {
		cosRad = 0.0;
		sinRad = -1.0;
	}

	CartesianPoint n = _axis.getUnit();

	// rotation matrix
	Matrix m(3, 3, 0.0);
	m[0][0] = tRad * n[0] * n[0] + cosRad;
	m[0][1] = tRad * n[0] * n[1] - sinRad * n[2];
	m[0][2] = tRad * n[0] * n[2] + sinRad * n[1];
	m[1][0] = tRad * n[0] * n[1] + sinRad * n[2];
	m[1][1] = tRad * n[1] * n[1] + cosRad;
	m[1][2] = tRad * n[1] * n[2] - sinRad * n[0];
	m[2][0] = tRad * n[0] * n[2] - sinRad * n[1];
	m[2][1] = tRad * n[1] * n[2] + sinRad * n[0];
	m[2][2] = tRad * n[2] * n[2] + cosRad;

	return m;
}

Matrix CartesianGeometry::getXRotationMatrix(double degrees) const {

	double radiants = degrees * M_PI / 180;
	double cosRad = cos(radiants);
	double sinRad = sin(radiants);

	// to avoid rounding errors
	if (degrees == 0.0) {
		cosRad = 1.0;
		sinRad = 0.0;
	} else if (degrees == 90.0) {
		cosRad = 0.0;
		sinRad = 1.0;
	} else if (degrees == 180.0) {
		cosRad = -1.0;
		sinRad = 0.0;
	} else if (degrees == 270.0) {
		cosRad = 0.0;
		sinRad = -1.0;
	}

	// rotation matrix
	Matrix m(3, 3, 0.0);
	m[0][0] = 1;
	m[1][1] = cosRad;
	m[2][2] = cosRad;
	m[1][2] = sinRad;
	m[2][1] = -sinRad;
	return m;

}
Matrix CartesianGeometry::getYRotationMatrix(double degrees) const {

	double radiants = degrees * M_PI / 180;
	double cosRad = cos(radiants);
	double sinRad = sin(radiants);

	// to avoid rounding errors
	if (degrees == 0.0) {
		cosRad = 1.0;
		sinRad = 0.0;
	} else if (degrees == 90.0) {
		cosRad = 0.0;
		sinRad = 1.0;
	} else if (degrees == 180.0) {
		cosRad = -1.0;
		sinRad = 0.0;
	} else if (degrees == 270.0) {
		cosRad = 0.0;
		sinRad = -1.0;
	}

	// rotation matrix
	Matrix m(3, 3, 0.0);
	m[0][0] = cosRad;
	m[1][1] = 1;
	m[2][2] = cosRad;
	m[0][2] = sinRad;
	m[2][0] = -sinRad;
	return m;

}
Matrix CartesianGeometry::getZRotationMatrix(double degrees) const {

	double radiants = degrees * M_PI / 180;
	double cosRad = cos(radiants);
	double sinRad = sin(radiants);

	// to avoid rounding errors
	if (degrees == 0.0) {
		cosRad = 1.0;
		sinRad = 0.0;
	} else if (degrees == 90.0) {
		cosRad = 0.0;
		sinRad = 1.0;
	} else if (degrees == 180.0) {
		cosRad = -1.0;
		sinRad = 0.0;
	} else if (degrees == 270.0) {
		cosRad = 0.0;
		sinRad = -1.0;
	}

	// rotation matrix
	Matrix m(3, 3, 0.0);
	m[0][0] = cosRad;
	m[0][1] = sinRad;
	m[1][0] = -sinRad;
	m[1][1] = cosRad;
	m[2][2] = 1;
	return m;

}

CartesianPoint CartesianGeometry::projection(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2) const {
	// projection of the point on a line
	
	if (_axis1.distance(_axis2) == 0) {
		// ERROR HANDLING HERE!!!
		cerr << "ERROR 49321, arguments _axis1 and _axis2 are identical in CartesianPoint CartesianGeometry::projection(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2) const" << endl;
		exit(49321);
	}

	CartesianPoint direction = _axis2 - _axis1;

	//double t = ((_p - _axis1) * direction) / pow(direction.length(), 2);
	double t = ((_p - _axis1) * direction) / (direction * direction);
	return CartesianPoint(_axis1 + direction * t);
}

double CartesianGeometry::distanceFromLine(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2) const {
	CartesianPoint proj = projection(_p, _axis1, _axis2);
	return _p.distance(proj);
}


double CartesianGeometry::distanceFromSegment(const CartesianPoint & _p, const CartesianPoint & _center, CartesianPoint & _axis) const {
	CartesianPoint p = (*CartesianGeometry::instance()).projection(_p,_center,_axis) - _center;
	double dp = p * _axis;
	if (dp < 0.0) {
		return _p.distance(_center);
	} else if (dp > _axis * _axis) {
		return _p.distance(_center + _axis);
	} else {
		return (*CartesianGeometry::instance()).distanceFromLine(_p, _center, _axis);
	}
		
}

CartesianPoint CartesianGeometry::normalToPlane(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3) const {
	return (_p1 - _p2).cross(_p1 - _p3).getUnit();
}

double CartesianGeometry::planeAngle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _q1, const CartesianPoint & _q2, const CartesianPoint & _q3) const {
	return normalToPlane(_p1, _p2, _p3).angle(normalToPlane(_q1, _q2, _q3));
}
