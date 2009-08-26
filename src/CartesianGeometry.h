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

#ifndef CARTESIANGEOMETRY_H
#define CARTESIANGEOMETRY_H

//#define M_PI = 3.14159265358979323846;

#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include "CartesianPoint.h"
#include "Matrix.h"

using namespace std;

//class CartesianPoint;
//class Matrix;

class CartesianGeometry {

	public:
		static CartesianGeometry * instance();
		CartesianPoint addCartesianPoints(const CartesianPoint & _point, const CartesianPoint & _translation); // Move _point along x, y, and z of _translation
		CartesianPoint matrixTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix); // About origin
		CartesianPoint matrixTransposeTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix);



		double distance(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const;
		double distance2(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const;
		vector<double> distanceDerivative(CartesianPoint & _firstCartesianPoint,CartesianPoint & _secondCartesianPoint);
		vector<double> distanceNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint, const double _deltaSize=0.01) const;

		// Are there faster ways to do angle and dihedral?  Any way without needing to take arc's?
		double angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const;
		double angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const;
		double angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const;
		double angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const;
		double cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) const;
		double cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) const;
		vector<double> angleDerivative( CartesianPoint & _firstCartesianPoint,  CartesianPoint & _center,  CartesianPoint & _secondCartesianPoint) ;
		vector<double> angleNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint, const double _deltaSize=0.01) const;

		double dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const;
		double dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const;
		double cosDihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) const;
		vector<double> dihedralNumericalDerivative(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, const double _deltaSize=0.01) const;

		CartesianPoint build(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral); // build an atom from a distance, angle, and a dihedral with the given points
		CartesianPoint buildRadians(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral); // build an atom from a distance, angle, and a dihedral with the given points
		void seed(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle); // Places three atoms at origin, on x axis, and in xy plane, respectively
		void seedRadians(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle); // Places three atoms at origin, on x axis, and in xy plane, respectively

		double radiansToDegrees(double _rad) { return _rad * 180 / M_PI; };

		Matrix getRotationMatrix(double degrees, const CartesianPoint & _axis) const;
		Matrix getXRotationMatrix(double degrees) const;
		Matrix getYRotationMatrix(double degrees) const;
		Matrix getZRotationMatrix(double degrees) const;

		CartesianPoint projection(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2=CartesianPoint(0.0, 0.0, 0.0)) const; // projection of the point on a line
		double distanceFromLine(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2=CartesianPoint(0.0, 0.0, 0.0)) const; // projection of the point on a line
		double distanceFromSegment(const CartesianPoint & _p, const CartesianPoint & _center, CartesianPoint & _axis) const; 

		double planeAngle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _q1, const CartesianPoint & _q2, const CartesianPoint & _q3) const;
		CartesianPoint normalToPlane(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3) const;


	protected:

		CartesianGeometry();
		CartesianGeometry(const CartesianGeometry & theInstance);
		void operator= (const CartesianGeometry & theInstance);

};

#endif
