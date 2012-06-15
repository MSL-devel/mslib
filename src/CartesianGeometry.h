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

#ifndef CARTESIANGEOMETRY_H_
#define CARTESIANGEOMETRY_H_

#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include "CartesianPoint.h"
#include "Matrix.h"

namespace MSL {
     namespace CartesianGeometry {

		CartesianPoint addCartesianPoints(const CartesianPoint & _point, const CartesianPoint & _translation); // Move _point along x, y, and z of _translation
		CartesianPoint matrixTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix); // About origin
		CartesianPoint matrixTransposeTimesCartesianPoint(const CartesianPoint & _point, const Matrix & _rotationMatrix);

		double distance(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint);
		double distance2(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint);
		double distanceDerivative(CartesianPoint & _firstCartesianPoint, CartesianPoint & _secondCartesianPoint, std::vector<double>* grad = NULL); // return distance
		double distanceNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint, std::vector<double>* grad = NULL, const double _deltaSize=0.01);// return distance


		double angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint);
		double angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint);
		double angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint);
		double angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint);
		double cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint);
		double cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint);
		double angleDerivative( CartesianPoint & _firstCartesianPoint,  CartesianPoint & _center,  CartesianPoint & _secondCartesianPoint, std::vector<double>* grad = NULL) ; // return angle in radians
		double angleNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint, std::vector<double>* grad = NULL, const double _deltaSize=0.01); // return angle in radians

		double dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4);
		double dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4);
		double cosDihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4);
		double dihedralNumericalCosDerivative(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, std::vector<double>* grad = NULL, const double _deltaSize=0.01);
		double dihedralCosDerivative(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4, std::vector<double>* grad = NULL);
		double dihedralDerivative(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4, std::vector<double> *grad = NULL);

		CartesianPoint build(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral); // build an atom from a distance, angle, and a dihedral with the given points
		CartesianPoint buildRadians(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral); // build an atom from a distance, angle, and a dihedral with the given points
		void seed(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle); // Places three atoms at origin, on x axis, and in xy plane, respectively
		void seedRadians(CartesianPoint & _originCartesianPoint, CartesianPoint & _distAtom, CartesianPoint & _angleAtom, const double _distance12, const double _distance23, const double _angle); // Places three atoms at origin, on x axis, and in xy plane, respectively

		double radiansToDegrees(double _rad);

		Matrix getRotationMatrix(double degrees, const CartesianPoint & _axis);
		Matrix getXRotationMatrix(double degrees);
		Matrix getYRotationMatrix(double degrees);
		Matrix getZRotationMatrix(double degrees);

		CartesianPoint projection(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2=CartesianPoint(0.0, 0.0, 0.0)); // projection of the point on a line
		double distanceFromLine(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2=CartesianPoint(0.0, 0.0, 0.0)); // projection of the point on a line
		double distanceFromSegment(const CartesianPoint & _p, const CartesianPoint & _center, CartesianPoint & _axis); 

		double planeAngle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _q1, const CartesianPoint & _q2, const CartesianPoint & _q3);
		CartesianPoint normalToPlane(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3);
     }
}

#endif
