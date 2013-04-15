/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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

#include "CartesianGeometry.h"

using namespace MSL;
using namespace std;

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

double CartesianGeometry::radiansToDegrees(double _rad) { return _rad * 180 / M_PI; };

double CartesianGeometry::distance(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint)
{
	CartesianPoint difference = _firstCartesianPoint - _secondCartesianPoint;
	return sqrt(difference * difference);
}

double CartesianGeometry::distance2(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint)
{
	CartesianPoint difference = _firstCartesianPoint - _secondCartesianPoint;
	return difference * difference;
}

double CartesianGeometry::distanceNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint, vector<double>* _partialDerivatives, const double _deltaSize)
{
	CartesianPoint p1 = _firstCartesianPoint;
	CartesianPoint p2 = _secondCartesianPoint;
	
	(*_partialDerivatives).resize(6,0);
	// change p1.x +/- 0.01
	p1[0] += _deltaSize;
	double dx1a = distance(p1,p2);

	p1[0] -= _deltaSize*2;
	double dx1b = distance(p1,p2);

	p1[0] += _deltaSize;

	double dx1 =  (dx1a - dx1b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[0] =  dx1;
	
	// change p1.y +/- 0.01
	p1[1] += _deltaSize;
	double dy1a = distance(p1,p2);

	p1[1] -= _deltaSize*2;
	double dy1b = distance(p1,p2);

	p1[1] += _deltaSize;

	double dy1 =  (dy1a - dy1b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[1] =  dy1;

	// change p1.z +/- 0.01
	p1[2] += _deltaSize;
	double dz1a = distance(p1,p2);

	p1[2] -= _deltaSize*2;
	double dz1b = distance(p1,p2);

	p1[2] += _deltaSize;

	double dz1 =  (dz1a - dz1b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[2] = dz1;


	// change p2.x +/- 0.01
	p2[0] += _deltaSize;
	double dx2a = distance(p1,p2);

	p2[0] -= _deltaSize*2;
	double dx2b = distance(p1,p2);

	p2[0] += _deltaSize;

	double dx2 =  (dx2a - dx2b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[3] = dx2;


	// change p2.y +/- 0.01
	p2[1] += _deltaSize;
	double dy2a = distance(p1,p2);

	p2[1] -= _deltaSize*2;
	double dy2b = distance(p1,p2);

	p2[1] += _deltaSize;

	double dy2 =  (dy2a - dy2b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[4] = dy2;

	// change p2.z +/- 0.01
	p2[2] += _deltaSize;
	double dz2a = distance(p1,p2);

	p2[2] -= _deltaSize*2;
	double dz2b = distance(p1,p2);

	p2[2] += _deltaSize;

	double dz2 =  (dz2a - dz2b) / (2 *_deltaSize) ;
	(*_partialDerivatives)[5] = dz2;

	return distance(p1,p2);

}

double CartesianGeometry::distanceDerivative(CartesianPoint & _firstCartesianPoint, CartesianPoint & _secondCartesianPoint, vector<double>* grad){

	(*grad).resize(6,0);
	CartesianPoint difference = _firstCartesianPoint - _secondCartesianPoint;
	double dist = sqrt(difference * difference);

	if (grad != NULL) {
	        double EPS = 0.0000000000000001; // 10 ^ -15
        	if ( dist < EPS ) {

	                (*grad)[0] = (_firstCartesianPoint[0] < _secondCartesianPoint[0] ? -1 : 1);
        	        (*grad)[1] = (_firstCartesianPoint[1] < _secondCartesianPoint[1] ? -1 : 1);
                	(*grad)[2] = (_firstCartesianPoint[2] < _secondCartesianPoint[2] ? -1 : 1);

	                (*grad)[3] = (_secondCartesianPoint[0] < _firstCartesianPoint[0] ? -1 : 1);
        	        (*grad)[4] = (_secondCartesianPoint[1] < _firstCartesianPoint[1] ? -1 : 1);
                	(*grad)[5] = (_secondCartesianPoint[2] < _firstCartesianPoint[2] ? -1 : 1);


	        }  else {

                	(*grad)[0] = (_firstCartesianPoint[0] - _secondCartesianPoint[0]) / dist;
        	        (*grad)[1] = (_firstCartesianPoint[1] - _secondCartesianPoint[1]) / dist;
	                (*grad)[2] = (_firstCartesianPoint[2] - _secondCartesianPoint[2]) / dist;

        	        (*grad)[3] = - (*grad)[0];
                	(*grad)[4] = - (*grad)[1];
	                (*grad)[5] = - (*grad)[2];

	        }
	}
	return dist;
}


double CartesianGeometry::angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) {
	return angleRadians(_firstCartesianPoint, _secondCartesianPoint) * 180.0 / M_PI;
}

double CartesianGeometry::angle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) {
	return angleRadians(_firstCartesianPoint, _center, _secondCartesianPoint) * 180.0 / M_PI;
}

double CartesianGeometry::angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) {
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

double CartesianGeometry::cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _secondCartesianPoint) {
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

double CartesianGeometry::angleRadians(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) {
	CartesianPoint diff1 = _firstCartesianPoint - _center;
	CartesianPoint diff2 = _secondCartesianPoint - _center;
	return angleRadians(diff1, diff2);
}
double CartesianGeometry::cosAngle(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint) {
	CartesianPoint diff1 = _firstCartesianPoint - _center;
	CartesianPoint diff2 = _secondCartesianPoint - _center;
	return cosAngle(diff1, diff2);
}

double CartesianGeometry::angleNumericalDerivative(const CartesianPoint & _firstCartesianPoint, const CartesianPoint & _center, const CartesianPoint & _secondCartesianPoint, vector<double>* _partialDerivatives, const double _deltaSize) 
{
	(*_partialDerivatives).resize(9,0);
	CartesianPoint p1 = _firstCartesianPoint;
	CartesianPoint pC = _center;
	CartesianPoint p2 = _secondCartesianPoint;

	// change p1.x +/- 0.01
	p1[0] += _deltaSize;
	double Ax1a = angleRadians(p1,pC,p2);
		
	p1[0] -= 2 * _deltaSize;

	double Ax1b = angleRadians(p1,pC,p2);
	
	p1[0] += _deltaSize;

	double Ax1 = (Ax1a - Ax1b) / (2 * _deltaSize);
	(*_partialDerivatives)[0] = Ax1;

	// change p1.y +/- 0.01
	p1[1] += _deltaSize;
	double Ay1a = angleRadians(p1,pC,p2);
		
	p1[1] -= 2 * _deltaSize;

	double Ay1b = angleRadians(p1,pC,p2);
	
	p1[1] += _deltaSize;

	double Ay1 = (Ay1a - Ay1b) / (2 * _deltaSize);
	(*_partialDerivatives)[1] = Ay1;



	// change p1.z +/- 0.01
	p1[2] += _deltaSize;
	double Az1a = angleRadians(p1,pC,p2);
		
	p1[2] -= 2 * _deltaSize;

	double Az1b = angleRadians(p1,pC,p2);
	
	p1[2] += _deltaSize;

	double Az1 = (Az1a - Az1b) / (2 * _deltaSize);
	(*_partialDerivatives)[2] = Az1;

	// change pC.x +/- 0.01
	pC[0] += _deltaSize;
	double Ax3a = angleRadians(p1,pC,p2);
		
	pC[0] -= 2 * _deltaSize;

	double Ax3b = angleRadians(p1,pC,p2);
	
	pC[0] += _deltaSize;

	double Ax3 = (Ax3a - Ax3b) / (2 * _deltaSize);
	(*_partialDerivatives)[3] = Ax3;

	// change pC.y +/- 0.01
	pC[1] += _deltaSize;
	double Ay3a = angleRadians(p1,pC,p2);
		
	pC[1] -= 2 * _deltaSize;

	double Ay3b = angleRadians(p1,pC,p2);
	
	pC[1] += _deltaSize;

	double Ay3 = (Ay3a - Ay3b) / (2 * _deltaSize);
	(*_partialDerivatives)[4] = Ay3;



	// change pC.z +/- 0.01
	pC[2] += _deltaSize;
	double Az3a = angleRadians(p1,pC,p2);
		
	pC[2] -= 2 * _deltaSize;

	double Az3b = angleRadians(p1,pC,p2);
	
	pC[2] += _deltaSize;

	double Az3 = (Az3a - Az3b) / (2 * _deltaSize);
	(*_partialDerivatives)[5] = Az3;


	// change p2.x +/- 0.01
	p2[0] += _deltaSize;
	double Ax2a = angleRadians(p1,pC,p2);
		
	p2[0] -= 2 * _deltaSize;

	double Ax2b = angleRadians(p1,pC,p2);
	
	p2[0] += _deltaSize;

	double Ax2 = (Ax2a - Ax2b) / (2 * _deltaSize);
	(*_partialDerivatives)[6] = Ax2;

	// change p2.y +/- 0.01
	p2[1] += _deltaSize;
	double Ay2a = angleRadians(p1,pC,p2);
		
	p2[1] -= 2 * _deltaSize;

	double Ay2b = angleRadians(p1,pC,p2);
	
	p2[1] += _deltaSize;

	double Ay2 = (Ay2a - Ay2b) / (2 * _deltaSize);
	(*_partialDerivatives)[7] = Ay2;



	// change p2.z +/- 0.01
	p2[2] += _deltaSize;
	double Az2a = angleRadians(p1,pC,p2);
		
	p2[2] -= 2 * _deltaSize;

	double Az2b = angleRadians(p1,pC,p2);
	
	p2[2] += _deltaSize;

	double Az2 = (Az2a - Az2b) / (2 * _deltaSize);
	(*_partialDerivatives)[8] = Az2;

	return angleRadians(p1,pC,p2);
}

double CartesianGeometry::angleDerivative( CartesianPoint & _p1,  CartesianPoint & _p2,  CartesianPoint & _p3, vector<double>* _partialDerivatives ) {
	
	(*_partialDerivatives).resize(9,0);
	CartesianPoint r1 = _p1 - _p2;
	CartesianPoint r2 = _p3 - _p2;

	double L1 = r1.length();
	double L2 = r2.length();

	double p  = r1 * r2;
	
	double EPS = 0.0000000000000001;// 10 ^ -15

	if (fabs(fabs(p) - L1*L2) < EPS) {
		if (sqrt(L1*L2) < EPS) {
			cerr << "POINTS ARE ON TOP OF EACH OTHER, ill defined angular derivative."<<endl;
			(*_partialDerivatives)[0] = (0.0);
			(*_partialDerivatives)[1] = (0.0);
			(*_partialDerivatives)[2] = (0.0);

			(*_partialDerivatives)[3] =(0.0);
			(*_partialDerivatives)[4] =(0.0);
			(*_partialDerivatives)[5] =(0.0);

			(*_partialDerivatives)[6] =(0.0);
			(*_partialDerivatives)[7] =(0.0);
			(*_partialDerivatives)[8] =(0.0);

		} else {

			cerr << "SPECIAL CASE HIT\n";
			(*_partialDerivatives)[0] = ( (1/L2)*sin(acos(r2[0] / L2)) ) ;
			(*_partialDerivatives)[0] *=  (r2-r1)*(CartesianPoint(1,0,0)) >= 0 ? 1 : -1;

			(*_partialDerivatives)[1] = (1/L2)*sin(acos(r2[1] / L2)) ;
			(*_partialDerivatives)[1] *=  (r2-r1)*(CartesianPoint(0,1,0)) >= 0 ? 1: -1;

			(*_partialDerivatives)[2] = ( (1/L2)*sin(acos(r2[2] / L2)) );
			(*_partialDerivatives)[2] *=  (r2-r1)*(CartesianPoint(0,0,1)) >= 0 ? 1: -1;

			/*
			(*_partialDerivatives)[3] = (0.0);
			(*_partialDerivatives)[4] = (0.0);
			(*_partialDerivatives)[5] = (0.0);
			*/

			(*_partialDerivatives)[6] = ( (1/L1)*sin(acos(r1[0] / L1)) );
			(*_partialDerivatives)[6] *=  (r1-r2)*(CartesianPoint(1,0,0)) > 0 ? 1 : -1;

			(*_partialDerivatives)[7] = ( (1/L1)*sin(acos(r1[1] / L1)) );
			(*_partialDerivatives)[7] *=  (r1-r2)*(CartesianPoint(0,1,0)) > 0 ? 1 : -1;

			(*_partialDerivatives)[8] = ( (1/L1)*sin(acos(r1[2] / L1)) );
			(*_partialDerivatives)[8] *=  (r1-r2)*(CartesianPoint(0,0,1)) > 0 ? 1 : -1;

		
			(*_partialDerivatives)[3] = -((*_partialDerivatives)[0] + (*_partialDerivatives)[6]);
			(*_partialDerivatives)[4] = -((*_partialDerivatives)[1] + (*_partialDerivatives)[7]);
			(*_partialDerivatives)[5] = -((*_partialDerivatives)[2] + (*_partialDerivatives)[8]);
		}


	} else {
		double d = p / (L1 * L2);
		double c  = (1 / (sqrt(1 - d*d) * L1*L1*L2*L2));

		(*_partialDerivatives)[0] = ( -c * ( r2[0]*L1*L2 - p*L2*r1[0]/L1 ));
		(*_partialDerivatives)[1] = ( -c * ( r2[1]*L1*L2 - p*L2*r1[1]/L1 ));
		(*_partialDerivatives)[2] = ( -c * ( r2[2]*L1*L2 - p*L2*r1[2]/L1 ));
		/*
		(*_partialDerivatives)[3] = (0.0);
		(*_partialDerivatives)[4] = (0.0);
		(*_partialDerivatives)[5] = (0.0);
		*/

		(*_partialDerivatives)[6] = ( -c * ( r1[0]*L1*L2 - p*L1*r2[0]/L2 ));
		(*_partialDerivatives)[7] = ( -c * ( r1[1]*L1*L2 - p*L1*r2[1]/L2 ));
		(*_partialDerivatives)[8] = ( -c * ( r1[2]*L1*L2 - p*L1*r2[2]/L2 ));

		(*_partialDerivatives)[3] = -((*_partialDerivatives)[0] + (*_partialDerivatives)[6]);
		(*_partialDerivatives)[4] = -((*_partialDerivatives)[1] + (*_partialDerivatives)[7]);
		(*_partialDerivatives)[5] = -((*_partialDerivatives)[2] + (*_partialDerivatives)[8]);
		
	}
	return angleRadians(_p1,_p2,_p3);

}

double CartesianGeometry::dihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) {
	return dihedralRadians(_p1, _p2, _p3, _p4) * 180.0 / M_PI;
}

double CartesianGeometry::dihedralRadians(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) {
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
double CartesianGeometry::cosDihedral(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4) {
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
double CartesianGeometry::dihedralNumericalCosDerivative(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _p4, vector<double>* _partialDerivatives, const double _deltaSize)
{
	
	(*_partialDerivatives).resize(12,0);
	CartesianPoint p1 = _p1;
	CartesianPoint p2 = _p2;
	CartesianPoint p3 = _p3;
	CartesianPoint p4 = _p4;
	

	// change p1.x +/- 0.01
	p1[0] += _deltaSize;
	double Ax1a = cosDihedral(p1,p2,p3,p4);
		
	p1[0] -= 2 * _deltaSize;

	double Ax1b = cosDihedral(p1,p2,p3,p4);
	
	p1[0] += _deltaSize;

	double Ax1 = (Ax1a - Ax1b) / (2 * _deltaSize);
	(*_partialDerivatives)[0] = Ax1;

	// change p1.y +/- 0.01
	p1[1] += _deltaSize;
	double Ay1a = cosDihedral(p1,p2,p3,p4);
		
	p1[1] -= 2 * _deltaSize;

	double Ay1b = cosDihedral(p1,p2,p3,p4);
	
	p1[1] += _deltaSize;

	double Ay1 = (Ay1a - Ay1b) / (2 * _deltaSize);
	(*_partialDerivatives)[1] = Ay1;



	// change p1.z +/- 0.01
	p1[2] += _deltaSize;
	double Az1a = cosDihedral(p1,p2,p3,p4);
		
	p1[2] -= 2 * _deltaSize;

	double Az1b = cosDihedral(p1,p2,p3,p4);
	
	p1[2] += _deltaSize;

	double Az1 = (Az1a - Az1b) / (2 * _deltaSize);
	(*_partialDerivatives)[2] = Az1;



	// change p2.x +/- 0.01
	p2[0] += _deltaSize;
	double Ax2a = cosDihedral(p1,p2,p3,p4);
		
	p2[0] -= 2 * _deltaSize;

	double Ax2b = cosDihedral(p1,p2,p3,p4);
	
	p2[0] += _deltaSize;

	double Ax2 = (Ax2a - Ax2b) / (2 * _deltaSize);
	(*_partialDerivatives)[3] = Ax2;

	// change p2.y +/- 0.01
	p2[1] += _deltaSize;
	double Ay2a = cosDihedral(p1,p2,p3,p4);
		
	p2[1] -= 2 * _deltaSize;

	double Ay2b = cosDihedral(p1,p2,p3,p4);
	
	p2[1] += _deltaSize;

	double Ay2 = (Ay2a - Ay2b) / (2 * _deltaSize);
	(*_partialDerivatives)[4] = Ay2;



	// change p2.z +/- 0.01
	p2[2] += _deltaSize;
	double Az2a = cosDihedral(p1,p2,p3,p4);
		
	p2[2] -= 2 * _deltaSize;

	double Az2b = cosDihedral(p1,p2,p3,p4);
	
	p2[2] += _deltaSize;

	double Az2 = (Az2a - Az2b) / (2 * _deltaSize);
	(*_partialDerivatives)[5] = Az2;



	// change p3.x +/- 0.01
	p3[0] += _deltaSize;
	double Ax3a = cosDihedral(p1,p2,p3,p4);
		
	p3[0] -= 2 * _deltaSize;

	double Ax3b = cosDihedral(p1,p2,p3,p4);
	
	p3[0] += _deltaSize;

	double Ax3 = (Ax3a - Ax3b) / (2 * _deltaSize);
	(*_partialDerivatives)[6] = Ax3;

	// change p3.y +/- 0.01
	p3[1] += _deltaSize;
	double Ay3a = cosDihedral(p1,p2,p3,p4);
		
	p3[1] -= 2 * _deltaSize;

	double Ay3b = cosDihedral(p1,p2,p3,p4);
	
	p3[1] += _deltaSize;

	double Ay3 = (Ay3a - Ay3b) / (2 * _deltaSize);
	(*_partialDerivatives)[7] = Ay3;



	// change p3.z +/- 0.01
	p3[2] += _deltaSize;
	double Az3a = cosDihedral(p1,p2,p3,p4);
		
	p3[2] -= 2 * _deltaSize;

	double Az3b = cosDihedral(p1,p2,p3,p4);
	
	p3[2] += _deltaSize;

	double Az3 = (Az3a - Az3b) / (2 * _deltaSize);
	(*_partialDerivatives)[8] = Az3;



	// change p4.x +/- 0.01
	p4[0] += _deltaSize;
	double Ax4a = cosDihedral(p1,p2,p3,p4);
		
	p4[0] -= 2 * _deltaSize;

	double Ax4b = cosDihedral(p1,p2,p3,p4);
	
	p4[0] += _deltaSize;

	double Ax4 = (Ax4a - Ax4b) / (2 * _deltaSize);
	(*_partialDerivatives)[9] = Ax4;


	// change p4.y +/- 0.01
	p4[1] += _deltaSize;
	double Ay4a = cosDihedral(p1,p2,p3,p4);
		
	p4[1] -= 2 * _deltaSize;

	double Ay4b = cosDihedral(p1,p2,p3,p4);
	
	p4[1] += _deltaSize;

	double Ay4 = (Ay4a - Ay4b) / (2 * _deltaSize);
	(*_partialDerivatives)[10] = Ay4;


	// change p4.z +/- 0.01
	p4[2] += _deltaSize;
	double Az4a = cosDihedral(p1,p2,p3,p4);
		
	p4[2] -= 2 * _deltaSize;

	double Az4b = cosDihedral(p1,p2,p3,p4);
	
	p4[2] += _deltaSize;

	double Az4 = (Az4a - Az4b) / (2 * _deltaSize);
	(*_partialDerivatives)[11] = Az4;

	return cosDihedral(p1,p2,p3,p4);

}

double CartesianGeometry::dihedralCosDerivative(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4, vector<double>* _partialDerivatives){

	(*_partialDerivatives).resize(12,0);
	CartesianPoint r12 = _p1 - _p2;
	CartesianPoint r13 = _p1 - _p3;
	CartesianPoint r21 = _p2 - _p1;
	CartesianPoint r23 = _p2 - _p3;
	CartesianPoint r24 = _p2 - _p4;
	CartesianPoint r31 = _p3 - _p1;
	CartesianPoint r32 = _p3 - _p2;
	CartesianPoint r34 = _p3 - _p4;
	CartesianPoint r42 = _p4 - _p2;
	CartesianPoint r43 = _p4 - _p3;

	double a1 = r12[1]*r32[2] - r12[2]*r32[1];
	double a2 = r43[1]*r32[2] - r43[2]*r32[1];
	double b1 = r12[2]*r32[0] - r12[0]*r32[2];
	double b2 = r43[2]*r32[0] - r43[0]*r32[2];
	double c1 = r12[0]*r32[1] - r12[1]*r32[0];
	double c2 = r43[0]*r32[1] - r43[1]*r32[0];
	

	double L1 = sqrt(a1*a1 + b1*b1 + c1*c1);
	double L2 = sqrt(a2*a2 + b2*b2 + c2*c2);

	double p = a1*a2 + b1*b2 + c1*c2;

	double L1L2 = L1*L2;


	double EPS = 0.0000000000000001;// 10 ^ -15

	if (sqrt(L1*L2) < EPS) { 
		cout << "SPECIAL CASE"<<endl;
		
	} else {

		double LL2 = (1 / L1L2) * (1 / L1L2);
		(*_partialDerivatives)[0] = (   LL2 * ( ( r23[2]*b2 + r32[1]*c2 )*L1L2 - p*(L2/L1)*(r23[2]*b1 + r32[1]*c1)) );
		(*_partialDerivatives)[1] = (   LL2 * ( ( r32[2]*a2 + r23[0]*c2 )*L1L2 - p*(L2/L1)*(r32[2]*a1 + r23[0]*c1)) );
		(*_partialDerivatives)[2] = (   LL2 * ( ( r23[1]*a2 + r32[0]*b2 )*L1L2 - p*(L2/L1)*(r23[1]*a1 + r32[0]*b1)) );

		(*_partialDerivatives)[3] = (   LL2 * ( ( r31[2]*b2 + r34[2]*b1 + r13[1]*c2 + r43[1]*c1)*L1L2 - p*(  (L2/L1) * (r31[2]*b1 + r13[1]*c1) + (L1/L2) * (r34[2]*b2 + r43[1]*c2) ) ) );
		(*_partialDerivatives)[4] = (   LL2 * ( ( r13[2]*a2 + r43[2]*a1 + r31[0]*c2 + r34[0]*c1)*L1L2 - p*(  (L2/L1) * (r13[2]*a1 + r31[0]*c1) + (L1/L2) * (r43[2]*a2 + r34[0]*c2) ) ) );
		(*_partialDerivatives)[5] = (   LL2 * ( ( r31[1]*a2 + r34[1]*a1 + r13[0]*b2 + r43[0]*b1)*L1L2 - p*(  (L2/L1) * (r31[1]*a1 + r13[0]*b1) + (L1/L2) * (r34[1]*a2 + r43[0]*b2) ) ) );

	
		(*_partialDerivatives)[6] = (   LL2 * ( ( r12[2]*b2 + r42[2]*b1 + r21[1]*c2 + r24[1]*c1)*L1L2 - p*(  (L2/L1) * (r12[2]*b1 + r21[1]*c1) + (L1/ L2) * (r42[2]*b2 + r24[1]*c2) ) ) );
		(*_partialDerivatives)[7] = (   LL2 * ( ( r21[2]*a2 + r24[2]*a1 + r12[0]*c2 + r42[0]*c1)*L1L2 - p*(  (L2/L1) * (r21[2]*a1 + r12[0]*c1) + (L1/ L2) * (r24[2]*a2 + r42[0]*c2) ) ) );
		(*_partialDerivatives)[8] = (   LL2 * ( ( r12[1]*a2 + r42[1]*a1 + r21[0]*b2 + r24[0]*b1)*L1L2 - p*(  (L2/L1) * (r12[1]*a1 + r21[0]*b1) + (L1/ L2) * (r42[1]*a2 + r24[0]*b2) ) ) );


		(*_partialDerivatives)[9] = (   LL2 * ( ( r23[2]*b1 + r32[1]*c1 )*L1L2 - p*(L1/L2)*(r23[2]*b2 + r32[1]*c2)) );
		(*_partialDerivatives)[10] = (   LL2 * ( ( r32[2]*a1 + r23[0]*c1 )*L1L2 - p*(L1/L2)*(r32[2]*a2 + r23[0]*c2)) );
		(*_partialDerivatives)[11] = (   LL2 * ( ( r23[1]*a1 + r32[0]*b1 )*L1L2 - p*(L1/L2)*(r23[1]*a2 + r32[0]*b2)) );
	}
   
	return cosDihedral(_p1,_p2,_p3,_p4);
}

double CartesianGeometry::dihedralDerivative(CartesianPoint & _p1, CartesianPoint & _p2, CartesianPoint & _p3, CartesianPoint & _p4, vector<double>* grad) {
	(*grad).resize(12,0);
        double chi = dihedralRadians(_p1, _p2, _p3, _p4);
	double f = -sin(chi);

	double EPS = 0.0000000000000001;// 10 ^ -15
	if (fabs(f) < EPS) {
		// at exactly zero and pi the gradient is difficult to compute exactly, so assume some reasonable value
		for (uint i = 0; i < 12; i++) {
			(*grad)[i] = 1.0;
		}
	} else {
		dihedralCosDerivative(_p1, _p2, _p3, _p4,grad);
		double fi = 1/f;
		for (int i = 0; i < (*grad).size(); i++) {
			(*grad)[i] *= fi;
		}
	}
	return chi;
}

CartesianPoint CartesianGeometry::build(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral) {
	return buildRadians(_distAtom, _angleAtom, _dihedralAtom, _distance, _angle * M_PI / 180, _dihedral * M_PI / 180);
}

CartesianPoint CartesianGeometry::buildRadians(const CartesianPoint & _distAtom, const CartesianPoint & _angleAtom, const CartesianPoint & _dihedralAtom, const double & _distance, const double & _angle, const double & _dihedral) {
	/***************************************************
	 * This function sets the coordinates of a cartesian
	 * point A based on:
	 *  - the position of atoms B C D 
	 *  - the distance of A from B
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
	 * Arguments:
	 *   returned CartesianPoint = atom A 
	 *   _distAtom               = atom B 
	 *   _angleAtom              = atom C 
	 *   _dihedralAtom           = atom D 
	 *   _distance               = A-B distance
	 *   _angle                  = A-B-C angle
	 *   _dihedral               = A-B-C-D dihedral
	 *
	 * NOTE: no check points are coded but the distance and the
	 * angle should not be zero and the atoms should not
	 * be overlapping or B-C-D be a 180 angle
	 *
	 ***************************************************/

	// unit vector from _distAtom to _angleAtom (B - C)
	CartesianPoint uCB = (_distAtom - _angleAtom).getUnit();

	// distance from _angleAtom to _dihedralAtom (C - D)
	CartesianPoint dDC = _angleAtom - _dihedralAtom;

	double angle2 = M_PI - _angle;
	double dihe2 = M_PI + _dihedral;
	double rsin = _distance * sin(angle2);
	double rcos = _distance * cos(angle2);
	double rsinsin = rsin * sin(dihe2);
	double rsincos = rsin * cos(dihe2);

	/****************************************
	* The following creates the resulting position by
	* adding three orthogonal components:
	*
	* - Set a component in the B-C direction ...
	*        (uCB * rcos) + ...
	*
	* - ... add a component on the B-C-D plane (note * denotes dot product used between two vectors) ...
	*        ... + (dDC - (uCB * (dDC * uCB))).getUnit() * rsincos + ...
	*
	* - ... add a component orgogonal to the B-C-D plane
	*        ... + (uCB.cross(dDC)).getUnit() * rsinsin + ...
	*
	* - ... finally, translate the point by the position of atom B
	*        ... + _distAtom
	*
	****************************************/
	return (uCB * rcos)  +  ((dDC - (uCB * (dDC * uCB))).getUnit() * rsincos)  +  ((uCB.cross(dDC)).getUnit() * rsinsin)  +  _distAtom;

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

Matrix CartesianGeometry::getRotationMatrix(double degrees, const CartesianPoint & _axis) {
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

Matrix CartesianGeometry::getXRotationMatrix(double degrees) {

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
	m[1][2] = -sinRad;
	m[2][1] = sinRad;
	return m;

}
Matrix CartesianGeometry::getYRotationMatrix(double degrees) {

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
Matrix CartesianGeometry::getZRotationMatrix(double degrees) {

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
	m[0][1] = -sinRad;
	m[1][0] = sinRad;
	m[1][1] = cosRad;
	m[2][2] = 1;
	return m;

}

CartesianPoint CartesianGeometry::projection(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
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

double CartesianGeometry::distanceFromLine(const CartesianPoint & _p, const CartesianPoint & _axis1, const CartesianPoint & _axis2) {
	CartesianPoint proj = projection(_p, _axis1, _axis2);
	return _p.distance(proj);
}


double CartesianGeometry::distanceFromSegment(const CartesianPoint & _p, const CartesianPoint & _center, CartesianPoint & _axis) {
	CartesianPoint p = projection(_p,_center,_axis) - _center;
	double dp = p * _axis;
	if (dp < 0.0) {
		return _p.distance(_center);
	} else if (dp > _axis * _axis) {
		return _p.distance(_center + _axis);
	} else {
		return distanceFromLine(_p, _center, _axis);
	}
		
}

CartesianPoint CartesianGeometry::normalToPlane(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3) {
	return (_p1 - _p2).cross(_p1 - _p3).getUnit();
}

double CartesianGeometry::planeAngle(const CartesianPoint & _p1, const CartesianPoint & _p2, const CartesianPoint & _p3, const CartesianPoint & _q1, const CartesianPoint & _q2, const CartesianPoint & _q3) {
	return normalToPlane(_p1, _p2, _p3).angle(normalToPlane(_q1, _q2, _q3));
}
