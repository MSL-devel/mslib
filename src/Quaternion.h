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

#ifndef _QUATERNION_H
#define _QUATERNION_H

/*
  Author: dwkulp
  Date  : 2/22/2007
  Note  : Code developed from variety of sources
          including asenes's code, csoto's code, koders.com, 'Game Physics' by Eberly and SireMol, 
 */

// STL Includes
#include <iostream>
#include <vector>
#include <math.h>

// MSL Includes
#include "Matrix.h"
#include "CartesianPoint.h"
#include "AtomPointerVector.h"

// GSL Includes
#ifdef __GSL__
  #include <gsl/gsl_eigen.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
#endif



namespace MSL { 
class Quaternion {  

	public:
		Quaternion();
		Quaternion(double _a, double _v0, double _v1, double _v2);
		Quaternion(double _a, std::vector<double> _v);
 
		~Quaternion();

		void operator+=(const Quaternion & _q);
		void operator-=(const Quaternion & _q);
		void operator*=(const Quaternion & _q);
		void operator*=(double _b);
		void operator/=(double _b);
		Quaternion operator+(const Quaternion & _q);
		Quaternion operator-(const Quaternion & _q);
		Quaternion operator*(const Quaternion & _q);
		Quaternion operator*(double _b);
		Quaternion operator/(double _b);
	        friend std::ostream & operator<<(std::ostream &_os, Quaternion & _q) {_os << _q.toString(); return _os;};


		bool convertToRotationMatrix(Matrix &matrix);                             // Converts this quaternion to  a rotation matrix
		bool convertToQuaternion(Matrix &matrix);                                 // Makes a quaternion (this object) from a rotation matrix
									       	      									 
		bool makeQuaternion(CartesianPoint &axis, const double theta);            // Makes a quaternion (this object) from an axis and an angle
#ifdef __GSL__
		bool makeQuaternion(AtomPointerVector &_align,  AtomPointerVector &_ref);
#endif

		bool rotatePoint(CartesianPoint &pt);                                     // Rotate point by this quaternion store results in pt
		bool rotatePoint(CartesianPoint &pt,CartesianPoint &pt_rot);              // Rotate point by this quaternion store results in pt_rot

		void setInverse();
		Quaternion getInverse() const;

		double norm() const;
		
		Quaternion getConiugate() const;
		void setConiugate();
		
		double getX() { return v[0]; }
		double getY() { return v[1]; }
		double getZ() { return v[2]; }
		double getW() { return a; }

		void setX(double _x) { v[0] = _x; }
		void setY(double _y) { v[1] = _y; }
		void setZ(double _z) { v[2] = _z; }
		void setW(double _w) { a    = _w; }   

	
	        std::string toString();
		
	private:
		std::vector<double> multiplyVec(std::vector<double> _v, double _b) const;
		double dotVec(std::vector<double> _v1, std::vector<double> _v2) const;
		std::vector<double> crossVec(std::vector<double> _v1, std::vector<double> _v2) const;
#ifdef __GSL__
		void getPrincipalAxes(std::vector<std::vector<double> > &mat);   // Replaces 'mat' with eigenvectors
#endif

		// Aren't these great variable names? .. Please change.
		double a;           // angle
		std::vector<double> v;   // axis

};

}

#endif
