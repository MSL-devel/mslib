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
#include "AtomVector.h"

// GSL Includes
#ifdef __GSL__
  #include <gsl/gsl_eigen.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
#endif


using namespace std;

class Quaternion {  

	public:
		Quaternion();
		Quaternion(double _a, double _v0, double _v1, double _v2);
		Quaternion(double _a, vector<double> _v);
 
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
	        friend ostream & operator<<(ostream &_os, Quaternion & _q) {_os << _q.toString(); return _os;};


		bool convertToRotationMatrix(Matrix &matrix);                             // Converts this quaternion to  a rotation matrix
		bool convertToQuaternion(Matrix &matrix);                                 // Makes a quaternion (this object) from a rotation matrix
									       	      									 
		bool makeQuaternion(CartesianPoint &axis, const double theta);            // Makes a quaternion (this object) from an axis and an angle
		bool makeQuaternion(AtomVector &_align,  AtomVector &_ref);

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

	
	        string toString();
		
	private:
		vector<double> multiplyVec(vector<double> _v, double _b) const;
		double dotVec(vector<double> _v1, vector<double> _v2) const;
		vector<double> crossVec(vector<double> _v1, vector<double> _v2) const;
		void getPrincipalAxes(vector<vector<double> > &mat);   // Replaces 'mat' with eigenvectors

		// Aren't these great variable names? .. Please change.
		double a;           // angle
		vector<double> v;   // axis

};

#endif
