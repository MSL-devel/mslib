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

#ifndef SPHERICALPOINT_H
#define SPHERICALPOINT_H


// STL Includes
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>

// MSL Includes
#include "math.h"
#include "Real.h"
#include "Matrix.h"


// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
//#include <boost/serialization/is_abstract.hpp>
#endif


// Forward Declarations
namespace MSL { 
class SphericalGeometry;

// Namespaces


class SphericalPoint {
	public:
		SphericalPoint();
		SphericalPoint(double _radius, double _sigma, double _theta);
		~SphericalPoint();

		void setCoor(double _radius, double _sigma, double _theta);
		
		double getRadius();
		double getSigma();
		double getTheta();

	private:
		double radius;
		double sigma;
		double theta;

};

inline SphericalPoint::SphericalPoint(){ setCoor(0.0,0.0,0.0); }
inline SphericalPoint::~SphericalPoint(){ }
inline SphericalPoint::SphericalPoint(double _radius, double _sigma, double _theta){
	setCoor(_radius,_sigma,_theta);
}
inline void SphericalPoint::setCoor(double _radius, double _sigma, double _theta){
	radius = _radius;
	sigma  = _sigma;
	theta  = _theta;
}
inline double SphericalPoint::getRadius() { return radius; }
inline double SphericalPoint::getSigma() { return sigma; }
inline double SphericalPoint::getTheta() { return theta; }
}

#endif
