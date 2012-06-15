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

#ifndef CARTESIANPOINT_H
#define CARTESIANPOINT_H


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
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#endif


// Forward Declarations
namespace MSL { 

class CartesianPoint {
	public:
		CartesianPoint();
		CartesianPoint(std::string _oXYZ);  // "X" makes (1,0,0), etc.
		CartesianPoint(Real _x, Real _y, Real _z);
		CartesianPoint(std::vector<Real> _vec); // std::vector must be of length 3, containing x, y, and z
		CartesianPoint(const CartesianPoint & _point); // copy constructor

		~CartesianPoint();  // deconstructor

		bool operator==(const CartesianPoint & _point) const { return ((x == _point.x) && (y == _point.y) && (z == _point.z)); }; // same point?
		bool operator!=(const CartesianPoint & _point) const { return ((x != _point.x) || (y != _point.y) || (z != _point.z)); }; // different point?
                
		void operator=(const CartesianPoint & _point) { x = _point.x; y = _point.y; z = _point.z; }; // assignment
		void operator+=(const CartesianPoint & _point) { x += _point.x; y += _point.y; z += _point.z; }; // add _point coordinates to this point
		void operator-=(const CartesianPoint &  _point) { x -= _point.x; y -= _point.y; z -= _point.z; }; // subtract _point coordinates to this point
		CartesianPoint operator- (const CartesianPoint &  _point) const { return CartesianPoint((x - _point.x), (y - _point.y), (z - _point.z)); };
		CartesianPoint operator+ (const CartesianPoint &  _point) const{ return CartesianPoint((x + _point.x), (y + _point.y), (z + _point.z)); };
		double operator* (const CartesianPoint &  _point) const { return ((x*_point.x)+(y*_point.y)+(z*_point.z)); };
		CartesianPoint operator* (double _factor) const { return CartesianPoint((x*_factor), (y*_factor), (z*_factor)); };
		CartesianPoint operator/ (double _factor) const;
		void operator*=(double _factor) { x = x*_factor; y = y*_factor; z = z*_factor; }; // multiply this point by _factor
		void operator/=(double _factor); // divide this point by _factor
		CartesianPoint operator* (const Matrix & _rotation) const;
		void operator*=(const Matrix & _rotation); // rotate this point by _rotation rotation matrix
		Real & operator[](size_t n);
		friend std::ostream & operator<<(std::ostream &_os, const CartesianPoint & _vec) {_os << _vec.toString(); return _os;};

		double distance(CartesianPoint _p) const;
		double distance2(CartesianPoint _p) const;
		double angle(CartesianPoint _p) const; // center is the origin
		double angle(CartesianPoint _center, CartesianPoint _third) const;
		double angleRadians(CartesianPoint _p) const; // center is the origin
		double angleRadians(CartesianPoint _center, CartesianPoint _third) const;
		double dihedral(CartesianPoint _second, CartesianPoint _third, CartesianPoint _fourth) const;
		double dihedralRadians(CartesianPoint _second, CartesianPoint _third, CartesianPoint _fourth) const;

		std::string toString() const { char c [100]; sprintf(c, "[%10.3f %10.3f %10.3f]", x, y, z); return (std::string)c; };

		Real getX() const {return x;};
		Real getY() const {return y;};
		Real getZ() const {return z;};
		void setX(Real _x) {x = _x;};
		void setY(Real _y) {y = _y;};
		void setZ(Real _z) {z = _z;};

		Real* getXptr() {return &x;};
		Real* getYptr() {return &y;};
		Real* getZptr() {return &z;};

		void setCoor(Real _x, Real _y, Real _z) {x = _x; y = _y; z = _z;}; // set coordinates
		void setCoor(std::vector<Real> _point); // std::vector must be of length 3, containing x, y, and z
		void setCoor(const CartesianPoint & _point) { x = _point.x; y = _point.y; z = _point.z; };
		std::vector<Real> getCoor() const;

		double length() const { return sqrt(*this * *this); };
		CartesianPoint cross(CartesianPoint _second) const;
		CartesianPoint getUnit() const { double dist = length(); return (*this/dist); };


	protected:

		Real x;
		Real y;
		Real z;


		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:
		void save_checkpoint(std::string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(std::string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}

	private:
		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			using boost::serialization::make_nvp;
			ar & make_nvp("x",x);
			ar & make_nvp("y",y);
			ar & make_nvp("z",z);
		}
#endif


};

// A custom compare function so that CartesianPoint can be used as the
// key in std::map.
class CartesianPointCompare {
public:
    bool operator()(const CartesianPoint &pt1, const CartesianPoint &pt2) {
        if (pt1.getX() < pt2.getX())
            return true;
        else if (pt1.getX() > pt2.getX())
            return false;

        if (pt1.getY() < pt2.getY())
            return true;
        else if (pt1.getY() > pt2.getY())
            return false;

        if (pt1.getZ() < pt2.getZ())
            return true;

        return false;
    }
};

}

#endif
