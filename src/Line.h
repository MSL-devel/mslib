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

#ifndef LINE_H
#define LINE_H

// STL Includes
#include <string>
#include <vector>
#include <math.h>
#include <iostream>


// MSL Includes
#include "CartesianPoint.h"
#include "CartesianGeometry.h"
#include "Matrix.h"


// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#endif


using namespace std;

/*! \brief Line in 3D space (a center and a direction vector) */

class Line {

	public:
		Line();
		Line(double _x1, double _y1, double _z1, double _x2, double _y2, double _z2);
		Line(const CartesianPoint & _center, const CartesianPoint & _direction);
		Line(const Line * _line);
		Line(const Line & _line);
		~Line();

		void  operator=(const Line& _line);
		void  operator=(Line*  _line);
		bool  operator!=(const Line& _line) const;
		bool  operator!=(const Line* _line) const;
		bool  operator==(const Line& _line) const;
		bool  operator==(const Line *_line) const;
		friend ostream & operator<<(ostream &_os, Line & _line3D) {_os << _line3D.toString(); return _os;};

		CartesianPoint getCenter() const;
		CartesianPoint getDirection() const;
		void setCenter(const CartesianPoint & theCenter);
		void setDirection(const CartesianPoint & theDirection);

		bool isParallel(const Line & line) const;
		bool isOpposite(const Line & line) const;


		double distance(const CartesianPoint & point) const;
		double distance(const Line & line) const;

		CartesianPoint projection(const CartesianPoint & point) const;
		CartesianPoint pointOfMinDistanceFromLine(const Line & line) const;
		CartesianPoint pointOfMinDistanceToLine(const Line & line) const;
		double segmentDistance(const Line & _otherSegment); 
		double segmentDistance(const CartesianPoint & _vec); 
		double segmentDihedral(const Line & _otherSegment);
		vector<CartesianPoint> segmentsClosestPoints(const Line & _otherSegment);
	

		double angle(const Line & line) const;
		double angle(const CartesianPoint & direction) const;
		
		void translate(const CartesianPoint & vec);
		void xRotate(double degrees, const CartesianPoint & rotationCenter);
		void yRotate(double degrees, const CartesianPoint & rotationCenter);
		void zRotate(double degrees, const CartesianPoint & rotationCenter);
		void rotateAroundAxis(double degrees, const CartesianPoint & center, const CartesianPoint & axis);
		void rotateAroundAxis(double degrees, const Line & line);

		string toString() ;
		
		void setName(string _lineName) { name = _lineName; }
		string getName() { return name; }

		void setColor(vector<double> &_vec) { color = _vec; }
		void setColor(double _red, double _green, double _blue);
		void setColor(double _red, double _green, double _blue, double _red2, double _green2, double _blue2);
		CartesianPoint getColor() { return color; }

		void setOutputFormat(string _format) { outputFormat = _format; }
		string getOutputFormat() { return outputFormat; }

	protected:
		vector<CartesianPoint> pointsOfMinDistance(const Line & line) const;
		vector<double> solve3x3System(vector<vector<double> > sys) const;
		CartesianPoint closestPointWithinSegment(CartesianPoint _point) const;

		CartesianPoint center;
		CartesianPoint direction;

		string name;
		vector<double> color;
		string outputFormat;


		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(string filename) const{
			std::ofstream fout(filename.c_str());
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & center;
			ar & direction;
			ar & name;
			ar & color;
			ar & outputFormat;
		}
#endif

};

#endif
