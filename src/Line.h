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



/*! \brief Line in 3D space (a center and a direction std::vector) */

namespace MSL { 
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
		friend std::ostream & operator<<(std::ostream &_os, Line & _line3D) {_os << _line3D.toString(); return _os;};

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
		std::vector<CartesianPoint> segmentsClosestPoints(const Line & _otherSegment);
	

		double angle(const Line & line) const;
		double angle(const CartesianPoint & direction) const;
		
		void translate(const CartesianPoint & vec);
		void xRotate(double degrees, const CartesianPoint & rotationCenter);
		void yRotate(double degrees, const CartesianPoint & rotationCenter);
		void zRotate(double degrees, const CartesianPoint & rotationCenter);
		void rotateAroundAxis(double degrees, const CartesianPoint & center, const CartesianPoint & axis);
		void rotateAroundAxis(double degrees, const Line & line);

		std::string toString() ;
		
		void setName(std::string _lineName) { name = _lineName; }
		std::string getName() { return name; }

		void setColor(std::vector<double> &_vec) { color = _vec; }
		void setColor(double _red, double _green, double _blue);
		void setColor(double _red, double _green, double _blue, double _red2, double _green2, double _blue2);
		CartesianPoint getColor() { return color; }

		void setOutputFormat(std::string _format) { outputFormat = _format; }
		std::string getOutputFormat() { return outputFormat; }

	protected:
		std::vector<CartesianPoint> pointsOfMinDistance(const Line & line) const;
		std::vector<double> solve3x3System(std::vector<std::vector<double> > sys) const;
		CartesianPoint closestPointWithinSegment(CartesianPoint _point) const;

		CartesianPoint center;
		CartesianPoint direction;

		std::string name;
		std::vector<double> color;
		std::string outputFormat;


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
			ar & center;
			ar & direction;
			ar & name;
			ar & color;
			ar & outputFormat;
		}
#endif

};

}

#endif
