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



#ifndef VECTORPAIR_H
#define VECTORPAIR_H

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#endif


// STL Includes
#include <map>
#include <vector>
#include <string>
#include <algorithm>


#include "CartesianPoint.h"

using namespace MSL;
using namespace std;



namespace MSL { 
  class VectorPair {
    public:
        VectorPair();
        VectorPair(CartesianPoint &_a1,CartesianPoint &_a2,CartesianPoint &_b1,CartesianPoint &_b2,string _vectorAid="", string _vectorBid="");
        VectorPair(const VectorPair &_vp);
        VectorPair& operator=(const VectorPair &_vp);
	void copy(const VectorPair &_vp);

        ~VectorPair();


	void calcAll();
	double calcDistance();
	double calcAngle1();
	double calcAngle2();
	double calcTorsion();

	double getDistance();
	double getAngle1();
	double getAngle2();
	double getTorsion();

	bool operator< ( const VectorPair &rhs ) const;

	string getVectorAid();
	string getVectorBid();
	string getVectorPairId();
	string getVectorPairIdReverse();

	CartesianPoint getA1() const { return a1; }
	CartesianPoint getA2() const { return a2; }
	CartesianPoint getB1() const { return b1; }
	CartesianPoint getB2() const { return b2; }
	
	string toString() const;

	// Boost related, is ok if no BOOST libraries are being used, just a std::string.
	void setArchiveType(std::string _type) { archiveType = _type; }
	std::string getArchiveType() { return archiveType; }

	
    private:
	CartesianPoint a1;
	CartesianPoint a2;

	CartesianPoint b1;
	CartesianPoint b2;
	
	double distance;
	double angle1;
	double angle2;
	double torsion;
	
	string vectorAid;
	string vectorBid;

	std::string archiveType;


		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(std::string filename) const{

			if (archiveType == "binary"){
				std::ofstream fout(filename.c_str(),std::ios::binary);
				boost::archive::binary_oarchive oa(fout);
				oa << (*this);
			} else if (archiveType == "xml"){
				std::ofstream fout(filename.c_str());
				boost::archive::xml_oarchive oa(fout);
				oa << boost::serialization::make_nvp("Matrix",*this);
			} else {
				std::ofstream fout(filename.c_str());
				boost::archive::text_oarchive oa(fout);
				oa << (*this);
			}

		}

		void load_checkpoint(std::string filename){

			if (archiveType == "binary"){
				std::ifstream fin(filename.c_str(), std::ios::binary);
				boost::archive::binary_iarchive ia(fin);
				ia >> (*this);
			} else if (archiveType == "xml"){
				std::ifstream fin(filename.c_str());
				boost::archive::xml_iarchive ia(fin);
				ia >> boost::serialization::make_nvp("Matrix",*this);
			} else {
				std::ifstream fin(filename.c_str());
				boost::archive::text_iarchive ia(fin);
				ia >> (*this);
			}
		}


	private:
		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			using boost::serialization::make_nvp;

			//ar & make_nvp("distanceData",distanceData);
			//			ar & make_nvp("atoms",atoms);

			ar & make_nvp("a1",a1);
			ar & make_nvp("a2",a2);
			ar & make_nvp("b1",b1);
			ar & make_nvp("b2",b2);

			ar & make_nvp("distance",distance);
			ar & make_nvp("angle1",angle1);
			ar & make_nvp("angle2",angle2);
			ar & make_nvp("torsion",torsion);

			ar & make_nvp("vectorAid",vectorAid);
			ar & make_nvp("vectorBid",vectorBid);

		}
#else
	public:
		void save_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF SAVE_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
		void load_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF LOAD_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}		
#endif

  };
}
#endif
