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
	double calcDistance1();
	double calcDistance2();
	double calcAngle1();
	double calcAngle2();
	double calcAngle3();
	double calcAngle4();
	double calcTorsion1();
	double calcTorsion2();

	double getDistance1();
	double getDistance2();
	double getAngle1();
	double getAngle2();
	double getAngle3();
	double getAngle4();
	double getTorsion1();
	double getTorsion2();

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
	
	double distance1;
	double distance2;
	double angle1;
	double angle2;
	double angle3;
	double angle4;
	double torsion1;
	double torsion2;
	
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

			ar & make_nvp("distance1",distance1);
			ar & make_nvp("distance2",distance2);
			ar & make_nvp("angle1",angle1);
			ar & make_nvp("angle2",angle2);
			ar & make_nvp("angle3",angle3);
			ar & make_nvp("angle4",angle4);
			ar & make_nvp("torsion1",torsion1);
			ar & make_nvp("torsion2",torsion2);

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
