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

#ifndef ATOMVECTOR_H
#define ATOMVECTOR_H

// STL Includes
#include <vector>
#include <sstream>
#include <sys/types.h>

// MSL Includes
#include "Atom.h"
#include "CartesianPoint.h"

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

// Forward Declarations
class Atom;


class AtomPointerVector : public vector<Atom *> {
	public:
		AtomPointerVector();
		AtomPointerVector(unsigned int _size, Atom * _pointer=NULL);
		AtomPointerVector(const AtomPointerVector & _atoms);
		~AtomPointerVector();

		void operator=(const AtomPointerVector & _atoms);
		AtomPointerVector operator+(const AtomPointerVector & _atoms);
		AtomPointerVector operator-(const AtomPointerVector & _atoms);
		void operator+=(const AtomPointerVector & _atoms);
		void operator-=(const AtomPointerVector & _atoms);

		Atom & operator()(size_t _n);

		friend ostream & operator<<(ostream &_os, const AtomPointerVector & _av)  {_os << _av.toString(); return _os;};

		// Setter/Getter Functions
		void   setName(string _name);
		string getName() const;
		CartesianPoint getGeometricCenter() const;

		// Geometric Center
	        void updateGeometricCenter(unsigned int _updateStamp=0);

		double rmsd(const AtomPointerVector &_av) const;

	        void translate(double _x, double _y, double _z);         
	        void translate(const CartesianPoint &_vec);         
		void rotate(const Matrix &_rotMat);


		/***************************************************
		 *  Saved sets:
		 *
		 *  coordinates can be saved to named buffers, and copied back
		 *  from them
		 *
		 *  The difference between saved coord and alt coor is that
		 *  the saved coord are never active, they can only be
		 *  used to store and copy back coordinates
		 ***************************************************/
		void saveCoor(string _coordName);
		bool applySavedCoor(string _coordName);
		void clearSavedCoor();		

		//void addAtomRanking(string _atomKey, double _val);
		//void sort();

		string toString() const;

		void deletePointers();


		// Boost related, is ok if no BOOST libraries are being used, just a string.
		void setArchiveType(string _type) { archiveType = _type; }
		string getArchiveType() { return archiveType; }

	private:
		string name;
		CartesianPoint geometricCenter;	

//		map<string, double> atomRankings;
		unsigned int updateStamp;

		string archiveType;



		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(string filename) const{

			if (archiveType == "binary"){
				std::ofstream fout(filename.c_str(),std::ios::binary);
				boost::archive::binary_oarchive oa(fout);
				oa << (*this);
			} else if (archiveType == "xml"){
				std::ofstream fout(filename.c_str());
				boost::archive::xml_oarchive oa(fout);
				oa << boost::serialization::make_nvp("AtomPointerVector",*this);
			} else {
				std::ofstream fout(filename.c_str());
				boost::archive::text_oarchive oa(fout);
				oa << (*this);
			}

		}

		void load_checkpoint(string filename){

			if (archiveType == "binary"){
				std::ifstream fin(filename.c_str(), std::ios::binary);
				boost::archive::binary_iarchive ia(fin);
				ia >> (*this);
			} else if (archiveType == "xml"){
				std::ifstream fin(filename.c_str());
				boost::archive::xml_iarchive ia(fin);
				ia >> boost::serialization::make_nvp("AtomPointerVector",*this);
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
			//ar & boost::serialization::base_object<vector<Atom *> >(*this);
			//ar & geometricCenter;

			ar & make_nvp("name",name);
			ar & make_nvp("atoms",boost::serialization::base_object<vector<Atom *> >(*this));
		}
#else
	public:
		void save_checkpoint(string filename) const{
			cout << "NO IMPLEMENTATION OF SAVE_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
		void load_checkpoint(string filename) const{
			cout << "NO IMPLEMENTATION OF LOAD_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
#endif
};

// INLINE FUNCTIONS
inline Atom & AtomPointerVector::operator()(size_t _n) { return *((*this)[_n]); }
inline void AtomPointerVector::setName(string _name) { name = _name; }
inline string AtomPointerVector::getName() const { return name;  }
inline CartesianPoint AtomPointerVector::getGeometricCenter() const { return geometricCenter; }
inline string AtomPointerVector::toString() const {
	stringstream ss;
	for (AtomPointerVector::const_iterator k = begin(); k!=end(); k++) {
		ss << **k << endl;
	}
	/*
	for (uint i = 0; i < (*this).size();i++){
		ss << "I: "<<i<<endl;
		ss << (*this)(i)<<endl;
	}
	*/
	return ss.str();
}
#endif
