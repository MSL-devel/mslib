#ifndef DISTANCEMATRIXDATABASE_H
#define DISTANCEMATRIXDATABASE_H

// STL Includes
#include <iostream>
#include <string>
#include <map>
#include <math.h>

// MSL Includes
#include "DistanceMatrix.h"
#include "MslTools.h"


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

//class DistanceMatrix;
//class MatrixWindow;

using namespace std;
namespace MSL{
using namespace MslTools;

class DistanceMatrixDatabase {  

	public:
		DistanceMatrixDatabase();
		DistanceMatrixDatabase(const DistanceMatrixDatabase &_dm);
		~DistanceMatrixDatabase();

		void addDistanceMatrix(DistanceMatrix *_dm) { listDMs.push_back(_dm);}
		vector<DistanceMatrix *> & getDistanceMatrixList() { return listDMs; }

		void setName(string _name){ name = _name; }
		string getName() { return name;}

		void setArchiveType(string _archiveType){ archiveType = _archiveType; }
		string getArchiveType() { return archiveType;}
		
	private:
		string name;
		string archiveType;
		vector<DistanceMatrix *> listDMs;

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
				oa << boost::serialization::make_nvp("DistanceMatrixDatabase",*this);
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
				ia >> boost::serialization::make_nvp("DistanceMatrixDatabase",*this);
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
			ar & make_nvp("archiveType",archiveType);
			ar & make_nvp("listDMs",listDMs);

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
}
#endif
