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

#ifndef DISTANCEHASHING_H
#define DISTANCEHASHING_H

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


#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "VectorPair.h"
#include "System.h"
#include "Chain.h"
#include "PDBTopology.h"

using namespace MSL;
using namespace std;



namespace MSL { 
class VectorHashing {

	public:
		VectorHashing();
		VectorHashing(const VectorHashing &_copyThis);
		~VectorHashing();


		void setFilterVectorPairs(bool _flag);
		bool getFilterVectorPairs();

		bool filterVectorPair(VectorPair &_vp, std::string residueName1="", std::string residueName2="");
		bool addToVectorHash(System &_sys,string _id,bool _printIt=false);

		vector<map<string,vector<string> > > searchForVectorMatchAll(VectorHashing &_vh, int _numAcceptableEdges);

		string getHashKey(VectorPair &_vp);

	private:

		bool addToVectorHash(Chain &_ch,string _id,bool _printIt=false);
		bool getPositionId(VectorPair &_vp, string &_posId1, string &_posId2);

		// Store a list of position ids for each chain that has been added
		vector<vector<string> > positionIds;

		map<string,VectorPair> pairPositionHash;
		map<string,vector<VectorPair*> > geometricHash;


		struct CandidateCycle {
			/*
			  inner pairs:
			      PositionId,CanididatePositionId
			  middle pair: 
			     first  = inner pairs(PositionId,CandidatePositionId)
			     second = inner pairs(PositionId2,CandidatePositionId2)
			  outer pair:
			      first = middle pairs
			      second = vector of rotamer pairs to connect CanidatePositionId and CandidatePositionId2
			 */
			vector< pair< pair< pair<string, string> , pair<string, string> > , vector<string> > > cycles;
			
		};
		
		double distanceGridSize;
		double angleGridSize;
		double dihedralGridSize;
		
		bool filterVectorPairs;

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

			ar & make_nvp("positionIds",positionIds);
			ar & make_nvp("pairPositionHash",pairPositionHash);
			//ar & make_nvp("geometricHash",geometricHash);
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
