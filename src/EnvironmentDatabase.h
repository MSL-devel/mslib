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

#ifndef ENVIRONMENTDATABASE_H
#define ENVIRONMENTDATABASE_H


// MSL Includes
#include "System.h"
#include "EnvironmentDescriptor.h"

// STL Includes

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#endif



namespace MSL { 
class EnvironmentDatabase {
	public:

		EnvironmentDatabase();
		EnvironmentDatabase(std::string _type);
		EnvironmentDatabase(EnvironmentDatabase &_ed);
		~EnvironmentDatabase();	

		void operator=(EnvironmentDatabase &_ed);

		void createDatabase(System &_sys, std::string _systemName);

		std::vector<EnvironmentDescriptor*>& getAllDescriptors();
		int getNumberDescriptors();

		bool searchForEnvironment(EnvironmentDescriptor &_ed,std::string _envType);
		std::vector<EnvironmentDescriptor*>& getSearchResults();

		

	private:
		void copy(EnvironmentDatabase &_ed);
		
		std::vector<EnvironmentDescriptor *> descriptors;
		std::map<std::string, std::vector<EnvironmentDescriptor *> > lookupTable;
		std::string type;

		std::string lastSearchedKey;

		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__

	public:

		void save_checkpoint(std::string filename) const{
			std::ofstream fout(filename.c_str());
			//boost::archive::binary_oarchive oa(fout);
			boost::archive::text_oarchive oa(fout);
			oa << (*this);
		}

		void load_checkpoint(std::string filename){
			std::ifstream fin(filename.c_str(), std::ios::binary);
			//boost::archive::binary_iarchive ia(fin);
			boost::archive::text_iarchive ia(fin);
			ia >> (*this);
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){
			ar & descriptors;
			ar & lookupTable;
			ar & type;
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

inline int EnvironmentDatabase::getNumberDescriptors() { return descriptors.size(); }

}

#endif
