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

#ifndef MOLECULEINTERFACEDATABASE
#define MOLECULEINTERFACEDATABASE

// MSL Includes
#include "InterfaceResidueDescriptor.h"

// STL Includes
#include <vector>
#include <iostream>

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
#include <boost/serialization/vector.hpp>
#endif



/**
 * This class holds a database of molecular interactions.
 * These interactions are held in objects of the 
 * InterfaceResidueDescriptor type.  File I/O is 
 * provided so that a database can be read in as
 * a text or binary format, and also saved out to file.
 *
 * @see InterfaceResidueDescriptor
 */
namespace MSL { 
class MoleculeInterfaceDatabase {
	public:
		MoleculeInterfaceDatabase();
		MoleculeInterfaceDatabase(std::string _flatFile);
		MoleculeInterfaceDatabase(MoleculeInterfaceDatabase &_pid);
		~MoleculeInterfaceDatabase();

		void operator=(MoleculeInterfaceDatabase &_pid);

		void setDatabaseFile(std::string _dbFile);
		std::string getDatabaseFile();

		//void addInterafaceResidueDescriptor(InterfaceResidueDescriptor &_ird);
		void addInterafaceResidueDescriptor(InterfaceResidueDescriptor *_ird);
		std::vector<InterfaceResidueDescriptor *> & getInterfaceResidueDescriptors();

		size_t size() const;
		InterfaceResidueDescriptor* operator[](size_t _n);

		// MOVED HERE FROM BOOST PRIVATE TO ALLOW COMPILATION WITHOUT BOOST
		void setArchiveType(std::string _type) { archiveType = _type; }
		std::string getArchiveType() { return archiveType; }


	private:

		void copy(MoleculeInterfaceDatabase &_pid);
		std::string databaseFile;
		std::vector<InterfaceResidueDescriptor *> residueDescriptors;

		// MOVED HERE FROM BOOST PRIVATE TO ALLOW COMPILATION WITHOUT BOOST
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
				oa << boost::serialization::make_nvp("MoleculeInterfaceDatabase",*this);
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
				ia >> boost::serialization::make_nvp("MoleculeInterfaceDatabase",*this);
			} else {
				std::ifstream fin(filename.c_str());
				boost::archive::text_iarchive ia(fin);
				ia >> (*this);
			}

		}



	private:
		friend class boost::serialization::access;		

		template<class Archive> void serialize(Archive & ar, const unsigned int version){

//			if (archiveType == "xml"){
				using boost::serialization::make_nvp;
				ar & make_nvp("DatabaseFile",databaseFile);
				ar & make_nvp("ResidueDescriptors",residueDescriptors);
/* 			} else { */
/* 				ar & databaseFile; */
/* 				ar & residueDescriptors; */
/* 			} */
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

// INLINES
/**
 * The basic constructor.
 */
inline MoleculeInterfaceDatabase::MoleculeInterfaceDatabase() { archiveType = "text";}
/**
 * A constructor that allows the user to specify a 
 * text file which should be read to populate the database.
 * 
 * @param _flatFile        Th e name of the text file 
 */
inline MoleculeInterfaceDatabase::MoleculeInterfaceDatabase(std::string _flatFile) { databaseFile = _flatFile; archiveType = "text";}
/**
 * A constructor that will clone the information found
 * in another MoleculeInterfaceDatabase.
 *
 * @param _pid      The MoleculeInterfaceDatabase to clone.
 */
inline MoleculeInterfaceDatabase::MoleculeInterfaceDatabase(MoleculeInterfaceDatabase &_pid) { copy(_pid);}
/**
 * Overload of the = operator so that we can easily clone
 * MoleculInterfaceDatabases.
 */
inline void MoleculeInterfaceDatabase::operator=(MoleculeInterfaceDatabase &_pid) { copy(_pid);}
/**
 * This method allows the user to set the filename
 * that should be used for this database.  This can be 
 * used instead of the constructor which takes the file name
 * as input.
 *
 *@param _dbFile       The name of the database file.
 */
inline void MoleculeInterfaceDatabase::setDatabaseFile(std::string _dbFile) { databaseFile = _dbFile; }
/**
 * This method will return the name of the database file.
 */
inline std::string MoleculeInterfaceDatabase::getDatabaseFile() { return databaseFile; }
/**
 * This method will return a std::vector of all of the 
 * interactions found in the database.
 *
 * @see InterfaceResidueDescriptor
 */
inline std::vector<InterfaceResidueDescriptor *> & MoleculeInterfaceDatabase::getInterfaceResidueDescriptors() { return residueDescriptors; }
/**
 * This method will return the number of entries
 * in the database.
 */
inline size_t  MoleculeInterfaceDatabase::size() const { return residueDescriptors.size(); }
/** 
 * Overload of the [] operator.  This allows us to retrieve 
 * the interaction in our database at a given index.
 *
 * @param _n     The index of the interaction in our database to be retrieved.
 */
inline InterfaceResidueDescriptor* MoleculeInterfaceDatabase::operator[](size_t _n){ 
    /// @todo ADD ERROR CHECK
	return residueDescriptors[_n]; 
}
}

#endif
