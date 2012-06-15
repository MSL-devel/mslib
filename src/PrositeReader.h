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

#ifndef PROSITEREADER_H
#define PROSITEREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

// Storage Formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"


// STL Includes
#include <vector>
#include <map>

/**
 * This class will provide an object which is able
 * to read in and interpret CRD files.
 */
namespace MSL { 
class PrositeReader : public Reader {

	public:
		// Constructors/Destructors
		PrositeReader();
		PrositeReader(const std::string &_filename);
		PrositeReader(const PrositeReader & _reader);
		PrositeReader(std::stringstream &_stream);
		virtual ~PrositeReader();
		
		bool read();
		bool read(std::string &_inputString);

		std::map<std::string,std::map<std::string, std::string> >& getPrositeData();

	protected:		
	private:
		void setup();

		// ID,DE,PA  (pattern id, pattern description, pattern regex string)
		std::map<std::string, std::map<std::string, std::string> > prosite_data;

};

//Inlines go HERE
/**
 * Simple constructor.
 */
inline PrositeReader::PrositeReader() : Reader() { setup(); }
/**
 * With this constructor the user specifies the filename
 * of the Prosite to be read.
 *
 * @param _filename  The name of the Prosite file to be read.
 */
inline PrositeReader::PrositeReader(const std::string &_filename) : Reader(_filename) { setup();}

/**
 * A copy constructor.  All of the atoms from the given PrositeReader are
 * copied into the new PrositeReader.
 *
 * @param _reader The PrositeReader to be copied.
 */
inline PrositeReader::PrositeReader(const PrositeReader & _reader) {
  prosite_data = _reader.prosite_data;
}
/**
 * A constructor which will read input data from a std::stringstream.
 *
 * @param _ss The std::stringstream to get data from.
 */
inline PrositeReader::PrositeReader(std::stringstream &_ss) : Reader(_ss)     {
	setup();
	read();
}

inline void PrositeReader::setup() {

  prosite_data.clear();
}

/**
 * The deconstructor.  All data will be deleted, so any Atom pointers
 * that were previously saved off will no longer be valid after the PrositeReader
 * object has been destroyed.
 */
inline PrositeReader::~PrositeReader() { close();}





inline bool PrositeReader::read(std::string &_inputString){
	fileName = "string";
	stringStreamPtr = new std::stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}

 inline std::map<std::string,std::map<std::string, std::string> >& PrositeReader::getPrositeData(){
   return prosite_data;
 } 
}


#endif

