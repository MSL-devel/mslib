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

#ifndef READER_H
#define READER_H

/*
  A simple base class for readers.
 */

// MSL Includes
#include "File.h"
#include "MslTools.h"

// STL Includes
#include <iostream>
#include <vector>

namespace MSL { 
class Reader : public File {

	public:
		Reader();
		Reader(const std::string &_filename);
		Reader(const std::string &_filename, int &_mode);
		Reader(std::stringstream &_ss);
		Reader(const Reader &_anotherReader);
		std::vector<std::string> getAllLines();

		void operator=(const Reader &_anotherReader);
		
		virtual ~Reader() ;


		virtual bool read(){ std::cout <<"Reader::read() does nothing\n"; return false;}; //  default implementation.
		virtual bool read(std::string &_inputString);
		std::string getLine();
	protected:		
	private:
		void copy(const Reader &_anotherReader);
		
		
};

//INLINES GO HERE
inline Reader::Reader():File((const std::string)"",0) {}
inline Reader::Reader(const std::string &_filename) : File(_filename,0) {}
inline Reader::Reader(const std::string &_filename, int &_mode) : File(_filename,0) {}
inline Reader::Reader(std::stringstream &_ss) : File(_ss) {}
inline Reader::Reader(const Reader &_anotherReader) : File(_anotherReader) { copy(_anotherReader); }
inline void Reader::operator=(const Reader &_anotherReader) { copy(_anotherReader); }


inline Reader::~Reader() {}
}

#endif
