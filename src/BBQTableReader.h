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

#ifndef BBQ_TABLE_READER_H
#define BBQ_TABLE_READER_H

//MSL Includes
#include "Reader.h"
#include "MslTools.h"

// Storage Formats
#include "BBQTable.h"

/**
 * This class will provide an object which is able
 * to read in and interpret BBQTable files.
 */
namespace MSL { 
class BBQTableReader : public Reader {
    public:
		// Constructors/Destructors
		BBQTableReader();
		BBQTableReader(const std::string &_filename);
		BBQTableReader(const BBQTableReader & _reader);
		BBQTableReader(std::stringstream &_stream);
		virtual ~BBQTableReader();

		bool read() { return read(bbqTable); };
                bool read(BBQTable &_bbqTable);

		BBQTable &getBBQTable()  { return bbqTable; }

		void reset();

	protected:
	private:
		BBQTable bbqTable;
};

//Inlines go HERE
/**
 * Simple constructor.
 */
inline BBQTableReader::BBQTableReader() : Reader() {}
/**
 * With this constructor, the user specifies the filename
 * of the BBQ Table to be read.
 *
 * @param _filename  The name of the BBQ Table file to be read.
 */
inline BBQTableReader::BBQTableReader(const std::string &_filename) : Reader(_filename) {}
/**
 * A copy constructor.  All of the data from the given BBQTableReader are
 * copied into the new BBQTableReader.
 *
 * @param _reader The BBQTableReader to be copied.
 */
inline BBQTableReader::BBQTableReader(const BBQTableReader & _reader) {  }
/**
 * A constructor which will read input data from a std::stringstream.
 *
 * @param _ss The std::stringstream to get data from.
 */
inline BBQTableReader::BBQTableReader(std::stringstream &_ss) : Reader(_ss)     {read();}

/**
 * The deconstructor.
 *
 */
inline BBQTableReader::~BBQTableReader() { close();}

}

#endif // BBQ_TABLE_READER_H
