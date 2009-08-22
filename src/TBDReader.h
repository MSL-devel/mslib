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

#ifndef TBDREADER_H
#define TBDREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"


// Storage Formats
class TwoBodyDistanceDependentPotentialTable;


// STL Includes
#include <vector>
using namespace std;

/**
 * This class will provide an object which is able
 * to read in and interpret TBD files.
 */
class TBDReader : public Reader {

	public:
		// Constructors/Destructors
		TBDReader();
		TBDReader(const string &_filename);
		TBDReader(const TBDReader & _reader);
		virtual ~TBDReader();


		bool read();
		bool read(TwoBodyDistanceDependentPotentialTable *_tbd);
		


		void reset();


	protected:		
	private:
		void parseTBDLine(string _tbdline);

};

//Inlines go HERE
/**
 * Simple constructor.
 */
inline TBDReader::TBDReader() : Reader() {}
/**
 * With this constructor the user specifies the filename
 * of the TBD to be read.
 *
 * @param _filename  The name of the TBD file to be read.
 */
inline TBDReader::TBDReader(const string &_filename) : Reader(_filename) { }
/**
 * A copy constructor.  All of the atoms from the given TBDReader are
 * copied into the new TBDReader.
 *
 * @param _reader The TBDReader to be copied.
 */
inline TBDReader::TBDReader(const TBDReader & _reader) {

}


/**
 * The deconstructor.  All data will be deleted, so any Atom pointers
 * that were previously saved off will no longer be valid after the TBDReader
 * object has been destroyed.
 */
inline TBDReader::~TBDReader() { }
/**
 * This method will delete all data held in the TBDReader.  All
 * Atom pointers that were previously saved off will no longer be valid
 * after this reset.
 */
inline void TBDReader::reset() { }




#endif
