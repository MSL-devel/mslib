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

#ifndef PHIPSIWRITER_H
#define PHIPSIWRITER_H
/*
  This class handles writing of PhiPsi files.
 */

// MSL Includes
#include "Writer.h"
#include "MslTools.h"
#include "PhiPsiStatistics.h"

// STL Includes
#include <vector>

namespace MSL { 
class PhiPsiWriter : public Writer {

	public:
		// Constructors/Destructors
		PhiPsiWriter();
		PhiPsiWriter(const std::string &_filename);
		virtual ~PhiPsiWriter();

		// Get/Set

		// Member Functions
		bool write(PhiPsiStatistics &_stat);

		bool open();               // There is a default implementation
		bool open(const std::string &_filename); // There is a default implementation
		bool open(const std::string &_filename, int mode); // There is a default implementation
		bool open(std::stringstream &_ss);
		void close();


	protected:		
	private:

};

//Inlines go HERE
inline PhiPsiWriter::PhiPsiWriter() : Writer() {}
inline PhiPsiWriter::PhiPsiWriter(const std::string &_filename) : Writer(_filename) {}
inline PhiPsiWriter::~PhiPsiWriter() {}
inline bool PhiPsiWriter::open() {bool success = Writer::open();  return success;}
inline bool PhiPsiWriter::open(const std::string &_filename) {bool success = Writer::open(_filename); return success;}
inline bool PhiPsiWriter::open(const std::string &_filename, int mode) {bool success = Writer::open(_filename, mode); return success;}
inline bool PhiPsiWriter::open(std::stringstream &_ss) {fileHandler = stringstyle; bool success = Writer::open(_ss); return success;}
inline void PhiPsiWriter::close() { Writer::close(); }

}

#endif
