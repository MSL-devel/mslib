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

#ifndef CRDWRITER_H
#define CRDWRITER_H
/*
  This class handles writing of CRD files.
 */

// MSL Includes
#include "Writer.h"
#include "MslTools.h"
#include "CRDFormat.h"

// Storage formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"

// STL Includes
#include <vector>

namespace MSL { 
class CRDWriter : public Writer {

	public:
		// Constructors/Destructors
		CRDWriter();
		CRDWriter(const std::string &_filename);
		virtual ~CRDWriter();

		// Get/Set

		// Member Functions
		bool write(AtomPointerVector &_av, bool _writeRemarks=true);
		void writeREMARKS();
		bool open();               // There is a default implementation
		bool open(const std::string &_filename); // There is a default implementation
		bool open(const std::string &_filename, int mode); // There is a default implementation
		bool open(std::stringstream &_ss);
		void close();

		// Operators
		friend CRDWriter& operator<<(CRDWriter& _crdWriter, std::vector<CartesianPoint> &_cv) { return _crdWriter;};

	protected:		
	private:


};

//Inlines go HERE
inline CRDWriter::CRDWriter() : Writer() {}
inline CRDWriter::CRDWriter(const std::string &_filename) : Writer(_filename) {}
inline CRDWriter::~CRDWriter() {}
inline bool CRDWriter::open() {bool success = Writer::open(); return success;}
inline bool CRDWriter::open(const std::string &_filename) {bool success = Writer::open(_filename); return success;}
inline bool CRDWriter::open(const std::string &_filename, int mode) {bool success = Writer::open(_filename, mode); return success;}
inline bool CRDWriter::open(std::stringstream &_ss) {fileHandler = stringstyle; bool success = Writer::open(_ss); return success;}
inline void CRDWriter::close() { Writer::close(); }

}

#endif
