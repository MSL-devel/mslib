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

#ifndef PDBWRITER_H
#define PDBWRITER_H
/*
  This class handles writing of PDB files.
 */

// MSL Includes
#include "Writer.h"
#include "MslTools.h"
#include "PDBFormat.h"

// Storage formats
#include "CartesianPoint.h"
#include "AtomVector.h"

// STL Includes
#include <vector>
using namespace std;

class PDBWriter : public Writer {

	public:
		// Constructors/Destructors
		PDBWriter();
		PDBWriter(const string &_filename);
		virtual ~PDBWriter();

		// Get/Set

		// Member Functions
		bool write(vector<CartesianPoint> &_cv);
		bool write(AtomVector &_av, bool _addTerm=true, bool _noHydrogens=false,bool _writeAsModel=false);
		void writeREMARKS();
		bool open();               // There is a default implementation
		bool open(const string &_filename); // There is a default implementation
		bool open(const string &_filename, int mode); // There is a default implementation
		bool open(stringstream &_ss);
		void close();

		// Operators
		friend PDBWriter& operator<<(PDBWriter& pdbWriter, vector<CartesianPoint> &_cv) { return pdbWriter;};

	protected:		
	private:


};

//Inlines go HERE
		inline PDBWriter::PDBWriter() : Writer() {}
		inline PDBWriter::PDBWriter(const string &_filename) : Writer(_filename) {}
		inline PDBWriter::~PDBWriter() {}
        inline bool PDBWriter::open() {bool success = Writer::open(); if(success) writeREMARKS(); return success;}
        inline bool PDBWriter::open(const string &_filename) {bool success = Writer::open(_filename); if(success) writeREMARKS(); return success;}
        inline bool PDBWriter::open(const string &_filename, int mode) {bool success = Writer::open(_filename, mode); if(success) writeREMARKS(); return success;}
        inline bool PDBWriter::open(stringstream &_ss) {bool success = Writer::open(_ss); if(success) writeREMARKS(); return success;}
        inline void PDBWriter::close() { string end = "END"; writeln(end); Writer::close(); }

#endif
