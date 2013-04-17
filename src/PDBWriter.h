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
#include "AtomPointerVector.h"
#include "FormatConverter.h"

// STL Includes
#include <vector>



namespace MSL { 

class PDBWriter : public Writer {

	public:
		// Constructors/Destructors
		PDBWriter();
		PDBWriter(const std::string &_filename);
		virtual ~PDBWriter();

		// Get/Set

		// Member Functions
		bool write(std::vector<CartesianPoint> &_cv);
		bool write(AtomPointerVector &_av, bool _addTerm=true, bool _noHydrogens=false,bool _writeAsModel=false);
		void writeREMARKS();
		bool open();               // There is a default implementation
		bool open(const std::string &_filename); // There is a default implementation
		bool open(const std::string &_filename, int mode); // There is a default implementation
		bool open(std::stringstream &_ss);
		void close();

		void setConvertFormat(std::string _orig, std::string _target); // for example from PDB3 to CHARMM22

		// Operators
		friend PDBWriter& operator<<(PDBWriter& pdbWriter, std::vector<CartesianPoint> &_cv) { return pdbWriter;};

	protected:		
	private:
		FormatConverter fc;
		bool convert_flag;
		std::string origFormat;
		std::string targetFormat;

};

//Inlines go HERE
inline bool PDBWriter::open() {bool success = Writer::open(); if(success) writeREMARKS(); return success;}
inline bool PDBWriter::open(const std::string &_filename) {bool success = Writer::open(_filename); if(success) writeREMARKS(); return success;}
inline bool PDBWriter::open(const std::string &_filename, int mode) {bool success = Writer::open(_filename, mode); if(success) writeREMARKS(); return success;}
inline bool PDBWriter::open(std::stringstream &_ss) {fileHandler = stringstyle; bool success = Writer::open(_ss); if(success) writeREMARKS(); return success;}
inline void PDBWriter::close() { std::string end = "END"; writeln(end); Writer::close(); }
inline void PDBWriter::setConvertFormat(std::string _orig, std::string _target) {
	if (_orig != "" and _target != "") {
		origFormat = _orig;
		targetFormat = _target;
		convert_flag = true;
		fc.setNamespaces(_orig, _target);
	} else {
		convert_flag = false;
	}
}

}

#endif
