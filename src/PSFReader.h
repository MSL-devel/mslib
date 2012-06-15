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

#ifndef PSFREADER_H
#define PSFREADER_H


//MSL Includes
#include "Reader.h"
#include "MslTools.h"


// Storage Includes
#include "AtomPointerVector.h"

// STL Includes
#include <vector>

/**
 * This class will provide an object which is able
 * to read in and interpret PSF files.
 */
namespace MSL { 
class PSFReader : public Reader {

	public:
		// Constructors/Destructors
		PSFReader();
		//PSFReader(const std::string &_filename);
		//PSFReader(const PSFReader & _reader);
		//PSFReader(std::stringstream &_stream);
		virtual ~PSFReader();

		// this function assigns coordinates from the atoms of the
		// PSF to an external AtomPointerVector as long as chainId, resnum, resname
		// and atom name are identical.  No errors are assigned for mismatches
		//void assignCoordinates(AtomPointerVector & _av);

		bool read();
		//bool read(std::string &_inputString);

		// Get/Set

		//AtomPointerVector & getAtomPointers() {return atoms;};

		//size_t size() const {return atoms.size();};

		//Atom * operator[](size_t _n) {return atoms[_n];};

		//void reset();

	protected:		
	private:
		//void deletePointers();
		AtomPointerVector atoms;

		std::vector<std::vector<int> > bonds;     //inner std::vector size 2
		std::vector<std::vector<int> > angles;    //inner std::vector size 3
		std::vector<std::vector<int> > dihedrals; //inner std::vector size 4
		std::vector<std::vector<int> > impropers; //inner std::vector size 4

		
};
inline PSFReader::PSFReader() : Reader() {}
inline PSFReader::~PSFReader() {}
}

#endif
