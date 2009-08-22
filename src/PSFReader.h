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

#ifndef PSFREADER_H
#define PSFREADER_H


//MSL Includes
#include "Reader.h"
#include "MslTools.h"


// Storage Includes
#include "AtomVector.h"

// STL Includes
#include <vector>
using namespace std;

/**
 * This class will provide an object which is able
 * to read in and interpret PSF files.
 */
class PSFReader : public Reader {

	public:
		// Constructors/Destructors
		PSFReader();
		//PSFReader(const string &_filename);
		//PSFReader(const PSFReader & _reader);
		//PSFReader(stringstream &_stream);
		virtual ~PSFReader();

		// this function assigns coordinates from the atoms of the
		// PSF to an external AtomVector as long as chainId, resnum, resname
		// and atom name are identical.  No errors are assigned for mismatches
		//void assignCoordinates(AtomVector & _av);

		bool read();
		//bool read(string &_inputString);

		// Get/Set

		//AtomVector & getAtoms() {return atoms;};

		//size_t size() const {return atoms.size();};

		//Atom * operator[](size_t _n) {return atoms[_n];};

		//void reset();

	protected:		
	private:
		//void deletePointers();
		AtomVector atoms;

		vector<vector<int> > bonds;     //inner vector size 2
		vector<vector<int> > angles;    //inner vector size 3
		vector<vector<int> > dihedrals; //inner vector size 4
		vector<vector<int> > impropers; //inner vector size 4

		
};
inline PSFReader::PSFReader() : Reader() {}
inline PSFReader::~PSFReader() {}
#endif
