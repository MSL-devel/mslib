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
#ifndef CRYSTALLATTICE_H
#define CRYSTALLATTICE_H

// STL Includes
#include <string>
#include <map>

// MSL Includes
#include "PDBReader.h"
#include "PDBWriter.h"

// Namespaces


namespace MSL { 
class CrystalLattice {

	public:
		CrystalLattice();
		CrystalLattice(std::string _pdbFile);
		~CrystalLattice();
		

		std::string getPdbFile();
		void setPdbFile(std::string _pdbFile);

		void generateCrystal();

		void writeCrystalUnits(std::string _pathAndPrefix, bool closeContactsOnly=true, bool singleFile=false, std::string _renameChains="", bool _nmrStyleFile=false);
		

	private:
		void readPdb();
		void copyAtoms(AtomPointerVector * _atoms, AtomPointerVector *newAts);

		std::string pdbFile;
		bool pdbFileRead;

		PDBReader pin;
		PDBWriter pout;
		std::map<std::string, AtomPointerVector *> crystalUnits;


};

inline std::string CrystalLattice::getPdbFile() { return pdbFile; }
inline void CrystalLattice::setPdbFile(std::string _pdbFile) { pdbFile = _pdbFile;}
}

#endif
