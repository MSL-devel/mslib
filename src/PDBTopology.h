/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2011 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

#ifndef PDBTOPOLOGY_H
#define PDBTOPOLOGY_H

// STL Includes
#include <vector>
#include <map>

//MSL Includes
#include "CharmmTopologyReader.h"
#include "RotamerLibrary.h"
#include "MslTools.h"
#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "Residue.h"

namespace MSL { 

class PDBTopology {

	public:
		PDBTopology();
		PDBTopology(const PDBTopology & _top);
		~PDBTopology();

		void operator=(const PDBTopology & _top);

		bool readCharmmTopology(std::string _charmmToplogyFile);
		bool readRotamerLibrary(std::string _rotamerLibrary);

		bool residueExists(std::string _name);

		AtomContainer getResidue(std::string _identityId);
		AtomContainer getResidue(std::string _identityId, AtomPointerVector &_backboneSeedAtoms, int _numRotamers);

		static Atom* getPseudoCbeta(Residue &_glycine);

		void reset();

		std::map<std::string, std::vector< std::vector<std::string> > > & getChis();
		
	        void setAddAtomsFromRotLib(bool _flag);
		bool getAddAtomsFromRotLib();

		AtomPointerVector getBackboneAtoms(Residue &_res);
	private:
		void setup();
		void copy(const PDBTopology & _top);

		void buildRotamers(AtomContainer &_newResidue, std::string _resName, int _numRotamers);

		std::map<std::string,std::vector<std::string> > atoms;

		std::map<std::string,bool> backboneAtoms;

		// Stores atom names for each residue, each chi angle
		std::map<std::string, std::vector< std::vector<std::string> > > chis;

		CharmmTopologyReader charmmReader;
		RotamerLibrary rotlib;

		bool addAtomsFromRotLib;
		

};
inline void PDBTopology::setAddAtomsFromRotLib(bool _flag) { addAtomsFromRotLib = _flag;}
inline bool PDBTopology::getAddAtomsFromRotLib() { return addAtomsFromRotLib; }
inline std::map<std::string, std::vector< std::vector<std::string> > > & PDBTopology::getChis() { return chis; }

}

#endif
