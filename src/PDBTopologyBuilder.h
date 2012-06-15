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


#ifndef PDBTOPOLOGYBUILDER_H
#define PDBTOPOLOGYBUILDER_H

#include <iostream>
#include <vector>

#include "System.h"
#include "PDBReader.h"
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "PolymerSequence.h"
//#include "RandomNumberGenerator.h"


namespace MSL { 
class PDBTopologyBuilder {
	public:
		PDBTopologyBuilder();
		PDBTopologyBuilder(System & _system, std::string _topologyFile);
		PDBTopologyBuilder(const PDBTopologyBuilder & _sysBuild);
		~PDBTopologyBuilder();

		void operator=(const PDBTopologyBuilder & _sysBuild);

		void setSystem(System & _system);

		bool readTopology(std::string _topologyFile);

		bool buildSystem(const PolymerSequence & _sequence);
		bool buildSystemFromPDB(std::string _fileName); // build from a PDB in CHARMM name format

		/**************************************************
		 * Add one or more new identities to a position,
		 * if a series of backbone atoms are specified
		 * (_bbAtoms) then the coordinates are passed to
		 * the atoms with the same name in the new residues
		 * to preserve the position of the backbone
		 **************************************************/
		bool addIdentity(std::string _positionId, std::string _resName, std::vector<std::string> _bbAtoms=std::vector<std::string>()); // id "A,37"
		bool addIdentity(std::string _positionId, const std::vector<std::string> & _resNames, std::vector<std::string> _bbAtoms=std::vector<std::string>());
		bool addIdentity(Position & _pos, std::string _resName, std::vector<std::string> _bbAtoms=std::vector<std::string>());
		bool addIdentity(Position & _pos, const std::vector<std::string> & _resNames, std::vector<std::string> _bbAtoms=std::vector<std::string>());
		// same functions but the bb atoms are passed as a space-separated list such as "N CA C O H"
		bool addIdentity(std::string _positionId, std::string _resName, std::string _bbAtoms="N CA C O H"); // id "A,37"
		bool addIdentity(std::string _positionId, const std::vector<std::string> & _resNames, std::string _bbAtoms="N CA C O H");
		bool addIdentity(Position & _pos, std::string _resName, std::string _bbAtoms="N CA C O H");
		bool addIdentity(Position & _pos, const std::vector<std::string> & _resNames, std::string _bbAtoms="N CA C O H");
		
		CharmmTopologyReader * getCharmmTopologyReader();

		bool fail() const; // return false if reading toppar failed

	private:
		void setup();
		void copy(const PDBTopologyBuilder & _sysBuild);
		void deletePointers();
		void reset();
		void getAtomPointersFromMulti(std::string _name, std::vector<Atom*> & _out, std::vector<CharmmTopologyResidue*> & _position, std::vector<std::map<std::string, Atom*> > & _atomMap);
		std::vector<Atom*> getAtomPointers(std::string _name, std::vector<std::vector<std::vector<CharmmTopologyResidue*> > >::iterator & _chItr, std::vector<std::vector<CharmmTopologyResidue*> >::iterator & _posItr, std::vector<CharmmTopologyResidue*>::iterator & _idItr);

		std::vector<std::vector<std::vector<CharmmTopologyResidue*> > > polymerDefi;
		std::vector<std::vector<std::vector<std::map<std::string, Atom*> > > > atomMap;

		System * pSystem;

		CharmmTopologyReader * pTopReader;

		bool fail_flag;
		bool wasBuilt_flag; // if a System is built from PDBTopology this is true

};
inline void PDBTopologyBuilder::setSystem(System & _system) {
	reset();
	pSystem = &_system;
}
inline bool PDBTopologyBuilder::readTopology(std::string _topologyFile) {
	pTopReader->reset();
	if (!pTopReader->open(_topologyFile)) {
		return false;
	}
	bool out = false;
	if (pTopReader->read()) {
		out = true;
	}
	pTopReader->close();
	fail_flag = !out;
	return out;
}
inline void PDBTopologyBuilder::reset() {
	polymerDefi.clear();
	atomMap.clear();
	pSystem = NULL;
	//vdwRescalingFactor = 1.0;
	//buildNonBondedInteractions = true;
	//elec14factor = 1;
	//dielectricConstant = 1;
	//useRdielectric = true;
	//useSolvation = false;
	//useSolvation = false;
	//solvent = pEEF1ParReader->getDefaultSolvent();
}

inline CharmmTopologyReader  * PDBTopologyBuilder::getCharmmTopologyReader(){return pTopReader; }
inline bool PDBTopologyBuilder::fail() const { return fail_flag;}

}

#endif
