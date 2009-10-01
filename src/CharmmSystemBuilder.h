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

#ifndef CHARMMSYSTEMBUILDER_H
#define CHARMMSYSTEMBUILDER_H

#include <iostream>
#include <vector>

#include "System.h"
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "CharmmParameterReader.h"
#include "PolymerSequence.h"

using namespace std;

class CharmmSystemBuilder {
	public:
		CharmmSystemBuilder();
		CharmmSystemBuilder(string _topologyFile, string _parameterFile);
		CharmmSystemBuilder(const CharmmSystemBuilder & _sysBuild);
		~CharmmSystemBuilder();

		void operator=(const CharmmSystemBuilder & _sysBuild);

		bool readTopology(string _topologyFile);
		bool readParameters(string _parameterFile);

		bool buildSystem(System & _pSys, const PolymerSequence & _sequence);
		
		bool getBuildNonBondedInteractions();
		void setBuildNonBondedInteractions(bool _flag);

		CharmmTopologyReader * getCharmmTopologyReader();
		CharmmParameterReader * getCharmmParameterReader();

	private:
		void setup();
		void copy(const CharmmSystemBuilder & _sysBuild);
		void deletePointers();
		void getAtomPointersFromMulti(string _name, vector<Atom*> & _out, vector<CharmmTopologyResidue*> & _position, vector<map<string, Atom*> > & _atomMap);
		vector<Atom*> getAtomPointers(string _name, vector<vector<vector<CharmmTopologyResidue*> > >::iterator & _chItr, vector<vector<CharmmTopologyResidue*> >::iterator & _posItr, vector<CharmmTopologyResidue*>::iterator & _idItr);

		vector<vector<vector<CharmmTopologyResidue*> > > polymerDefi;
		vector<vector<vector<map<string, Atom*> > > > atomMap;

		CharmmTopologyReader * pTopReader;
		CharmmParameterReader * pParReader;

		bool buildNonBondedInteractions;

};

inline bool CharmmSystemBuilder::readTopology(string _topologyFile) {pTopReader->reset(); if (!pTopReader->open(_topologyFile)) {return false;} return pTopReader->read();}
inline bool CharmmSystemBuilder::readParameters(string _parameterFile) {pParReader->reset(); if (!pParReader->open(_parameterFile)) {return false;} return pParReader->read();}

inline bool CharmmSystemBuilder::getBuildNonBondedInteractions()  { return buildNonBondedInteractions; }
inline void CharmmSystemBuilder::setBuildNonBondedInteractions(bool _flag) { buildNonBondedInteractions = _flag; }

inline CharmmTopologyReader  * CharmmSystemBuilder::getCharmmTopologyReader(){return pTopReader; }
inline CharmmParameterReader * CharmmSystemBuilder::getCharmmParameterReader(){return pParReader; }

#endif
