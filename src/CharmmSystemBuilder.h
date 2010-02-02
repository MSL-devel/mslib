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

		bool buildSystem(System & _system, const PolymerSequence & _sequence);
		bool updateNonBonded(System & _system, double _ctonnb=0.0, double _ctofnb=0.0, double _cutnb=0.0);
		
		bool getBuildNonBondedInteractions();
		void setBuildNonBondedInteractions(bool _flag);

		CharmmTopologyReader * getCharmmTopologyReader();
		CharmmParameterReader * getCharmmParameterReader();

		void setVdwRescalingFactor(double _factor);
		double setVdwRescalingFactor() const;

		// rescaling of the 1-4 electostatic interactions.  should be 1 for charmm 22
		// and 0.6 for charmm 19
		void setElec14factor(double _e14);
		double getElec14factor() const;

		// the dielectric constant
		void setDielectricConstant(double _diel);
		double getDielectricConstant() const;

		// use a distance dependent dielectric
		void setUseRdielectric(bool _flag);
		bool getUseRdielectric() const;

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

		double vdwRescalingFactor;
		
		double elec14factor;
		double dielectricConstant;
		bool useRdielectric;

};

inline bool CharmmSystemBuilder::readTopology(string _topologyFile) {pTopReader->reset(); if (!pTopReader->open(_topologyFile)) {return false;} return pTopReader->read();}
inline bool CharmmSystemBuilder::readParameters(string _parameterFile) {pParReader->reset(); if (!pParReader->open(_parameterFile)) {return false;} return pParReader->read();}

inline bool CharmmSystemBuilder::getBuildNonBondedInteractions()  { return buildNonBondedInteractions; }
inline void CharmmSystemBuilder::setBuildNonBondedInteractions(bool _flag) { buildNonBondedInteractions = _flag; }

inline CharmmTopologyReader  * CharmmSystemBuilder::getCharmmTopologyReader(){return pTopReader; }
inline CharmmParameterReader * CharmmSystemBuilder::getCharmmParameterReader(){return pParReader; }
inline void CharmmSystemBuilder::setVdwRescalingFactor(double _factor) {vdwRescalingFactor = _factor;}
inline double CharmmSystemBuilder::setVdwRescalingFactor() const {return vdwRescalingFactor;}
inline void CharmmSystemBuilder::setElec14factor(double _e14) {elec14factor = _e14;}
inline double CharmmSystemBuilder::getElec14factor() const {return elec14factor;}
inline void CharmmSystemBuilder::setDielectricConstant(double _diel) {dielectricConstant = _diel;}
inline double CharmmSystemBuilder::getDielectricConstant() const {return dielectricConstant;}
inline void CharmmSystemBuilder::setUseRdielectric(bool _flag) {useRdielectric = _flag;}
inline bool CharmmSystemBuilder::getUseRdielectric() const {return useRdielectric;}

#endif
