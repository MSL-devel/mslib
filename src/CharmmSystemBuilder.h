/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
Kulp DW et al. "Structural informatics, modeling and design with a open 
source Molecular Software Library (MSL)" (2012) J. Comp. Chem, in press
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


#ifndef CHARMMSYSTEMBUILDER_H
#define CHARMMSYSTEMBUILDER_H

#include <iostream>
#include <vector>

#include "System.h"
#include "PDBReader.h"
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "CharmmParameterReader.h"
#include "CharmmEEF1ParameterReader.h"
#include "PolymerSequence.h"
#include "CharmmVdwInteraction.h"
#include "CharmmBondInteraction.h"
#include "CharmmElectrostaticInteraction.h"
#include "CharmmUreyBradleyInteraction.h"
#include "CharmmAngleInteraction.h"
#include "CharmmDihedralInteraction.h"
#include "CharmmImproperInteraction.h"
#include "CharmmEEF1Interaction.h"
#include "CharmmEEF1RefInteraction.h"
#include "RandomNumberGenerator.h"


namespace MSL { 
class CharmmSystemBuilder {
	public:
		CharmmSystemBuilder();
		CharmmSystemBuilder(System & _system, std::string _topologyFile, std::string _parameterFile, std::string _solvationFile="");
	//	CharmmSystemBuilder(std::string _topologyFile, std::string _parameterFile, std::string _solvationFile=""); // DEPRECATED, pass System
		CharmmSystemBuilder(const CharmmSystemBuilder & _sysBuild);
		~CharmmSystemBuilder();

		void operator=(const CharmmSystemBuilder & _sysBuild);

		void setSystem(System & _system);

		bool readTopology(std::string _topologyFile);
		bool readParameters(std::string _parameterFile);
		bool readSolvation(std::string _solvationFile);
		void setSolvent(std::string _solvent);

		bool buildSystem(const PolymerSequence & _sequence);
		bool buildSystemFromPDB(std::string _fileName); // build from a PDB in CHARMM name format
	//	bool buildSystem(System & _system, const PolymerSequence & _sequence); // DEPRECATED, system in constructor
	//	bool updateNonBonded(System & _system, double _ctonnb=0.0, double _ctofnb=0.0, double _cutnb=0.0); // DEPRECATED!!!!
		bool updateNonBonded(double _ctonnb=0.0, double _ctofnb=0.0, double _cutnb=0.0, bool _ignoreNonVariable = false);

		/**************************************************
		 * Add one or more new identities to a position,
		 * if a series of backbone atoms are specified
		 * (_bbAtoms) then the coordinates are passed to
		 * the atoms with the same name in the new residues
		 * to preserve the position of the backbone
		 **************************************************/
		bool addIdentity(std::string _positionId, std::string _resName, std::string _bbAtoms="N CA C O HN"); // id "A,37"
		bool addIdentity(std::string _positionId, const std::vector<std::string> & _resNames, std::string _bbAtoms="N CA C O HN");
		bool addIdentity(Position & _pos, std::string _resName, std::string _bbAtoms="N CA C O HN");
		bool addIdentity(Position & _pos, const std::vector<std::string> & _resNames, std::string _bbAtoms="N CA C O HN");

		/**************************************************
		 *  MUTATION FUNCTIONS
		 * 
		 *  TODO: right now the new identity is placed at
		 *        the end of the identity vector in the Position.
		 *        It should take the place of the old identity
		 **************************************************/
		// replace the current identity
		/*
		bool mutate(std::string _positionId, std::string _oldResName, _newResname, std::string _bbAtoms="N CA C O HN"); // id "A,37"
		bool mutate(Position & _pos, , std::string _oldResName, _newResname, std::string _bbAtoms="N CA C O HN");
		// replace a specific identity
		bool mutate(std::string _positionId, std::string _oldResName, _newResname, std::string _bbAtoms="N CA C O HN"); // id "A,37"
		bool mutate(Position & _pos, , std::string _oldResName, _newResname, std::string _bbAtoms="N CA C O HN");
		*/
		bool removeIdentity(std::string _positionId, string _resName);
		bool removeIdentity(Position & _pos, string _resName);


		
		bool getBuildNonBondedInteractions();
		void setBuildNonBondedInteractions(bool _flag);

		CharmmTopologyReader * getCharmmTopologyReader();
		CharmmParameterReader * getCharmmParameterReader();
		CharmmEEF1ParameterReader * getCharmmEEF1ParameterReader();

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

		// use group-based non-bond cutoffs
		void setUseGroupCutoffs(bool _flag);
		bool getUseGroupCutoffs() const;

		bool fail() const; // return false if reading toppar failed

		// introduce the ability to disable the creation of the pairwise lookup table
		// pairInteractions, used only with on-the-fly objects
		void setCreatePairwiseTable(bool _flag); 

	private:
		void setup();
		void copy(const CharmmSystemBuilder & _sysBuild);
		void deletePointers();
		void reset();
		void getAtomPointersFromMulti(std::string _name, std::vector<Atom*> & _out, std::vector<CharmmTopologyResidue*> & _position, std::vector<std::map<std::string, Atom*> > & _atomMap);
		std::vector<Atom*> getAtomPointers(std::string _name, std::vector<std::vector<std::vector<CharmmTopologyResidue*> > >::iterator & _chItr, std::vector<std::vector<CharmmTopologyResidue*> >::iterator & _posItr, std::vector<CharmmTopologyResidue*>::iterator & _idItr);

		std::vector<std::vector<std::vector<CharmmTopologyResidue*> > > polymerDefi;
		std::vector<std::vector<std::vector<std::map<std::string, Atom*> > > > atomMap;

		System * pSystem;

		CharmmTopologyReader * pTopReader;
		CharmmParameterReader * pParReader;
		CharmmEEF1ParameterReader * pEEF1ParReader;

		bool buildNonBondedInteractions;

		double vdwRescalingFactor;
		
		double elec14factor;
		double dielectricConstant;
		bool useRdielectric;
		bool useSolvation;
		bool useGroupCutoffs;
		std::string solvent;

		bool fail_flag;

		// if false the pairwise atom intearction tables are not created in the EnergySet
		// this is a temporary hack for speeding up things when on-the-fly is not used
		bool createPairInteractions_flag; 

};
inline void CharmmSystemBuilder::setSystem(System & _system) {
	reset();
	pSystem = &_system;
}
inline bool CharmmSystemBuilder::readTopology(std::string _topologyFile) {
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
inline bool CharmmSystemBuilder::readParameters(std::string _parameterFile) {
	pParReader->reset();
	if (!pParReader->open(_parameterFile)) {
		return false;
	} 
	bool out = false;
	if (pParReader->read()) {
		out = true;
	}
	pParReader->close();
	fail_flag = !out;
	return out;
}
inline bool CharmmSystemBuilder::readSolvation(std::string _solvationFile) {
	useSolvation = true;
	pEEF1ParReader->reset();
	if (!pEEF1ParReader->open(_solvationFile)) {
		return false;
	} 
	bool out = false;
	if (pEEF1ParReader->read()) {
		out = true;
	}
	pEEF1ParReader->close();
	useRdielectric = true;

	fail_flag = !out;
	return out;
}
inline void CharmmSystemBuilder::reset() {
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

inline bool CharmmSystemBuilder::getBuildNonBondedInteractions()  { return buildNonBondedInteractions; }
inline void CharmmSystemBuilder::setBuildNonBondedInteractions(bool _flag) { buildNonBondedInteractions = _flag; }

inline CharmmTopologyReader  * CharmmSystemBuilder::getCharmmTopologyReader(){return pTopReader; }
inline CharmmParameterReader * CharmmSystemBuilder::getCharmmParameterReader(){return pParReader; }
inline CharmmEEF1ParameterReader * CharmmSystemBuilder::getCharmmEEF1ParameterReader(){return pEEF1ParReader; }
inline void CharmmSystemBuilder::setVdwRescalingFactor(double _factor) {vdwRescalingFactor = _factor;}
inline double CharmmSystemBuilder::setVdwRescalingFactor() const {return vdwRescalingFactor;}
inline void CharmmSystemBuilder::setElec14factor(double _e14) {elec14factor = _e14;}
inline double CharmmSystemBuilder::getElec14factor() const {return elec14factor;}
inline void CharmmSystemBuilder::setDielectricConstant(double _diel) {dielectricConstant = _diel;}
inline double CharmmSystemBuilder::getDielectricConstant() const {return dielectricConstant;}
inline void CharmmSystemBuilder::setUseRdielectric(bool _flag) {useRdielectric = _flag;}
inline bool CharmmSystemBuilder::getUseRdielectric() const {return useRdielectric;}
inline void CharmmSystemBuilder::setUseGroupCutoffs(bool _flag) { useGroupCutoffs = _flag;}
inline bool CharmmSystemBuilder::getUseGroupCutoffs() const { return useGroupCutoffs; }
inline bool CharmmSystemBuilder::fail() const { return fail_flag;}
inline void CharmmSystemBuilder::setSolvent(std::string _solvent) {solvent = _solvent;}
inline void CharmmSystemBuilder::setCreatePairwiseTable(bool _flag) { createPairInteractions_flag = _flag; }

}

#endif
