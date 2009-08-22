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

#ifndef ATOMICPAIRWISEENERGY
#define ATOMICPAIRWISEENERGY

/*
  This class takes two conformations and computes energies.  This class contains the
  logic to compute proper self,template and pair atomic energies.

  (including 1-2,1-3,1-4 exclusions)

*/

// STL Includes
#include <map>

// MSL Includes 
#include "System.h"
#include "AtomVector.h"
#include "Interaction.h"
#include "CharmmParameterReader.h"


class AtomicPairwiseEnergy {

	public:
		AtomicPairwiseEnergy(string _charmmParameterFile);
		~AtomicPairwiseEnergy(); 
		
		
		double calculateTotalEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2);

		// Components of Total energy
		double calculateSelfEnergy(System &_sys, int _position1, int _rotamer1);
		double calculateTemplateEnergy(System &_sys, int _position1, int _rotamer1);
		double calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2);

		// For quenching, etc., calculate the energy of a rotamer with the surrounding residues
		double calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, vector< vector< vector< vector<double> > > > & rotamerInteractions, vector<uint> & currentAllRotamers);

		// Utility function. Sets must be unique or equal. If equal set _sameSet=true.
		map<string,double> calculatePairwiseEnergy(System &_sys, AtomVector &_a, AtomVector &_b, bool _sameSet=false);

		map<string,double> calculatePairwiseNonBondedEnergy(System &_sys, AtomVector &_a, AtomVector &_b, bool _sameSet=false);

		double getComputedEnergy(string _energyType);
		map<string,double> & getAllComputedEnergiesByType();


		void setEnergyByType(bool _flag);
		bool getEnergyByType();

		void setEnergyByGroup(bool _flag);
		bool getEnergyByGroup();


		// Set CharmmParameterReader
		CharmmParameterReader * getCharmmParameterReader();

		// Get interactions when in debug mode...
		map<Interaction*, int> & getInteractions();

	private:
		// No calls to default constructor..we need charmm parameter file. 
		AtomicPairwiseEnergy();

		map<string, double> energiesByType;
		map<string, map<string,double> > energiesByGrouping;
		map<Interaction*, int> interactionsComputed;
		bool storeEneByType;
		bool storeEneByGroup;

		CharmmParameterReader *parReader;

};
inline map<Interaction*, int> & AtomicPairwiseEnergy::getInteractions() { return interactionsComputed; }

inline double AtomicPairwiseEnergy::calculateTotalEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2) {
	return (calculateTemplateEnergy(_sys, _position1, _rotamer1)  + 
         	calculateTemplateEnergy(_sys, _position2, _rotamer2)  + 
		calculateSelfEnergy(_sys, _position1, _rotamer1)  + 
		calculateSelfEnergy(_sys, _position2, _rotamer2)  + 
		calculatePairEnergy(_sys, _position1, _rotamer1, _position2,_rotamer2) );
}

inline void AtomicPairwiseEnergy::setEnergyByType(bool _flag) { storeEneByType = true; }
inline bool AtomicPairwiseEnergy::getEnergyByType() { return storeEneByType; }

inline void AtomicPairwiseEnergy::setEnergyByGroup(bool _flag) { storeEneByGroup = true; }
inline bool AtomicPairwiseEnergy::getEnergyByGroup() { return storeEneByGroup; }

inline double AtomicPairwiseEnergy::getComputedEnergy(string _energyType){
	map<string,double>::iterator it;
	it = energiesByType.find(_energyType);
	if (it == energiesByType.end()){
		cerr << "WARNING NO ENERGY TYPE '"<<_energyType<<"' computed."<<endl;
		return 0;
	}

	return it->second;
}


inline map<string,double>& AtomicPairwiseEnergy::getAllComputedEnergiesByType(){
	return energiesByType;
}

//inline void AtomicPairwiseEnergy::setCharmmParameterReader(CharmmParameterReader *_par) { parReader = _par; }
inline CharmmParameterReader * AtomicPairwiseEnergy::getCharmmParameterReader()         { return parReader; } 
#endif
