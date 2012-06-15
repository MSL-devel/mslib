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

#ifndef CHARMMENERGYCALCULATOR
#define CHARMMENERGYCALCULATOR

/*
  CharmmEnergyCalculator.  
      The purpose of this class is to compute energy between subsets of atoms, where the non-bonded interactions are
  computed "on-the-fly".  Thereby, the no "Interaction" objects are created for non-bonded energy calculations.  The class also
computes different subsets of energies: self,template and pair. Convient for protein design calculations. 

  Self     = interactions within a given position. 
  Template = interactions between a variable position and a fixed position.
  Pair     = interactions between a variable position and a variable position

*/

// STL Includes
#include <map>

// MSL Includes 
#include "System.h"
#include "AtomPointerVector.h"
#include "Interaction.h"
#include "CharmmParameterReader.h"
#include "RandomNumberGenerator.h"


namespace MSL { 
class CharmmEnergyCalculator {

	public:
		CharmmEnergyCalculator(std::string _charmmParameterFile);
		~CharmmEnergyCalculator(); 

		// Extract bonded terms from an EnergySet, otherwise we can query the Atom->isBoundTo(), Atom->isOneThree(), Atom->isOneFour() functions.
		bool extractBondedInteractions(EnergySet &_es);

		double calculateTotalEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2);

		// Components of Total energy
		double calculateSelfEnergy(System &_sys, int _position, int _rotamer);
		double calculateTemplateEnergy(System &_sys, int _position, int _rotamer, bool _calcAllForFixed=false);
		double calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2);

		// For quenching, etc., calculate the energy of a rotamer with the surrounding residues
		double calculateBackgroundEnergy(System &_sys, int _position1, int _rotamer1);
		double calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, std::vector< std::vector< std::vector< std::vector<double> > > > & rotamerInteractions, std::vector<uint> & currentAllRotamers);

		// Utility function. Sets must be unique or equal. If equal set _sameSet=true.
		std::map<std::string,double> calculatePairwiseEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet=false);
		std::map<std::string,double> calculatePairwiseNonBondedEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet=false);

		double getComputedEnergy(std::string _energyType);
		std::map<std::string,double> & getAllComputedEnergiesByType();
		void clearEnergiesByType();

		void setEnergyByType(bool _flag);
		bool getEnergyByType();

		void setEnergyByGroup(bool _flag);
		bool getEnergyByGroup();


		void   setVdwScale(double _vdwScale);
		double getVdwScale();
		// Set CharmmParameterReader
		CharmmParameterReader * getCharmmParameterReader();

		
		// Get interactions when in debug mode...
		std::map<Interaction*, int> & getInteractions();

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


		// Set NonBondedCutoffs On/Off for switching
		void setNonBondedCutoffs(double _fullyOn, double _fullyOff);
		double getNonBondedCutoffOn();
		double getNonBondedCutoffOff();


		bool getCollectNumberOfInteractions();
		void setCollectNumberOfInteractions(bool _flag);
		int getNumberInteractionsUsed(string _type);
		void resetNumberInteractionsUsed();



	private:
		// No calls to default constructor..we need charmm parameter file. 
		CharmmEnergyCalculator();

		std::map<std::string, double> energiesByType;
		std::map<std::string, std::map<std::string,double> > energiesByGrouping;
		std::map<Interaction*, int> interactionsComputed;
		std::map<std::string, int> numberOfInteractionsUsed;
		bool collectInteractionsFlag;

		bool storeEneByType;
		bool storeEneByGroup;

		CharmmParameterReader *parReader;

		double vdwRescalingFactor;
		
		double elec14factor;
		double dielectricConstant;
		bool useRdielectric;

		double nonBondCutoffOn;
		double nonBondCutoffOff;

		bool extractInteractions(EnergySet &_es, string type);
		pair<double,double> computeVdwElec(Atom *_a, Atom *_b, double _rmin, double _emin, double _Kq_q1_q2_rescal,int _stamp);

		std::vector<Interaction *> & getEnergyInteractions(Atom *a, Atom *b, std::string _termName);

		class AtomPair : public std::pair<Atom *, Atom *> {
		        public:
		            AtomPair() : std::pair<Atom *, Atom *>(NULL,NULL){}
		            AtomPair(Atom *a, Atom *b) : std::pair<Atom *, Atom *>(a,b) {};

		 };

		 struct cmpAtomPair {

			 bool operator()(const AtomPair &_apair, const AtomPair &_bpair) const{


				 bool val1 = (_apair.first == _bpair.first  && _apair.second == _bpair.second);
				 bool val2 = (_apair.first == _bpair.second && _apair.second == _bpair.first);

				 bool equal = val1 || val2;

				 if (equal){
					 return false;
				 } else if (_apair.first < _bpair.first){
					 return true;
				 } else if (_apair.first > _bpair.first){
					 return false;
				 } else {
					 return (_apair.second < _bpair.second);
				 }

			 }
		 };


		 typedef std::map<AtomPair ,  std::vector<Interaction *>, cmpAtomPair > atomPairMap;
		 typedef std::map<AtomPair ,  std::vector<Interaction *>, cmpAtomPair >::iterator atomPairMapIt;

		 std::map<std::string, atomPairMap> pairInteractions;
		 std::vector<Interaction *> blank; // Use as a return value in getEnegyInteractions, a hack I know..
		 void dumpAtomMap(string _term);
		 

};
inline std::map<Interaction*, int> & CharmmEnergyCalculator::getInteractions() { return interactionsComputed; }

inline double CharmmEnergyCalculator::calculateTotalEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2) {
	return (calculateTemplateEnergy(_sys, _position1, _rotamer1)  + 
         	calculateTemplateEnergy(_sys, _position2, _rotamer2)  + 
		calculateSelfEnergy(_sys, _position1, _rotamer1)  + 
		calculateSelfEnergy(_sys, _position2, _rotamer2)  + 
		calculatePairEnergy(_sys, _position1, _rotamer1, _position2,_rotamer2) );
}

inline void CharmmEnergyCalculator::setEnergyByType(bool _flag) { storeEneByType = true; }
inline bool CharmmEnergyCalculator::getEnergyByType() { return storeEneByType; }

inline void CharmmEnergyCalculator::setEnergyByGroup(bool _flag) { storeEneByGroup = true; }
inline bool CharmmEnergyCalculator::getEnergyByGroup() { return storeEneByGroup; }

inline double CharmmEnergyCalculator::getComputedEnergy(std::string _energyType){
	std::map<std::string,double>::iterator it;
	it = energiesByType.find(_energyType);
	if (it == energiesByType.end()){
		std::cerr << "WARNING NO ENERGY TYPE '"<<_energyType<<"' computed."<<std::endl;
		return 0;
	}

	return it->second;
}


inline std::map<std::string,double>& CharmmEnergyCalculator::getAllComputedEnergiesByType(){
	return energiesByType;
}

//inline void CharmmEnergyCalculator::setCharmmParameterReader(CharmmParameterReader *_par) { parReader = _par; }
inline CharmmParameterReader * CharmmEnergyCalculator::getCharmmParameterReader()         { return parReader; } 

//inline void   CharmmEnergyCalculator::setVdwScale(double _vdwScale) { vdwScale = _vdwScale; }
//inline double CharmmEnergyCalculator::getVdwScale() { return vdwScale; }
inline void  CharmmEnergyCalculator::clearEnergiesByType() {
	energiesByType.clear();
}
inline void CharmmEnergyCalculator::setVdwRescalingFactor(double _factor) {vdwRescalingFactor = _factor;}
inline double CharmmEnergyCalculator::setVdwRescalingFactor() const {return vdwRescalingFactor;}
inline void CharmmEnergyCalculator::setElec14factor(double _e14) {elec14factor = _e14;}
inline double CharmmEnergyCalculator::getElec14factor() const {return elec14factor;}
inline void CharmmEnergyCalculator::setDielectricConstant(double _diel) {dielectricConstant = _diel;}
inline double CharmmEnergyCalculator::getDielectricConstant() const {return dielectricConstant;}
inline void CharmmEnergyCalculator::setUseRdielectric(bool _flag) {useRdielectric = _flag;}
inline bool CharmmEnergyCalculator::getUseRdielectric() const {return useRdielectric;}
inline void CharmmEnergyCalculator::setNonBondedCutoffs(double _fullyOn, double _fullyOff) { 
	nonBondCutoffOn = _fullyOn;
	nonBondCutoffOff = _fullyOff;
}
inline double CharmmEnergyCalculator::getNonBondedCutoffOn() { return nonBondCutoffOn;}
inline double CharmmEnergyCalculator::getNonBondedCutoffOff(){ return nonBondCutoffOff;}

inline int CharmmEnergyCalculator::getNumberInteractionsUsed(string _type) {
  map<string,int>::iterator it = numberOfInteractionsUsed.find(_type); 
  if (it == numberOfInteractionsUsed.end()){
    return 0;
  } 
  return it->second;
}
inline void CharmmEnergyCalculator::resetNumberInteractionsUsed() { numberOfInteractionsUsed.clear();}
inline bool CharmmEnergyCalculator::getCollectNumberOfInteractions() { return collectInteractionsFlag; }
inline void CharmmEnergyCalculator::setCollectNumberOfInteractions(bool _flag) { collectInteractionsFlag = _flag;}
};

#endif
