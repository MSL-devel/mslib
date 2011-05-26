/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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


#ifndef ENERGYSET_H
#define ENERGYSET_H

#include <vector>
#include <string>
#include <iostream>

#include "Interaction.h"
#include "SpringConstraintInteraction.h"
//#include "CharmmVdwInteraction.h"
//#include "CharmmBondInteraction.h"
//#include "CharmmElectrostaticInteraction.h"
//#include "CharmmUreyBradleyInteraction.h"
//#include "CharmmAngleInteraction.h"
//#include "CharmmDihedralInteraction.h"
//#include "CharmmImproperInteraction.h"
//#include "CharmmEEF1Interaction.h"
//#include "CharmmEEF1RefInteraction.h"



#include "Atom.h"
#include "AtomPointerVector.h"

#include <iomanip>

/*************************************************
 *  TO DO:
 *
 *  COPY CONSTRUCTOR
 *
 *  Implement and update function that calls the
 *  update function of all the interaction to 
 *  re-lookup for all parameters (and partial
 *  charges in the electrostaic term
 *************************************************/



namespace MSL { 
class EnergySet {
	public:
		EnergySet();
		//EnergySet(const EnergySet & _set);
		~EnergySet();

		//enum InteractionTypes { CharmmVdw=0, CharmmElec=1, CharmmBond=2, CharmmAngle=3, CharmmUreyBradley=4, CharmmDihedral=5, CharmmImproper=6, CharmmEFF1=7 };

		void addInteraction(Interaction * _interaction);

		//unsigned int getTotalNumberOfInteractions(unsigned int _type);  // NOT DEFINED COMMENTED OUT
		unsigned int getTotalNumberOfInteractions(std::string _type);

		// WARNING THE NEXT TWO FUNCTIONS ARE NOT IMPLEMENTED!!!!
		void deleteInteractionsWithAtom(Atom & _a);
		void deleteInteractionsWithAtoms(AtomPointerVector & _atomVec);
		void clearAllInteractions(); // this does not DELETE memory!!!! Avoids duplicate deletes when the same interaction is part of multiple EnergySets.

		void eraseTerm(std::string _term); // remove all interactions for this term - DELETES MEMORY

		void setUseTerm(unsigned int _type);

		/* Calculate the energies */
		/***********************************************************
		 *  20 Jan 2010 (AS) NOTE, CHANGED API TO FIX A BUG:
		 *  Before it was:
		 *  	double calcEnergy(bool _activeOnly=true);
		 *  	double calcEnergy(std::string _selection, bool _activeOnly=true);
		 *  The problem was that if calcEnergy("somename") was given, the std::string "somename"
		 *  was interpreted as a bool instead of a selection name, therefore the function
		 *  called would be calcEnergy(true) instead of calcEnergy("somename", true).
		 *  Fixed by adding a separate set of functions (calcEnergyAllAtoms) to calculate 
		 *  the energy of all active and inactive atoms.
		 ***********************************************************/
		double calcEnergy();
		double calcEnergy(std::string _selection);
		double calcEnergy(std::string _selection1, std::string _selection2);

		void calcEnergyGradient(std::vector<double> &_gradients);
		double calcEnergyAndEnergyGradient(std::vector<double> &_gradients);

		/* Calculate the energies including the interactions that inlcude atoms that belong to inactive side chains */
		double calcEnergyAllAtoms();
		double calcEnergyAllAtoms(std::string _selection);
		double calcEnergyAllAtoms(std::string _selection1, std::string _selection2);

		double calcEnergyOfSubset(std::string _subsetName, bool _activeOnly = true);

		void saveEnergySubset(std::string _subsetName);
		void saveEnergySubset(std::string _subsetName, std::string _selection);
		void saveEnergySubset(std::string _subsetName, std::string _selection1, std::string _selection2);
		void saveEnergySubsetAllAtoms(std::string _subsetName);
		void saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection);
		void saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection1, std::string _selection2);

		void removeEnergySubset(std::string _subsetName);

		/********************************************************************
		 *
		 *  Terms can be set active or inactive (inactive won't be calcuated)
		 *  using their name or number:
		 *    
		 ********************************************************************/
		void setTermActive(std::string _termName, bool _active=true);
		void setAllTermsInactive();
		void setAllTermsActive();
		bool isTermActive(std::string _termName) const;
	
		/*************************************************
		 *   Get the data after a set of energies was
		 *   calculated (number of terms and energies by
		 *   term)
		 *************************************************/
		std::string getSummary () const;
		void printSummary() const;
		double getTotalEnergy() const;
		double getTermEnergy(std::string _name) const;
		// the following returns a pointer to the map, the key is the energy term
		// (i.e. CHARMM_VDW) and the second is a vector of interaction pointers
		std::map<std::string, std::vector<Interaction*> > * getEnergyTerms();
		std::map<std::string, std::map<std::string, std::vector<Interaction*> > > * getEnergyTermsSubsets();
		unsigned int getTotalNumberOfInteractionsCalculated() const;
		unsigned int getTermNumberOfInteractionsCalculated(std::string _name) const;


		/*
		  
		 */
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


		void setCheckForCoordinates(bool _flag);
		bool getCheckForCoordinates() const;

	private:
		void deletePointers();
		void setup();
		//void copy(const EnergySet & _set);

		double calculateEnergy(std::string _selection1, std::string _selection2, bool _noSelect, bool _activeOnly);
		void saveEnergySubset(std::string _subsetName, std::string _selection1, std::string _selection2, bool _noSelect, bool _activeOnly);

		bool checkForCoordinates_flag;

		std::map<std::string, std::vector<Interaction*> > energyTerms;
		std::map<std::string, std::map<std::string, std::vector<Interaction*> > > energyTermsSubsets;
		std::map<std::string, bool> activeEnergyTerms;
		std::map<std::string, unsigned int> interactionCounter;
		std::map<std::string, double> termTotal;
		double totalEnergy;
		unsigned int totalNumberOfInteractions;




		std::map<std::string, atomPairMap> pairInteractions;
		std::vector<Interaction *> blank; // Use as a return value in getEnegyInteractions, a hack I know..
		
		unsigned int stamp;


};

// INLINE 
inline std::map<std::string, std::vector<Interaction*> > * EnergySet::getEnergyTerms() {return & energyTerms;}
inline std::map<std::string, std::map<std::string,std::vector<Interaction*> > > * EnergySet::getEnergyTermsSubsets() {return & energyTermsSubsets;}
inline void EnergySet::setTermActive(std::string _termName, bool _active) {
	if (energyTerms.find(_termName) != energyTerms.end()) {
		activeEnergyTerms[_termName] = _active;
	//	std::cout << "UUU set " <<  _termName << " " << activeEnergyTerms[_termName] << std::endl;
	}
}
inline bool EnergySet::isTermActive(std::string _termName) const {
	std::map<std::string, bool>::const_iterator found = activeEnergyTerms.find(_termName);
	if (found != activeEnergyTerms.end()) {
		return found->second;
	}
	return false;
}
inline void EnergySet::setAllTermsInactive() {
	for (std::map<std::string, bool>::iterator k=activeEnergyTerms.begin(); k!=activeEnergyTerms.end(); k++) {
		k->second = false;
	}
}
inline void EnergySet::setAllTermsActive() {
	for (std::map<std::string, bool>::iterator k=activeEnergyTerms.begin(); k!=activeEnergyTerms.end(); k++) {
		k->second = true;
	}
}
inline unsigned int EnergySet::getTotalNumberOfInteractionsCalculated() const {
	return(totalNumberOfInteractions);
}

inline double EnergySet::getTotalEnergy() const {
	return(totalEnergy);
}
inline void EnergySet::setCheckForCoordinates(bool _flag) {checkForCoordinates_flag = _flag;}
inline bool EnergySet::getCheckForCoordinates() const {return checkForCoordinates_flag;}

inline unsigned int EnergySet::getTotalNumberOfInteractions(std::string _type){
	std::map<std::string,std::vector<Interaction*> >::iterator it;
	it = energyTerms.find(_type);

	if (it == energyTerms.end()){
		//return -1; // BUG: returning -1 as usigned int will return 4294967295!!!
		return 0;
	}
	
	return (it->second).size();
}

inline void EnergySet::removeEnergySubset(std::string _subsetName) {
	std::map<std::string, std::map<std::string, std::vector<Interaction*> > >::iterator found = energyTermsSubsets.find(_subsetName);
	if (found != energyTermsSubsets.end()) {
		energyTermsSubsets.erase(found);
	}
}

inline void EnergySet::printSummary() const {std::cout << getSummary();};

}

#endif
