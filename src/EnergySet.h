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

#ifndef ENERGYSET_H
#define ENERGYSET_H

#include <vector>
#include <string>
#include <iostream>

#include "Interaction.h"
#include "CharmmVdwInteraction.h"
#include "CharmmBondInteraction.h"
#include "CharmmElectrostaticInteraction.h"
#include "CharmmUreyBradleyInteraction.h"
#include "CharmmAngleInteraction.h"
#include "CharmmDihedralInteraction.h"
#include "CharmmImproperInteraction.h"
#include "UserDefinedInteraction.h"

#include "UserDefinedEnergy.h"


#include "Atom.h"
#include "AtomVector.h"

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


using namespace std;

class EnergySet {
	public:
		EnergySet();
		//EnergySet(const EnergySet & _set);
		~EnergySet();

		//enum InteractionTypes { CharmmVdw=0, CharmmElec=1, CharmmBond=2, CharmmAngle=3, CharmmUreyBradley=4, CharmmDihedral=5, CharmmImproper=6, CharmmEFF1=7 };

		void addInteraction(Interaction * _interaction);

		//unsigned int getTotalNumberOfInteractions(unsigned int _type);  // NOT DEFINED COMMENTED OUT
		unsigned int getTotalNumberOfInteractions(string _type);

		void deleteInteractionsWithAtom(Atom & _a);
		void deleteInteractionsWithAtoms(AtomVector & _atomVec);

		void setUseTerm(unsigned int _type);

		/* Calculate the energies */
		double calcEnergy(bool _activeOnly=true);
		double calcEnergy(string _selection, bool _activeOnly=true);
		double calcEnergy(string _selection1, string _selection2, bool _activeOnly=true);

		/********************************************************************
		 *
		 *  Terms can be set active or inactive (inactive won't be calcuated)
		 *  using their name or number:
		 *    
		 ********************************************************************/
		void setTermActive(string _termName, bool _active=true);
		void setAllTermsInactive();
		void setAllTermsActive();
		bool isTermActive(string _termName) const;
	
		/*************************************************
		 *   Get the data after a set of energies was
		 *   calculated (number of terms and energies by
		 *   term)
		 *************************************************/
		string getSummary () const;
		void printSummary() const {cout << getSummary();};
		double getTotalEnergy() const;
		double getTermEnergy(string _name) const;
		map<string, vector<Interaction*> > * getEnergyTerms();
		unsigned int getTotalNumberOfInteractionsCalculated() const;
		unsigned int getTermNumberOfInteractionsCalculated(string _name) const;


		/*
		  
		 */
		 vector<Interaction *> & getEnergyInteractions(Atom *a, Atom *b, string _termName);

		 class AtomPair : public pair<Atom *, Atom *> {
		        public:
		            AtomPair() : pair<Atom *, Atom *>(NULL,NULL){}
		            AtomPair(Atom *a, Atom *b) : pair<Atom *, Atom *>(a,b) {};

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


		 typedef map<AtomPair ,  vector<Interaction *>, cmpAtomPair > atomPairMap;
		 typedef map<AtomPair ,  vector<Interaction *>, cmpAtomPair >::iterator atomPairMapIt;


		void setCheckForCoordinates(bool _flag);
		bool getCheckForCoordinates() const;

	private:
		void deletePointers();
		void setup();
		//void copy(const EnergySet & _set);

		double calculateEnergy(string _selection1, string _selection2, bool _noSelect, bool _activeOnly=true);

		bool checkForCoordinates_flag;

		map<string, vector<Interaction*> > energyTerms;
		map<string, bool> activeEnergyTerms;
		map<string, unsigned int> interactionCounter;
		map<string, double> termTotal;
		double totalEnergy;
		unsigned int totalNumberOfInteractions;




		map<string, atomPairMap> pairInteractions;
		vector<Interaction *> blank; // Use as a return value in getEnegyInteractions, a hack I know..
		

		//map<Atom*, map<Atom*, AtomDistanceRelationship*> > atomDistanceRelationships;
		//map<Atom*, map<Atom*, map<Atom*, AtomAngleRelationship*> > > atomAngleRelationships;
		//map<Atom*, map<Atom*, map<Atom*, map<Atom*, AtomDihedralRelationship*> > > > atomDihedralRelationships;
		
		unsigned int stamp;


};

// INLINE 
inline map<string, vector<Interaction*> > * EnergySet::getEnergyTerms() {return & energyTerms;}
inline void EnergySet::setTermActive(string _termName, bool _active) {
	if (energyTerms.find(_termName) != energyTerms.end()) {
		activeEnergyTerms[_termName] = _active;
	//	cout << "UUU set " <<  _termName << " " << activeEnergyTerms[_termName] << endl;
	}
}
inline bool EnergySet::isTermActive(string _termName) const {
	map<string, bool>::const_iterator found = activeEnergyTerms.find(_termName);
	if (found != activeEnergyTerms.end()) {
		return found->second;
	}
	return false;
}
inline void EnergySet::setAllTermsInactive() {
	for (map<string, bool>::iterator k=activeEnergyTerms.begin(); k!=activeEnergyTerms.end(); k++) {
		k->second = false;
	}
}
inline void EnergySet::setAllTermsActive() {
	for (map<string, bool>::iterator k=activeEnergyTerms.begin(); k!=activeEnergyTerms.end(); k++) {
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

inline unsigned int EnergySet::getTotalNumberOfInteractions(string _type){
	map<string,vector<Interaction*> >::iterator it;
	it = energyTerms.find(_type);

	if (it == energyTerms.end()){
		//return -1; // BUG: returning -1 as usigned int will return 4294967295!!!
		return 0;
	}
	
	return (it->second).size();
}

#endif
