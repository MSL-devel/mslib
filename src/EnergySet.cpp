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

#include "EnergySet.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("EnergySet");

EnergySet::EnergySet() {
	setup();
}

/*
EnergySet::EnergySet(const EnergySet & _set) {
	setup();
	copy(_set);
}
*/

EnergySet::~EnergySet() {
	deletePointers();
}

void EnergySet::deletePointers() {
	/************************************
	 *  The way it is setup, the EnergySet receives pointers
	 *  and will take care of deleting them
	 *
	 *  Not ideal but it works
	 ************************************/
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		for (vector<Interaction*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			delete *l;
		}
	}
	energyTerms.clear();
	weights.clear();

	
}

void EnergySet::eraseTerm(string _term) {
	// remove all interactions for this term
	map<string,vector<Interaction*> >::iterator it = energyTerms.find(_term);
	if(it != energyTerms.end()) {
		for (vector<Interaction*>::iterator l=it->second.begin(); l!=it->second.end(); l++) {
			delete *l;
		}
		energyTerms.erase(it);
	}
	for (map<string, double>::iterator k=weights.begin(); k!=weights.end(); k++) {
		if (k->first == _term) {
			weights.erase(k);
			break;
		}
	}

}

void EnergySet::setup() {
	stamp = 0;
	totalEnergy = 0.0;
	checkForCoordinates_flag = false;
}


/*
void EnergySet::copy(const EnergySet & _set) {
	deletePointers();
	for (map<string, vector<Interaction*> >::const_iterator k=_set.energyTerms.begin(); k!=_set.energyTerms.end(); k++) {
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			switch ((*l)->getType()) {
				case 0:  // CHARMM VDW
					// curlies necessary because of the declaration within a switch statement
					{
						CharmmVdwInteraction * pInter =  static_cast<CharmmVdwInteraction*>(*l);
						energyTerms[k->first].push_back(new CharmmVdwInteraction(*pInter));
					}
					break;
				case 2:  // CHARMM BOND
					// curlies necessary because of the declaration within a switch statement
					{
					 	CharmmBondInteraction * pInter =  static_cast<CharmmBondInteraction*>(*l);
						energyTerms[k->first].push_back(new CharmmBondInteraction(*pInter));

					}
					break;
				default:
					cerr << "ERROR 34192: unknown interaction term (" << (*l)->getType() << ")" << endl;
					exit(34192);
					break;
			}
		}
	}
	stamp = _set.stamp;
	totalEnergy = _set.totalEnergy;
}
*/

void EnergySet::addInteraction(Interaction * _interaction) {
	// get the type from the pointer and push it back to the correct vector
	// in the map
	string name = _interaction->getName();
	energyTerms[name].push_back(_interaction);
	if (activeEnergyTerms.find(name) == activeEnergyTerms.end()) {
		activeEnergyTerms[name] = true;
	}
	if (weights.find(name) == weights.end()) {
		weights[name] = 1.0;
	}
}

/*   FUNCTIONS FOR ENERGY CALCULATION: 1) USE SELECTIONS   */

/* Public wrapper functions */
double EnergySet::calcEnergy() {
	return calculateEnergy("", "", true, true);
}

double EnergySet::calcEnergy(string _selection) {
	return calculateEnergy(_selection, _selection, false, true);
}

double EnergySet::calcEnergy(string _selection1, string _selection2) {
	return calculateEnergy(_selection1, _selection2, false, true);
}

double EnergySet::calcEnergyAllAtoms() {
	return calculateEnergy("", "", true, false);
}

double EnergySet::calcEnergyAllAtoms(string _selection) {
	return calculateEnergy(_selection, _selection, false, false);
}

double EnergySet::calcEnergyAllAtoms(string _selection1, string _selection2) {
	return calculateEnergy(_selection1, _selection2, false, false);
}

/* Private actual function */
double EnergySet::calculateEnergy(string _selection1, string _selection2, bool _noSelect, bool _activeOnly) {
	interactionCounter.clear();
	termTotal.clear();
	totalEnergy = 0.0;
	totalNumberOfInteractions = 0;
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		// for all the terms
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term, don't calculate it
			continue;
		}
		double tmpTermTotal = 0.0;
		unsigned int tmpTermCounter = 0;
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			// for all the interactions
			if ((!_activeOnly || (*l)->isActive()) && (_noSelect || (*l)->isSelected(_selection1, _selection2)) && (!checkForCoordinates_flag || (*l)->atomsHaveCoordinates())) {
				tmpTermCounter++;
				tmpTermTotal += (*l)->getEnergy(); 

				// Error checking code, add BOND Interactions to a hash, lookup later.
				//if ((*l)->getName() == "CHARMM_BOND"){
				//MSLOUT.getStaticLookup()[(*l)->toString()] = (*l)->getEnergy();
				//}
			}
		}
		interactionCounter[k->first] = tmpTermCounter;
		termTotal[k->first] = tmpTermTotal * weights[k->first];
		totalEnergy += termTotal[k->first];
		totalNumberOfInteractions += interactionCounter[k->first];
	}
	return totalEnergy;

}
double EnergySet::calcEnergyWithoutSwitchingFunction() {
	// only active terms
	interactionCounter.clear();
	termTotal.clear();
	totalEnergy = 0.0;
	totalNumberOfInteractions = 0;

	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		// for all the terms
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term, don't calculate it
			continue;
		}
		double tmpTermTotal = 0.0;
		unsigned int tmpTermCounter = 0;
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			// for all the interactions
			if ((*l)->isActive() && (!checkForCoordinates_flag || (*l)->atomsHaveCoordinates())) {
				tmpTermCounter++;
				// energy without applying the switching function 
				tmpTermTotal += (*l)->getEnergy((vector<double>*)NULL); 
				 
			}
		}
		interactionCounter[k->first] = tmpTermCounter;
		termTotal[k->first] = tmpTermTotal * weights[k->first];
		totalEnergy += termTotal[k->first];
		totalNumberOfInteractions += tmpTermCounter;
	}
	return totalEnergy;

}

/*   FUNCTIONS FOR ENERGY CALCULATION: 2) SAVE SUBSETS AND CALCULATE THEM FOR BETTER PERFORMANCE ON REPEATED CALCULATIONS ON THE SAME SELECTIONS   */

/* function to calculate the energy of a subset */
double EnergySet::calcEnergyOfSubset(string _subsetName, bool _activeOnly) {
	interactionCounter.clear();
	termTotal.clear();
	totalEnergy = 0.0;
	totalNumberOfInteractions = 0;

	map<string, map<string, vector<Interaction*> > >::iterator found = energyTermsSubsets.find(_subsetName);
	if (found == energyTermsSubsets.end()) {
		// the subset not found
		return 0.0;
	}

	for (map<string, vector<Interaction*> >::iterator k=found->second.begin(); k!=found->second.end(); k++) {
		// for all the terms
		double tmpTermTotal = 0.0;
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			// for all the interactions
			if(!_activeOnly || (*l)->isActive()) {
				tmpTermTotal += (*l)->getEnergy(); 
			}
		}
		interactionCounter[k->first] = k->second.size();
		termTotal[k->first] = tmpTermTotal * weights[k->first];
		totalEnergy += termTotal[k->first];
		totalNumberOfInteractions += interactionCounter[k->first];
	}
	return totalEnergy;
}

/* Public wrapper functions for saving subsets */
void EnergySet::saveEnergySubset(string _subsetName) {
	return saveEnergySubset(_subsetName, "", "", true, true);
}

void EnergySet::saveEnergySubset(string _subsetName, string _selection) {
	return saveEnergySubset(_subsetName, _selection, _selection, false, true);
}

void EnergySet::saveEnergySubset(string _subsetName, string _selection1, string _selection2) {
	return saveEnergySubset(_subsetName, _selection1, _selection2, false, true);
}

void EnergySet::saveEnergySubsetAllAtoms(string _subsetName) {
	return saveEnergySubset(_subsetName, "", "", true, false);
}

void EnergySet::saveEnergySubsetAllAtoms(string _subsetName, string _selection) {
	return saveEnergySubset(_subsetName, _selection, _selection, false, false);
}

void EnergySet::saveEnergySubsetAllAtoms(string _subsetName, string _selection1, string _selection2) {
	return saveEnergySubset(_subsetName, _selection1, _selection2, false, false);
}

/* Private actual function for saving subsets */
void EnergySet::saveEnergySubset(string _subsetName, string _selection1, string _selection2, bool _noSelect, bool _activeOnly) {

	energyTermsSubsets[_subsetName].clear(); // reset the subset if existing
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		// for all the terms
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term, don't calculate it
			continue;
		}
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			// for all the interactions
			if ((!_activeOnly || (*l)->isActive()) && (_noSelect || (*l)->isSelected(_selection1, _selection2)) && (!checkForCoordinates_flag || (*l)->atomsHaveCoordinates())) {
				// add the interaction to the subset
				energyTermsSubsets[_subsetName][k->first].push_back(*l);
			}
		}
	}
}

/*   FUNCTIONS FOR ENERGY CALCULATION: DONE   */

string EnergySet::getSummary(unsigned int _precision) const{
	if (_precision > 15) {
		_precision = 15;
	}
	ostringstream os;	
	os << setiosflags(ios::left);
	os << "================  ======================  ===============" << endl;
	os << setw(18) <<"Interaction Type"<< setw(24) <<"Energy" << setw(15) << "Interactions" << endl;
	os << "================  ======================  ===============" << endl;
	for (map<string, double>::const_iterator l = termTotal.begin(); l!=termTotal.end(); l++) {
		double E = getTermEnergy(l->first);
		if (E<1E+14 && E>-1E+14) {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(_precision) << E << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
//			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
//			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(15) << E << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
		} else {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
//			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
//			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(15) << "********************" << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
		}
	}
	os << "================  ======================  ===============" << endl;
	double E = getTotalEnergy();
	if (E<1E+14 && E>-1E+14) {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(_precision) << E << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
//		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
//		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(15) << E << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
	} else {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(_precision) << "********************" << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
//		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
//		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(15) << "********************" << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
	}
	os << "================  ======================  ===============" << endl;
	return (os.str());

}

unsigned int EnergySet::getTermNumberOfInteractionsCalculated(string _name) const {
	map<string, uint>::const_iterator k = interactionCounter.find(_name);
	if (k != interactionCounter.end()) {
		return k->second;
	}
	return(0);


	/*  Unsafe but fast way
	 * 
	 * return ((interactionCounter.find(_type))->second);
	 */
}

double EnergySet::getTermEnergy(string _name) const {
	map<string, double>::const_iterator k = termTotal.find(_name);
	if (k != termTotal.end()) {
		return k->second;
	}
	return(0.0);
	//return ((termTotal.find(_type))->second);
}



double EnergySet::calcEnergyAndEnergyGradient(vector<double> &_gradients){
	
	// TODO: should we use term weights in minimization?
	// For each interaction
	double energy = 0.0;
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {


		// Only compute active energy terms
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term
			continue;
		}

		// Loop over each interaction
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {


			// If interaction is active...
			if ((*l)->isActive()){  

				vector<Atom *> &ats = (*l)->getAtomPointers();
				/*
				cout << "Calculating gradient of "<<(*l)->getName()<<endl;
				for (uint a = 0; a < ats.size();a++){
					cout << "\t"<<ats[a]->toString()<<endl;
				}
				*/
				// TODO: combine partialDerivative and getEnergy with the gradients into a single function
				//cout << "Partials: "<<partials.first<<" "<<partials.second.size()<<endl;
				vector<double> gradient;
				double e  = (*l)->getEnergy(&gradient);
				energy += e;

				for (uint a = 0; a < ats.size();a++){
					
					if (ats[a]->getMinimizationIndex() == -1){
						continue;
					}
					int fullGradIndex = 3*ats[a]->getMinimizationIndex();
					int localGradIndex = 3*(a+1);
					
					_gradients[fullGradIndex-3] += gradient[localGradIndex-3];
					_gradients[fullGradIndex-2] += gradient[localGradIndex-2];
					_gradients[fullGradIndex-1] += gradient[localGradIndex-1];

				}
			}
		}
	}

	return energy;
}

void EnergySet::calcEnergyGradient(vector<double> &_gradients){

	// TODO: should we use term weights in minimization?
	// For each interaction

	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {


		// Only compute active energy terms
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term
			continue;
		}

		// Loop over each interaction
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {


			// If interaction is active...
			if ((*l)->isActive()){  

				vector<Atom *> &ats = (*l)->getAtomPointers();
				vector<double> gradient = (*l)->getEnergyGrad();		   


				for (uint a = 0; a < ats.size();a++){
					
					if (ats[a]->getMinimizationIndex() == -1){
						continue;
					}
					int fullGradIndex = 3*ats[a]->getMinimizationIndex();
					int localGradIndex = 3*(a+1);
					
					_gradients[fullGradIndex-3] += gradient[localGradIndex-3];
					_gradients[fullGradIndex-2] += gradient[localGradIndex-2];
					_gradients[fullGradIndex-1] += gradient[localGradIndex-1];
				}
			}
		}
	}



	
}


void EnergySet::clearAllInteractions(){
  /* ***************
    this does not DELETE memory!!!!
    
  */
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		k->second.clear();
	}
	energyTerms.clear();
	weights.clear();
}

void EnergySet::deleteInteractionsWithAtom(Atom & _a, string _type) {
	AtomPointerVector av(1, &_a);
	deleteInteractionsWithAtoms(av,_type);
}
void EnergySet::deleteInteractionsWithAtoms(AtomPointerVector & _atomVec, string _type) {
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		if(_type != "" && _type != k->first) {
			continue;
		}
		for (vector<Interaction*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			for (AtomPointerVector::iterator m=_atomVec.begin(); m!=_atomVec.end(); m++) {
				if ((*l)->hasAtom(*m)) {
				//	vector<Atom*> & atoms = (*l)->getAtomPointers();
				//	cerr << "UUUY11 " << k->first << " " << *l << " ";
				//	for (vector<Atom*>::iterator m=atoms.begin(); m!=atoms.end(); m++) {
				//		cerr << **m << "/";
				//	}
				//	cerr << endl;
				//	cout << "UUUU Deleting " << (*l)->toString() << endl;
					delete *l;
					k->second.erase(l);
					l--;
					break;
				}
			}
		}
	}
	energyTermsSubsets.clear();
}

