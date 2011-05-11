/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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

#include "EnergySet.h"

using namespace MSL;
using namespace std;


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
			//break;  // WHY WAS THERE A BREAK HERE???  IT WAS CAUSING A MEMORY LEAK!!!
		}
	}
	energyTerms.clear();


}

void EnergySet::resetTerm(string _term) {
	// remove all interactions for this term
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		if (k->first != _term) {
			continue;
		}
		for (vector<Interaction*>::iterator l=k->second.begin(); l!=k->second.end(); l++) {
			delete *l;
			//break;  // WHY WAS THERE A BREAK HERE???  IT WAS CAUSING A MEMORY LEAK!!!
		}
		k->second.clear();
		energyTerms.erase(k);
		break;
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
	energyTerms[_interaction->getName()].push_back(_interaction);
	if (activeEnergyTerms.find(_interaction->getName()) == activeEnergyTerms.end()) {
		activeEnergyTerms[_interaction->getName()] = true;
	}



	// Keep track of interactions by AtomPair key (Atom *, Atom *)
	vector<Atom *> &ats = _interaction->getAtomPointers();
	map<string,atomPairMap >::iterator pairIt;
	atomPairMapIt atIt;

	
	pairIt = pairInteractions.find(_interaction->getName());
	AtomPair ab(ats[0],ats[ats.size()-1]);		

	/*
	if (_interaction->getName() == "CHARMM_BOND"){
		cout << "ADDING INTERACTION: "<<ats[0]->getResidueNumber()<<" "<<ats[0]->getName()<<" "<<ats[ats.size()-1]->getResidueNumber()<<" "<<ats[ats.size()-1]->getName()<<endl;
	}
	*/

	if (pairIt == pairInteractions.end()){
		
		atomPairMap  tmp;
		tmp[ab].push_back(_interaction);
		
		pairInteractions[_interaction->getName()] = tmp;
	} else {

		atIt = (pairIt->second).find(ab);
		
		if (atIt == (pairIt->second).end()){
			
			(pairIt->second)[ab].push_back(_interaction);


		} else {


			// This output should not be suppressed, but it can be a lot of stuff to see.... debug/verbose flag?
			/*
			cerr << "WARNING 7724 EnergySet::addInteraction()... THIS ATOM PAIR ALREADY HAS THIS TYPE OF INTERACTION ("<<_interaction->getName()<<")"<<endl;

			for (uint a = 0; a < ats.size();a++){
				cerr << (*ats[a])<<endl;
			}

			// maybe print out the interaction found at 'atIt' ?
			cerr << "MATCHES WITH " <<(atIt->second)[0]->getName()<<endl;
			vector<Atom *> tmpAts = (atIt->second)[0]->getAtomPointers();
			for (uint a = 0; a < tmpAts.size();a++){
				cerr << (*tmpAts[a])<<endl;
			}


			cerr << "ADDING ANYWAYS, BUT MAKE SURE THIS IS NOT AN ERROR\n";

			cerr << " I HAVE SEEN THIS WITH PHE/TYR dihedrals: CG-CD2-CE2-CZ and CG-CD1-CE1-CZ\n";
			*/

			(pairIt->second)[ab].push_back(_interaction);

		}
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
			}
		}
		interactionCounter[k->first] = tmpTermCounter;
		termTotal[k->first] = tmpTermTotal;
		totalEnergy += tmpTermTotal;
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
		termTotal[k->first] = tmpTermTotal;
		totalEnergy += tmpTermTotal;
		totalNumberOfInteractions += k->second.size();
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

string EnergySet::getSummary() const{
	ostringstream os;	
	os << setiosflags(ios::left);
	os << "================  ======================  ===============" << endl;
	os << setw(20) <<"Interaction Type"<< setw(22) <<"Energy" << setw(15) << "Number of Terms" << endl;
	os << "================  ======================  ===============" << endl;
	for (map<string, double>::const_iterator l = termTotal.begin(); l!=termTotal.end(); l++) {
		double E = getTermEnergy(l->first);
		if (E<1E+14 && E>-1E+14) {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
		} else {
			os << resetiosflags(ios::right) << setw(20) << l->first << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getTermNumberOfInteractionsCalculated(l->first) << endl;
		}
	}
	os << "================  ======================  ===============" << endl;
	double E = getTotalEnergy();
	if (E<1E+14 && E>-1E+14) {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) <<setiosflags(ios::fixed)<< setprecision(6) << E << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
	} else {
		os << resetiosflags(ios::right) << setw(20) << "Total" << setw(20) << setiosflags(ios::right) << setiosflags(ios::fixed)<< setprecision(6) << "********************" << setw(15) << getTotalNumberOfInteractionsCalculated() << endl;
	}
	os << "================  ======================  ===============" << endl;
	return (os.str());

}

unsigned int EnergySet::getTermNumberOfInteractionsCalculated(string _name) const {
	map<string, uint>::const_iterator k = interactionCounter.find(_name);
	if (k != interactionCounter.end()) {
		return k->second;
	}
	return(0.0);


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
	return(-1.0);
	//return ((termTotal.find(_type))->second);
}


vector<Interaction *> & EnergySet::getEnergyInteractions(Atom *a, Atom *b, string _termName){

	// THIS FUNCTION COULD USE SOME COMMENTS TO EXPLAING WHAT IT IS DOING

	map<string, atomPairMap >::iterator typeIt;
	typeIt = pairInteractions.find(_termName);
	if (typeIt == pairInteractions.end()){
		return blank;
	}

	
	atomPairMapIt atomIt;
	
	AtomPair ab(a,b);

	atomIt = (typeIt->second).find(ab);
	if (atomIt == (typeIt->second).end()){

		// RRRRRRRR . I shouldn't have to look things up in reverse. The map<stl::pair, ...,cmpOperator> , yet cmpOperator doesn
		AtomPair ba(b,a);
		atomIt = (typeIt->second).find(ba);
		if (atomIt != (typeIt->second).end()){
			return atomIt->second;
		}
			

		return blank;
	}


	// We should have a set interactions between atoms a and b.
	return atomIt->second;
	
}


double EnergySet::calcEnergyAndEnergyGradient(vector<double> &_gradients){
	
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
				pair<double,vector<double> > partials  =  partialDerivative(ats);
				//cout << "Partials: "<<partials.first<<" "<<partials.second.size()<<endl;
				double e  = (*l)->getEnergy(partials.first,&partials.second);
				energy += e;

				for (uint a = 0; a < ats.size();a++){
					
					if (ats[a]->getMinimizationIndex() == -1){
						continue;
					}
					int fullGradIndex = 3*ats[a]->getMinimizationIndex();
					int localGradIndex = 3*(a+1);
					
					_gradients[fullGradIndex-3] += partials.second[localGradIndex-3];
					_gradients[fullGradIndex-2] += partials.second[localGradIndex-2];
					_gradients[fullGradIndex-1] += partials.second[localGradIndex-1];

				}
			}
		}
	}

	return energy;
}

void EnergySet::calcEnergyGradient(vector<double> &_gradients){


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

std::pair<double,std::vector<double> > EnergySet::partialDerivative(std::vector<Atom *> &ats){
        std::pair<double, std::vector<double> > partials;
	if (ats.size() == 2) {
		partials.first = CartesianGeometry::distanceDerivative(ats[0]->getCoor(), ats[1]->getCoor(),&(partials.second));
	}


	if (ats.size() == 3){
		partials.first = CartesianGeometry::angleDerivative(ats[0]->getCoor(),ats[1]->getCoor(),ats[2]->getCoor(),&(partials.second));
	}

	if (ats.size() == 4){
		partials.first = CartesianGeometry::dihedralDerivative(ats[0]->getCoor(),ats[1]->getCoor(),ats[2]->getCoor(),ats[3]->getCoor(),&(partials.second));
	}

	return partials;
}

void EnergySet::clearAllInteractions(){
  /* ***************
    this does not DELETE memory!!!!
    
  */
	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		k->second.clear();
	}
	energyTerms.clear();
}
