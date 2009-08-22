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

#include "EnergySet.h"

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
			break;
		}
	}
	energyTerms.clear();


}

void EnergySet::setup() {
	stamp = 0;
	totalEnergy = 0.0;
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
	vector<Atom *> &ats = _interaction->getAtoms();
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
			vector<Atom *> tmpAts = (atIt->second)[0]->getAtoms();
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

double EnergySet::calcEnergy(bool _activeOnly) {
	return calculateEnergy("", "", true, _activeOnly);
}

double EnergySet::calcEnergy(string _selection, bool _activeOnly) {
	return calculateEnergy(_selection, _selection, false, _activeOnly);
}

double EnergySet::calcEnergy(string _selection1, string _selection2, bool _activeOnly) {
	return calculateEnergy(_selection1, _selection2, false, _activeOnly);
}

double EnergySet::calculateEnergy(string _selection1, string _selection2, bool _noSelect, bool _activeOnly) {
	interactionCounter.clear();
	termTotal.clear();
	totalEnergy = 0.0;
	totalNumberOfInteractions = 0;


	for (map<string, vector<Interaction*> >::iterator k=energyTerms.begin(); k!=energyTerms.end(); k++) {
		if (activeEnergyTerms.find(k->first) == activeEnergyTerms.end() || !activeEnergyTerms[k->first]) {
			// inactive term
			continue;
		}
		double tmpTermTotal = 0.0;
		unsigned int tmpTermCounter = 0;
		for (vector<Interaction*>::const_iterator l=k->second.begin(); l!=k->second.end(); l++) {
			if ((!_activeOnly || (*l)->isActive()) && (_noSelect || (*l)->isSelected(_selection1, _selection2))) {
				tmpTermCounter++;
				tmpTermTotal += (*l)->getEnergy(); 
			//	cout << "UUU Interaction " << (*l)->toString() << endl;
				//interactionCounter[k->first]++;
				//termTotal[k->first] += (*l)->getEnergy(stamp);
				//totalNumberOfInteractions++;
				//cout << "Energy " << (*l)->getEnergy(stamp) << endl;
			}
		}
		interactionCounter[k->first] = tmpTermCounter;
		termTotal[k->first] = tmpTermTotal;
		totalEnergy += tmpTermTotal;
		totalNumberOfInteractions += tmpTermCounter;
	}
	//stamp++;
//	return 0.0;
	return totalEnergy;

}

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

