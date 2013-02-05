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


#include "CharmmEnergyCalculator.h"
#include "CharmmEnergy.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("CharmmEnergyCalculator");

CharmmEnergyCalculator::CharmmEnergyCalculator(string _charmmParameterFile){
	storeEneByType  = false;
	storeEneByGroup = false;

	parReader = new CharmmParameterReader();
	parReader->open(_charmmParameterFile);
	parReader->read();
	parReader->close();
	parReader->createVdwParamPairs();

	vdwRescalingFactor = 1.0;
		
	elec14factor = 1.0;
	dielectricConstant = 1.0;
	useRdielectric = true;

	// Setting these cutoff values to 0.0 means NO cutoffs, compute ALL pairwise interactions regardless of distance.
	nonBondCutoffOn  = 0.0;
	nonBondCutoffOff = 0.0;

	collectInteractionsFlag = false;
}

CharmmEnergyCalculator::~CharmmEnergyCalculator(){
	delete(parReader);

	for (map<string,atomPairMap>::iterator k = pairInteractions.begin(); k != pairInteractions.end();k++){
	    pairInteractions.erase(k);
	}

}

/*
  Calculate all pairwise energy between atoms within a residue
 */
double CharmmEnergyCalculator::calculateSelfEnergy(System &_sys, int _position, int _rotamer){



	//cout << "Calculating self for position "<<_position<<" and rotamer "<<_rotamer<<endl;

	// Dirty.. I set , but never reset.  There is no such function yet.
	//      we don't keep track of currentConformations in Residue...
	//      basically if we add a "currentConformation" vector in Residue and
	//      clear+then fill it each time we call setActiveRotamer we will be
	//      ok.... maybe at position level?

	// Set the active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);



	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position).getAtomPointers();


	//cout << "Atoms for "<<_position<<" "<<_rotamer<<" are: "<<atoms1.toString()<<endl;
	double energy = 0.0;

	// Calculate atom-atom energies of residue
	map<string,double> energies = calculatePairwiseEnergy(atoms1,atoms1,true);

	// Add energy by type to member variable "energiesByType"
	if (storeEneByType) {
		map<string,double>::iterator it;
		for (it = energies.begin();it != energies.end();it++){
			energiesByType[it->first] += it->second;
		}
	}

	energy = energies["TOTAL"];


	return energy;
	
}

/*
  Calculate energy of rotamer with positions that are fixed.  Much like calculateTemplateEnergy, but use only non-bonded energies
 */
double CharmmEnergyCalculator::calculateBackgroundEnergy(System &_sys, int _position, int _rotamer){
	//Set active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position).getAtomPointers();

	// Decide starting point for loop over the positions in our system
	int start = 0;

	// Loop over all or a subset of the positions in the system.
	double energy  = 0.0;
	for (uint i = start; i < _sys.positionSize();i++){

		Position &pos = _sys.getPosition(i);

		// Skip Variable Position Side chains
		if (pos.getTotalNumberOfRotamers() > 1){
			continue;
		}

		// Skip self..
		if (i == _position) {
			continue;
		}

		map<string,double> energies = calculatePairwiseNonBondedEnergy(atoms1, pos.getAtomPointers());

		energy += energies["TOTAL"];
	}

	return energy;
}

/*
  Calculating a given rotamer vs the "fixed" part of the system

  If the rotamer given is "fixed", then we compute starting at _position+1 (used when computing fixed portion of energy)
  If the rotmaer given is variable then we compute over all positions (start at 0).

 */
double CharmmEnergyCalculator::calculateTemplateEnergy(System &_sys, int _position, int _rotamer,bool _calcAllForFixed){

	//cout << "Calculating template for position "<<_position<<" and rotamer "<<_rotamer<<endl;

	//Set active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position).getAtomPointers();

	// Decide starting point for loop over the positions in our system
	int start = _position+1;
	if (_calcAllForFixed || _sys.getPosition(_position).getTotalNumberOfRotamers() > 1){
		start = 0;
	}

	// Loop over all or a subset of the positions in the system.
	double energy  = 0.0;
	for (uint i = start; i < _sys.positionSize();i++){

		Position &pos = _sys.getPosition(i);

		// Skip Variable Position Side chains
		if (pos.getTotalNumberOfRotamers() > 1){
			continue;
		}

		// Skip self..
		if (i == _position) {
			continue;
		}

		//cout << "\t TEMPLATE POSITION IS "<<i<<" "<<pos.getResidueNumber()<<" "<<pos.getChainId()<<" "<<pos.getResidueName()<<endl;
		map<string,double> energies = calculatePairwiseEnergy(atoms1, pos.getAtomPointers());

		if (storeEneByType) {
			map<string,double>::iterator it;
			for (it = energies.begin();it != energies.end();it++){
				energiesByType[it->first] += it->second;
			}
		}


		energy += energies["TOTAL"];

	}
	


	return energy;
}

/*
  Calculating a given rotamer vs the variable portion of the rest of the system
 */
double CharmmEnergyCalculator::calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, vector< vector< vector< vector<double> > > > & rotamerInteractions, vector<uint> & currentAllRotamers){

	//Set active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position).getAtomPointers();

	// Loop over all positions in the system.
	double energy  = 0.0;
	for (uint i = 0; i < _sys.positionSize();i++){

		Position &pos = _sys.getPosition(i);

		// Skip self..
		if (i == _position) {
			continue;
		}

		// Skip fixed position Side chains
		if (pos.getTotalNumberOfRotamers() == 1){
			continue;
		}

		if (rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]] == MslTools::doubleMax) {
			map<string,double> energies = calculatePairwiseNonBondedEnergy(atoms1, pos.getAtomPointers());
			energy += energies["TOTAL"];
			rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]] = energies["TOTAL"];
			rotamerInteractions[i][currentAllRotamers[i]][_position][_rotamer] = energies["TOTAL"];
		}
		else {
			energy += rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]];
		}
	}
	
	return energy;
}

/*
  Calcuate pairwise energy between to rotamers.  These must define unique atom sets (don't give _position1 = _position2 && _rotamer1 == _rotamer2, that would be self.
  
 */
double CharmmEnergyCalculator::calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2){

	// Set active rotamer
	_sys.getPosition(_position1).setActiveRotamer(_rotamer1);

	// Get atoms of rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position1).getAtomPointers();

	// Set active rotamer
	_sys.getPosition(_position2).setActiveRotamer(_rotamer2);

	// Get atoms of rotamer
	AtomPointerVector &atoms2 = _sys.getPosition(_position2).getAtomPointers();
	

	//cout << "Pairwise "<<_position1<<" "<<_position2<<endl;
	map<string,double> energies = calculatePairwiseEnergy(atoms1, atoms2);

	if (storeEneByType) {
		map<string,double>::iterator it;
		for (it = energies.begin();it != energies.end();it++){
			energiesByType[it->first] += it->second;
		}
	}

	return energies["TOTAL"];
}



/*
  3 Scenerios
  
  1. Non-intersection                2 sets
  2. Single set                      1 set
  3. Partial to full intersection    < 2 sets


  This code assumes you don't have #3. Unexpected results will arise. There should be some sort of check, but there is not.
	  
 */
map<string,double> CharmmEnergyCalculator::calculatePairwiseEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet){

  /*
           NOTE: 
           This function does CHARMM on-the-fly non-bonded calculations, but current needs pre-computed bonded terms.  The bonded interactions are stored in this object, extracted from an EnergySet object.  

	   An alternative method is to use Atom.isBoundTo(), Atom.isOneThree(), Atom.isOneFour() ;  However, these properties of the Atom are set in CharmmSystemBuilder, which creates the Interactions of the EnergySet in the first place.
	   If we had another mechanism for populating Atom.isBoundTo(), isOneThree,isOneFour(), then we could use these attributes and not ever worry about ever calling the CharmmSystemBuilder...for TRUE on-the-fly calculations.
  */


	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 **********************************************************************/
	RandomNumberGenerator rng;
	unsigned int stamp = rng.getRandomInt(1000000);

	/********************************************************************************
	 *  About the cutoffs:
	 *
	 *   - if nonBondCutoffOn is not zero, a distance cutoff is applied to exclude interactions between far atoms
	 *   - the nonBondCutoffOn is the cutoff in which the switching function is applied to bring the energy
	 *     smoothly to zero.  
	 *   - the energy goes to zero at nonBondCutoffOff
	 *  That is:
	 *   - between 0 and nonBondCutoffOn E = full energy
	 *   - between nonBondCutoffOn and nonBondCutoffOff E = energy * switching function
	 *   - between nonBondCutoffOff and infinity E = 0.0
	 ********************************************************************************/


        // Resulting map of energies by type
	map<string,double> energies;
	energies["CHARMM_BOND"] = 0.0;
	energies["CHARMM_ANGL"] = 0.0;
	energies["CHARMM_ANGL"] = 0.0;
	energies["CHARMM_DIHE"] = 0.0;
	energies["CHARMM_U-BR"] = 0.0;
	energies["CHARMM_IMPR"] = 0.0;
	energies["CHARMM_VDW"] = 0.0;
	energies["CHARMM_ELEC"] = 0.0;

	for (uint i = 0; i < _a.size(); i++){

		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		for (uint j = startJ; j < _b.size();j++){

		  
			// Get Bonded Interactions
			vector<Interaction *> &bonds  = getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_BOND");
			vector<Interaction *> &angles = getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_ANGL");
			vector<Interaction *> &dihes  = getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_DIHE");
			vector<Interaction *> &urebr  = getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_U-BR");
			vector<Interaction *> &imprs  = getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_IMPR");

			// Collect the numbers of each type of interaction
			if (collectInteractionsFlag){
			  numberOfInteractionsUsed["CHARMM_BOND"] += bonds.size();
			  numberOfInteractionsUsed["CHARMM_ANGL"] += angles.size();
			  numberOfInteractionsUsed["CHARMM_DIHE"] += dihes.size();
			  numberOfInteractionsUsed["CHARMM_U-BR"] += urebr.size();
			  numberOfInteractionsUsed["CHARMM_IMPR"] += imprs.size();
			}
			

			// Sum up all bond energies associated with these 2 atoms (hopefully only one, not handling double/triple bonds in this code)
			for (uint n = 0; n < bonds.size();n++){
			  energies["CHARMM_BOND"] += bonds[n]->getEnergy();

			  // ERROR checking code, lookup this Interaction in stored hash table...
			  /*
			  map<string,double>::iterator it;
			  it = MSLOUT.getStaticLookup().find(bonds[n]->toString());
			  if (it == MSLOUT.getStaticLookup().end()){
			    fprintf(stdout,"NEW INTERACTION: %s %8.15f\n",it->first.c_str(),it->second);
			  } else {
			    fprintf(stdout, "INTERACTION[%s] = %8.15f %8.15f",it->first.c_str(),it->second,bonds[n]->getEnergy());
			    if (fabs(it->second - bonds[n]->getEnergy()) > 0.000000000000000000001){
			      fprintf(stdout," ERRRRRRRRRRRRRRRORRRRRRRRR\n");
			    } else {
				fprintf(stdout,"\n");
			    }

			  }
			  */
			}
			


			// Sum up all angle energies associated with these 2 atoms [ atoms must be first and last in angle ]
			for (uint n = 0; n < angles.size();n++){
			  energies["CHARMM_ANGL"] += angles[n]->getEnergy();
			}

			// Sum up all urey-bradley energies associated with these 2 atoms [ atoms must be first and last  ]
			for (uint n = 0; n < urebr.size();n++){
			  energies["CHARMM_U-BR"] += urebr[n]->getEnergy();
			}
			
			// Sum up all dihedral energies associated with these 2 atoms [ atoms must be first and last in dihedral ]
			for (uint n = 0; n < dihes.size();n++){
			  energies["CHARMM_DIHE"] += dihes[n]->getEnergy();

			  // Compute 1-4 VDW/ELEC  ** OLD WAY, now check isOneFour() **
			  // n == 0 because aromatics have more than 1 dihedral.
			  //if (n == 0 && bonds.size() == 0 && angles.size() == 0) {
			  //}
			}

			// Sum up all improper dihedral energies associated with these 2 atoms [ atoms must be first and last in dihedral ]
			for (uint n = 0; n < imprs.size();n++){
			  energies["CHARMM_IMPR"] += imprs[n]->getEnergy();
			}


			// Non-bonded calculations here..
			pair<double,double> nonBondedE(0.0,0.0);
			if (_a(i).isBoundTo(&_b(j)) || _a(i).isOneThree(&_b(j))) { continue; }
			else if (_a(i).isOneFour(&_b(j))) {


			  // Get pre-computed parameters for these types (better than getting vdwParam for each atom type and doing a sqrt right here)
			  vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());
			  double Kq_q1_q2_rescal = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * elec14factor;

			  nonBondedE = computeVdwElec(_a[i],_b[j],vdwParam[3],vdwParam[2],Kq_q1_q2_rescal,stamp);
			  
			} else {

			  // Get pre-computed parameters for these types (better than getting vdwParam for each atom type and doing a sqrt right here)
			  vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());
			  double Kq_q1_q2_rescal = (CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge());

			  nonBondedE = computeVdwElec(_a[i],_b[j],vdwParam[1],vdwParam[0],Kq_q1_q2_rescal,stamp);
			}
			
			if (collectInteractionsFlag) numberOfInteractionsUsed["CHARMM_VDW"]++;
			if (collectInteractionsFlag) numberOfInteractionsUsed["CHARMM_ELEC"]++;
			energies["CHARMM_VDW"]  += nonBondedE.first;
			energies["CHARMM_ELEC"] += nonBondedE.second;			
		}
	  }
	

	// STORE TOTAL ENERGY 
	map<string,double>::iterator it;
	double energyTotal = 0.0;
	for (it = energies.begin();it != energies.end();it++){
		energyTotal += it->second;
	}
	energies["TOTAL"] = energyTotal;

	return energies;
}

map<string,double> CharmmEnergyCalculator::calculatePairwiseNonBondedEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet){


	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 **********************************************************************/
	RandomNumberGenerator rng;
	unsigned int stamp = rng.getRandomInt(1000000);

	/********************************************************************************
	 *  About the cutoffs:
	 *
	 *   - if nonBondCutoffOn is not zero, a distance cutoff is applied to exclude interactions between far atoms
	 *   - the nonBondCutoffOn is the cutoff in which the switching function is applied to bring the energy
	 *     smoothly to zero.  
	 *   - the energy goes to zero at nonBondCutoffOff
	 *  That is:
	 *   - between 0 and nonBondCutoffOn E = full energy
	 *   - between nonBondCutoffOn and nonBondCutoffOff E = energy * switching function
	 *   - between nonBondCutoffOff and infinity E = 0.0
	 ********************************************************************************/


	map<string,double> energies;
	//bool clash = false;
	for (int i = (_a.size() - 1); i >= 0; i--){


		// Extract VDW parameters for atom1
		//vector<double> vdwParam1  = parReader->vdwParam(_a[i]->getType());
		//vector<double> vdwParam1;
		//if (!parReader->vdwParam(vdwParam1, _a[i]->getType())) {
		//	cerr << "WARNING 49319: VDW parameters not found for type " << _a[i]->getType() << " map<string,double> CharmmEnergyCalculator::calculatePairwiseNonBondedEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet)" << endl;
		//	continue;
		//}

		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		for (int j = (_b.size() - 1); j >= startJ; j--){

			// Extract VDW parameters for atom2
			//vector<double> vdwParam2 = parReader->vdwParam(_b[j]->getType());
			//vector<double> vdwParam2;
			//if (!parReader->vdwParam(vdwParam2, _b[j]->getType())) {
			//	cerr << "WARNING 49319: VDW parameters not found for type " << _a[i]->getType() << " map<string,double> CharmmEnergyCalculator::calculatePairwiseNonBondedEnergy(AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet)" << endl;
			//	continue;
			//}


			// Non-bonded calculations here..
			pair<double,double> nonBondedE(0.0,0.0);
			if (_a(i).isBoundTo(&_b(j)) || _a(i).isOneThree(&_b(j))) { continue; }
			else if (_a(i).isOneFour(&_b(j))) {


			  // Get pre-computed parameters for these types (better than getting vdwParam for each atom type and doing a sqrt right here)
			  vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());
			  double Kq_q1_q2_rescal = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * elec14factor;

			  nonBondedE = computeVdwElec(_a[i],_b[j],vdwParam[3],vdwParam[2],Kq_q1_q2_rescal,stamp);
			  
			} else {

			  // Get pre-computed parameters for these types (better than getting vdwParam for each atom type and doing a sqrt right here)
			  vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());
			  double Kq_q1_q2_rescal = (CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge());

			  nonBondedE = computeVdwElec(_a[i],_b[j],vdwParam[1],vdwParam[0],Kq_q1_q2_rescal,stamp);
			}
			energies["CHARMM_VDW"]  += nonBondedE.first;
			energies["CHARMM_ELEC"] += nonBondedE.second;
		}
	}

	// STORE TOTAL ENERGY 
	map<string,double>::iterator it;
	double energyTotal = 0.0;
	for (it = energies.begin();it != energies.end();it++){
		energyTotal += it->second;
	}
	energies["TOTAL"] = energyTotal;

	return energies;
}

/*
  Extract bonded interactions from an EnergySet
 */
bool CharmmEnergyCalculator::extractBondedInteractions(EnergySet &_es) {

  bool status = extractInteractions(_es, "CHARMM_BOND");
  if (!status){
    cerr << "ERROR CharmmEnergyCalculator::extractBondedInteractions CHARMM_BOND\n";
    return status;
  }
  status = extractInteractions(_es, "CHARMM_ANGL");
  if (!status){
    cerr << "ERROR CharmmEnergyCalculator::extractBondedInteractions CHARMM_ANGL\n";
    return status;
  }
  status = extractInteractions(_es, "CHARMM_DIHE");
  if (!status){
    cerr << "ERROR CharmmEnergyCalculator::extractBondedInteractions CHARMM_DIHE\n";
    return status;
  }
  status = extractInteractions(_es, "CHARMM_U-BR");
  if (!status){
    cerr << "ERROR CharmmEnergyCalculator::extractBondedInteractions CHARMM_U-BR\n";
    return status;
  }
  status = extractInteractions(_es, "CHARMM_IMPR");
  if (!status){
    cerr << "ERROR CharmmEnergyCalculator::extractBondedInteractions CHARMM_IMPR\n";
    return status;
  }

  return true;
}
bool CharmmEnergyCalculator::extractInteractions(EnergySet &_es, string _term){

  std::map<std::string, std::vector<Interaction*> > *allEnergyTerms = _es.getEnergyTerms();
  std::map<std::string, std::vector<Interaction*> >::iterator it;

  /*
  map<string,atomPairMap>::iterator k = pairInteractions.find(_term);
  if (k != pairInteractions.end()){
    //k->second.clear(); //unnecessary
    pairInteractions.erase(k);
  }
  */

  it = allEnergyTerms->find(_term);
  if (it == allEnergyTerms->end()){
    return false;
  }
  //cout << "Total " <<_term<<" interactions = "<< it->second.size()<<endl;
  for (uint i = 0; i < it->second.size();i++){
    Interaction *theInteraction = it->second[i];


    // Keep track of interactions by AtomPair key (Atom *, Atom *)
    vector<Atom *> &ats = theInteraction->getAtomPointers();
    map<string,atomPairMap >::iterator pairIt;
    atomPairMapIt atIt;  

    pairIt = pairInteractions.find(_term);

    AtomPair ab(ats[0],ats[ats.size()-1]);		
    if (pairIt == pairInteractions.end()){
			
      atomPairMap  tmp;
      tmp[ab].push_back(theInteraction);
      pairInteractions[_term] = tmp;
    } else {

      atIt = (pairIt->second).find(ab);
      
      // Atom pair ab, does not have a listed interaction in this term yet.
      if (atIt == (pairIt->second).end()){
				
	(pairIt->second)[ab].push_back(theInteraction);
      } else {

	/*
	// This output should not be suppressed, but it can be a lot of stuff to see.... debug/verbose flag?
	cerr << "WARNING 7724 CharmmEnergyCacluator::extractBondedInteractions()... THIS ATOM PAIR ALREADY HAS THIS TYPE OF INTERACTION ("<<name<<")"<<endl;

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


	(pairIt->second)[ab].push_back(theInteraction);
      }
    }
  }
  
  return true;
}



/*
  This function returns all interactions of a given type that contain Atom a and Atom b (where Atom a and Atom b are the first or last atoms in the interaction)
 */
vector<Interaction *> & CharmmEnergyCalculator::getEnergyInteractions(Atom *a, Atom *b, string _termName){


        //  Check that the type of energy term is in the pairInteraction map, return blank.
	map<string, atomPairMap >::iterator typeIt;
	typeIt = pairInteractions.find(_termName);
	if (typeIt == pairInteractions.end()){
	  return blank; // blank =  empty std::vector<Interaction *>
	}

	
	// Create an atom pair from the two input atoms
	AtomPair ab(a,b);

	// Search the atom pair in the pairInteractions map; do we have Interactions of _termName for Atoms a,b?
	atomPairMapIt atomIt;
	atomIt = (typeIt->second).find(ab);
	if (atomIt == (typeIt->second).end()){

	        // Double check : look up reverse atoms
		// RRRRRRRR . I shouldn't have to look things up in reverse. The map<stl::pair, ...,cmpOperator> , yet cmpOperator doesn
		AtomPair ba(b,a);
		atomIt = (typeIt->second).find(ba);
		if (atomIt != (typeIt->second).end()){
			return atomIt->second;
		}
			

		// If we get here, that means no Interactions exist for Atoms a,b of type _termType; blank =  empty std::vector<Interaction *>
		return blank;
	}


	// We should have a set interactions between atoms a and b.
	return atomIt->second;
	
}

/*
  Helper function for seeing the contents of "pairInteractions" map.
 */
void CharmmEnergyCalculator::dumpAtomMap(string _term){
  std::map<std::string, atomPairMap>::iterator it;
  it =  pairInteractions.find(_term);
  if (it != pairInteractions.end()){
    for (atomPairMapIt it2 = it->second.begin();it2 != it->second.end();it2++){
      cout << "Atom Pair: "<<it2->first.first->getAtomId()<< " "<<it2->first.second->getAtomId()<<endl;
      for (uint i = 0; i < it2->second.size();i++){
	cout << "Interaction : "<<it2->second[i]->toString()<<endl;
      }
    }
  } else {
    cout << "No terms to dump for "<<_term<<endl;
  }
}


pair<double,double> CharmmEnergyCalculator::computeVdwElec(Atom *_a, Atom *_b, double _rmin, double _emin, double _Kq_q1_q2_rescal, int _stamp){

  // Atoms further away than cutoff
  if (nonBondCutoffOn != 0.0 && _a->groupDistance(*_b,_stamp) > nonBondCutoffOff) {
    return pair<double,double>(0.0,0.0);
  }

  pair<double,double> result(0.0,0.0); // first = vdw, second = elec;
  double dist = _a->distance(*_b);
  double Kq_q1_q2_rescal_over_diel = _Kq_q1_q2_rescal / dielectricConstant;
  if (useRdielectric){
    Kq_q1_q2_rescal_over_diel /= dist;
  }

  if (nonBondCutoffOn != 0.0){
    result.first  = CharmmEnergy::instance()->LJSwitched(dist, vdwRescalingFactor*_rmin,_emin, _a->groupDistance(*_b,_stamp),nonBondCutoffOn,nonBondCutoffOff);
    result.second = CharmmEnergy::instance()->coulombEnerPrecomputedSwitched(dist,Kq_q1_q2_rescal_over_diel, _a->groupDistance(*_b,_stamp),nonBondCutoffOn,nonBondCutoffOff);
  } else {
    result.first  = CharmmEnergy::instance()->LJ(dist, vdwRescalingFactor*_rmin,_emin);
    result.second = CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);
  }

  return result;
}

		  /*
			        // Get pre-computed parameters for these types (better than getting vdwParam for each atom type and doing a sqrt right here)
				vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());

				if (collectInteractionsFlag) numberOfInteractionsUsed["CHARMM_VDW"]++;

				double dist = _a(i).distance(_b(j));
				double vdw = CharmmEnergy::instance()->LJ(dist, vdwRescalingFactor*vdwParam[1],vdwParam[0]);
				energies["CHARMM_VDW"] += vdw;


				double Kq_q1_q2_rescal_over_diel = (CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * 1.0 / dielectricConstant);
				if (useRdielectric){
				  Kq_q1_q2_rescal_over_diel /= dist;
				}
				energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);
				if (collectInteractionsFlag) numberOfInteractionsUsed["CHARMM_ELEC"]++;
			  */
