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


#include "AtomicPairwiseEnergy.h"
#include "CharmmEnergy.h"

using namespace MSL;
using namespace std;


AtomicPairwiseEnergy::AtomicPairwiseEnergy(string _charmmParameterFile){
	storeEneByType  = false;
	storeEneByGroup = false;

	parReader = new CharmmParameterReader();
	parReader->open(_charmmParameterFile);
	parReader->read();
	parReader->close();
	parReader->createVdwParamPairs();

	// These kinds of parameters need to be set outside this object, otherwise
	//  all sorts of tests will start to fail.
	//CharmmEnergy::instance()->setDielectricConstant(80);


	vdwRescalingFactor = 1.0;
		
	elec14factor = 1.0;
	dielectricConstant = 1.0;
	useRdielectric = false;

	// Setting these cutoff values to 0.0 means NO cutoffs, compute ALL pairwise interactions regardless of distance.
	nonBondCutoffOn  = 0.0;
	nonBondCutoffOff = 0.0;
}

AtomicPairwiseEnergy::~AtomicPairwiseEnergy(){
	delete(parReader);

	
}

/*
  Calculate all pairwise energy between atoms within a residue
 */
double AtomicPairwiseEnergy::calculateSelfEnergy(System &_sys, int _position, int _rotamer){



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
	map<string,double> energies = calculatePairwiseEnergy(_sys, atoms1,atoms1,true);

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
double AtomicPairwiseEnergy::calculateBackgroundEnergy(System &_sys, int _position, int _rotamer){
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

		map<string,double> energies = calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtomPointers());

		energy += energies["TOTAL"];
	}

	return energy;
}

/*
  Calculating a given rotamer vs the "fixed" part of the system

  If the rotamer given is "fixed", then we compute starting at _position+1 (used when computing fixed portion of energy)
  If the rotmaer given is variable then we compute over all positions (start at 0).

 */
double AtomicPairwiseEnergy::calculateTemplateEnergy(System &_sys, int _position, int _rotamer,bool _calcAllForFixed){

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
		map<string,double> energies = calculatePairwiseEnergy(_sys, atoms1, pos.getAtomPointers());

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
double AtomicPairwiseEnergy::calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, vector< vector< vector< vector<double> > > > & rotamerInteractions, vector<uint> & currentAllRotamers){

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
			map<string,double> energies = calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtomPointers());
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
double AtomicPairwiseEnergy::calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2){

	// Set active rotamer
	_sys.getPosition(_position1).setActiveRotamer(_rotamer1);

	// Get atoms of rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position1).getAtomPointers();

	// Set active rotamer
	_sys.getPosition(_position2).setActiveRotamer(_rotamer2);

	// Get atoms of rotamer
	AtomPointerVector &atoms2 = _sys.getPosition(_position2).getAtomPointers();
	

	//cout << "Pairwise "<<_position1<<" "<<_position2<<endl;
	map<string,double> energies = calculatePairwiseEnergy(_sys, atoms1, atoms2);

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
map<string,double> AtomicPairwiseEnergy::calculatePairwiseEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet){

	/**********************************************************************
	 *  QUESTION: should atoms be checked for having coordinates?
	 **********************************************************************/

	// Get energy set for bonded terms. (Must have already been built somewhere else.. see CharmmSystemBuilder::buildSystem).
	EnergySet *eset = _sys.getEnergySet();


	map<string,double> energies;
	for (uint i = 0; i < _a.size(); i++){


		// Extract VDW parameters for atom1
		//vector<double> vdwParam1  = parReader->vdwParam(_a[i]->getType());

		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		for (uint j = startJ; j < _b.size();j++){

			// Extract VDW parameters for atom2
			//vector<double> vdwParam2 = parReader->vdwParam(_b[j]->getType());

			// Get Bonded Interactions
			vector<Interaction *> bonds  ;//= eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_BOND");
			vector<Interaction *> angles ;//= eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_ANGL");
			vector<Interaction *> dihes  ;//= eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_DIHE");
			vector<Interaction *> urebr  ;//= eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_U-BR");
			vector<Interaction *> imprs  ;//= eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_IMPR");

			// Comment this in for debuggingn
			//vector<Interaction *> vdwes  = eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_VDW");


			// Sum up all bond energies associated with these 2 atoms (hopefully only one, not handling double/triple bonds in this code)
			for (uint n = 0; n < bonds.size();n++){
				if (bonds[n] != NULL) {

					energies["CHARMM_BOND"] += bonds[n]->getEnergy();
					if (bonds[n]->getEnergy() > 100){
						cout << bonds[n]->toString()<<endl;
					}
					//
										
				}
			}


			// Sum up all angle energies associated with these 2 atoms [ atoms must be first and last in angle ]
			for (uint n = 0; n < angles.size();n++){
				if (angles[n] != NULL) {
					energies["CHARMM_ANGL"] += angles[n]->getEnergy();
				}
			}

			// Sum up all urey-bradley energies associated with these 2 atoms [ atoms must be first and last  ]
			for (uint n = 0; n < urebr.size();n++){
				if (urebr[n] != NULL) {
					energies["CHARMM_U-BR"] += urebr[n]->getEnergy();
				}
				
			}
			
			// Sum up all dihedral energies associated with these 2 atoms [ atoms must be first and last in dihedral ]
			for (uint n = 0; n < dihes.size();n++){
				if (dihes[n] != NULL) {
					energies["CHARMM_DIHE"] += dihes[n]->getEnergy();


					// Compute 1-4 VDW/ELEC
					// n == 0 because aromatics have more than 1 dihedral.
					if (n == 0 && bonds.size() == 0 && angles.size() == 0) {

						// Get pre-computed parameters for these types
						vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());

						double dist = _a(i).distance(_b(j));
						//double vdw = CharmmEnergy::instance()->LJ(dist, vdwParam1[3]+vdwParam2[3],sqrt(vdwParam1[2]*vdwParam2[2]));
						double vdw = CharmmEnergy::instance()->LJ(dist, vdwRescalingFactor*vdwParam[3],vdwParam[2]);
						/*
						if (vdw > 100){
							cout << "14Atoms have large VDW:"<<vdw<<"\n"<<_a(i).toString()<<"\n"<<_b(i).toString()<<endl;
						}
						*/
						energies["CHARMM_VDW"] += vdw;

						//double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * CharmmEnergy::instance()->getElec14factor() / CharmmEnergy::instance()->getDielectricConstant();
						double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * elec14factor / dielectricConstant;
						energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);


						
						// Comment this in for debugging.
						/*

						map<Interaction*,int>::iterator it;
						it = interactionsComputed.find(vdwes[0]);
						if (it ==  interactionsComputed.end()){
							interactionsComputed[vdwes[0]]=1;
						} else {
							cout << "WARNING.......... I've already computed this term.\n";
							cout << vdwes[0]<<" "<<vdwes[0]->toString()<<endl;
							exit(1);
						}
							
						*/


					}
				}
			}


			// Sum up all improper dihedral energies associated with these 2 atoms [ atoms must be first and last in dihedral ]
			for (uint n = 0; n < imprs.size();n++){
				if (imprs[n] != NULL) {
					energies["CHARMM_IMPR"] += imprs[n]->getEnergy();
				}
			}


			// NOW NON_BONDED... (2 sqrts!), I can remove them I think... LJ without using it.... pre-compute params...
			if (bonds.size() == 0 && angles.size() == 0 && dihes.size() == 0){

				// Get pre-computed parameters for these types
				vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());

				/*
				cout << "VDW PARAMS: ";
				for (uint p = 0; p < 4;p++){
					cout << vdwParam[p] <<" vs [ "<<vdwParam1[p]<<" , "<<vdwParam2[p]<<" ] === ";
				}
				cout <<endl;
				*/

				double dist = _a(i).distance(_b(j));
				//double vdw = CharmmEnergy::instance()->LJ(dist, vdwParam1[1]+vdwParam2[1],sqrt(vdwParam1[0]*vdwParam2[0]));
				double vdw = CharmmEnergy::instance()->LJ(dist, vdwRescalingFactor*vdwParam[1],vdwParam[0]);
				energies["CHARMM_VDW"] += vdw;
				/*
				if (vdw > 100){
					cout << "14Atoms have large VDW:"<<vdw<<"\n"<<_a(i).toString()<<"\n"<<_b(i).toString()<<endl;
				}
				*/

				double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * 1.0 / dielectricConstant;
				energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);


				/*
				map<Interaction*,int>::iterator it;
				it = interactionsComputed.find(vdwes[0]);
				if (it ==  interactionsComputed.end()){
					interactionsComputed[vdwes[0]]=1;
				} else {
					cout << "WARNING2.......... I've already computed this term.\n";
					cout << vdwes[0]<<" "<<vdwes[0]->toString()<<endl;
					exit(2);
				


				}
				*/



			}
			
			

			
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

map<string,double> AtomicPairwiseEnergy::calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet){


	// Get energy set for bonded terms. (Must have already been built somewhere else.. see CharmmSystemBuilder::buildSystem).
	//EnergySet *eset = _sys.getEnergySet();

	/**********************************************************************
	 * the stamp is a random number that is used to recall the center of each atom
	 * group to avoid to calculate it multiple times (essentially the center is calculated
	 * the first time a group distance is called and the value is cached and returned directly
	 * if the groupDistance function is called on the same group with the same stamp
	 **********************************************************************/
	//unsigned int stamp = MslTools::getRandomInt(1000000);
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
		vector<double> vdwParam1;
		if (!parReader->vdwParam(vdwParam1, _a[i]->getType())) {
			cerr << "WARNING 49319: VDW parameters not found for type " << _a[i]->getType() << " map<string,double> AtomicPairwiseEnergy::calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet)" << endl;
			continue;
		}

		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		for (int j = (_b.size() - 1); j >= startJ; j--){

			// Extract VDW parameters for atom2
			//vector<double> vdwParam2 = parReader->vdwParam(_b[j]->getType());
			vector<double> vdwParam2;
			if (!parReader->vdwParam(vdwParam2, _b[j]->getType())) {
				cerr << "WARNING 49319: VDW parameters not found for type " << _a[i]->getType() << " map<string,double> AtomicPairwiseEnergy::calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet)" << endl;
				continue;
			}

			// Atoms further away than cutoff
			if (nonBondCutoffOn != 0.0 && _a(i).groupDistance(_b(j),stamp) > nonBondCutoffOff) {
				  continue;
			}
			
			if (_a(i).isBoundTo(&_b(j)) || _a(i).isOneThree(&_b(j))) {}
                        else if (_a(i).isOneFour(&_b(j))) {

						double dist = _a(i).distance(_b(j));						    

						double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * elec14factor / dielectricConstant;
						if (nonBondCutoffOn != 0.0){
							energies["CHARMM_VDW"]  += CharmmEnergy::instance()->LJSwitched(dist, vdwParam1[3]+vdwParam2[3],sqrt(vdwParam1[2]*vdwParam2[2]),_a(i).groupDistance(_b(j),stamp),nonBondCutoffOn,nonBondCutoffOff);
							energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputedSwitched(dist,Kq_q1_q2_rescal_over_diel,_a(i).groupDistance(_b(j),stamp),nonBondCutoffOn,nonBondCutoffOff);
						} else {
							energies["CHARMM_VDW"]  += CharmmEnergy::instance()->LJ(dist, vdwParam1[3]+vdwParam2[3],sqrt(vdwParam1[2]*vdwParam2[2]));
							energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);
						}
			}
			else {
				double dist = _a(i).distance(_b(j));
				double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * 1.0 / dielectricConstant;
				if (nonBondCutoffOn != 0.0){
	    				energies["CHARMM_VDW"]  += CharmmEnergy::instance()->LJSwitched(dist, vdwParam1[1]+vdwParam2[1],sqrt(vdwParam1[0]*vdwParam2[0]),_a(i).groupDistance(_b(j),stamp),nonBondCutoffOn,nonBondCutoffOff);
					energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputedSwitched(dist,Kq_q1_q2_rescal_over_diel,_a(i).groupDistance(_b(j),stamp),nonBondCutoffOn,nonBondCutoffOff);
				} else {
					energies["CHARMM_VDW"]  += CharmmEnergy::instance()->LJ(dist, vdwParam1[1]+vdwParam2[1],sqrt(vdwParam1[0]*vdwParam2[0]));
					energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);
				}
			}
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



/*map<string,double> AtomicPairwiseEnergy::calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet){


	// Get energy set for bonded terms. (Must have already been built somewhere else.. see CharmmSystemBuilder::buildSystem).
	EnergySet *eset = _sys.getEnergySet();


	map<string,double> energies;
	//bool clash = false;
	for (int i = (_a.size() - 1); i >= 0; i--){


		// Extract VDW parameters for atom1
		vector<double> vdwParam1  = parReader->vdwParam(_a[i]->getType());

		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		for (int j = (_b.size() - 1); j >= startJ; j--){

			// Extract VDW parameters for atom2
			vector<double> vdwParam2 = parReader->vdwParam(_b[j]->getType());
			
			vector<Interaction *> bonds  = eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_BOND");
			vector<Interaction *> angles = eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_ANGL");
			vector<Interaction *> dihes  = eset->getEnergyInteractions(_a[i], _b[j], (string)"CHARMM_DIHE");

			// Sum up all dihedral energies associated with these 2 atoms [ atoms must be first and last in dihedral ]
			for (uint n = 0; n < dihes.size();n++){
				if (dihes[n] != NULL) {
					// Compute 1-4 VDW/ELEC
					// n == 0 because aromatics have more than 1 dihedral.
					if (n == 0 && bonds.size() == 0 && angles.size() == 0) {

						// Get pre-computed parameters for these types
						//vector<double> * vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());

						double dist = _a(i).distance(_b(j));
						double vdw = CharmmEnergy::instance()->LJ(dist, vdwParam1[3]+vdwParam2[3],sqrt(vdwParam1[2]*vdwParam2[2]));
						//double vdw = CharmmEnergy::instance()->LJ(dist, (*vdwParam)[3],(*vdwParam)[2]);
						energies["CHARMM_VDW"] += vdw;

						double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * CharmmEnergy::instance()->getElec14factor() / CharmmEnergy::instance()->getDielectricConstant();
						energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);

					}
				}
			}


			// NOW NON_BONDED... (2 sqrts!), I can remove them I think... LJ without using it.... pre-compute params...
			if (bonds.size() == 0 && angles.size() == 0 && dihes.size() == 0){

				// Get pre-computed parameters for these types
				//vector<double> vdwParam = parReader->vdwParamPair(_a[i]->getType(),_b[j]->getType());

				double dist = _a(i).distance(_b(j));
				//double vdw = CharmmEnergy::instance()->LJ(dist, (*vdwParam)[1],(*vdwParam)[0]);
				double vdw = CharmmEnergy::instance()->LJ(dist, vdwParam1[1]+vdwParam2[1],sqrt(vdwParam1[0]*vdwParam2[0]));
				energies["CHARMM_VDW"] += vdw;

				double Kq_q1_q2_rescal_over_diel = CharmmEnergy::Kq * _a[i]->getCharge() * _b[j]->getCharge() * 1.0 / CharmmEnergy::instance()->getDielectricConstant();
				energies["CHARMM_ELEC"] += CharmmEnergy::instance()->coulombEnerPrecomputed(dist, Kq_q1_q2_rescal_over_diel);
			}
			
			

			//if (energies["CHARMM_VDW"] > 1e3) { clash = true; }
		}
		//if (clash) { energies["TOTAL"] = MslTools::doubleMax; break; }
	}

	// STORE TOTAL ENERGY 
	map<string,double>::iterator it;
	double energyTotal = 0.0;
	for (it = energies.begin();it != energies.end();it++){
		energyTotal += it->second;
	}
	energies["TOTAL"] = energyTotal;

	return energies;
}*/




