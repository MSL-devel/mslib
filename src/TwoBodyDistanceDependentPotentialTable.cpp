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

#include "TwoBodyDistanceDependentPotentialTable.h"
#include "MslExceptions.h"

using namespace MSL;
using namespace std;



void TwoBodyDistanceDependentPotentialTable::readPotentialTable(string _fileName, string _potentialName){
	fileName      = _fileName;
	potentialName = _potentialName;

	
	reader.open(_fileName);
	reader.read(this);
	reader.close();


}

void TwoBodyDistanceDependentPotentialTable::readPotentialTable(string _fileName){
	fileName      = _fileName;
	
	reader.open(_fileName);
	reader.read(this);
	reader.close();


	cout << "Potential: "<<potentialName<<endl;
	cout << "Residue Skipping Number: "<<residueSkippingNumber<<endl;
	cout << "Bins: "<<endl;
	for (uint i = 0 ; i < distBins.size();i++){
		cout << "\t"<<distBins[i].index<<" [ "<<distBins[i].startDistance<<" , "<<distBins[i].endDistance<<" ] "<<endl;
	}

}


void TwoBodyDistanceDependentPotentialTable::addPotential(string _name1, string _name2, int _distBin, double _value){

	if (_distBin >=  distBins.size() ){
		stringstream ss;
		ss << "TwoBodyDistanceDependentPotentialTable::addPotential() , asking for a distance bin that is not specified: "<<_distBin<<" size of vector : "<<distBins.size();
		throw MslSizeException(ss.str());
	}

	
	stringstream key;
	key << _name1 <<":"<<_name2<<":"<<_distBin;

	PotentialTable::addPotential(key.str(),_value);
}

double TwoBodyDistanceDependentPotentialTable::getPotential(string _name1, string _name2, int _distBin){

	if (distBins.size() < _distBin){
		stringstream ss;
		ss << "TwoBodyDistanceDependentPotentialTable::getPotential() , asking for a distance bin that is not specified: "<<_distBin<<" size of vector : "<<distBins.size();
		throw MslSizeException(ss.str());
	}

	
	stringstream key;
	key << _name1 <<":"<<_name2<<":"<<_distBin;

	return PotentialTable::getPotential(key.str());
}

/*
  Calculate all pairwise energy between atoms within a residue
 */
double TwoBodyDistanceDependentPotentialTable::calculateSelfEnergy(System &_sys, int _position, int _rotamer){
	return 0;
}

/*
  Calculating a given rotamer vs the "fixed" part of the system

  If the rotamer given is "fixed", then we compute starting at _position+1 (used when computing fixed portion of energy)
  If the rotmaer given is variable then we compute over all positions (start at 0).
 */
double TwoBodyDistanceDependentPotentialTable::calculateTemplateEnergy(System &_sys, int _position, int _rotamer, bool _calcAllForFixed, bool _countLocalSCBB){
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

		energy += calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtomPointers(), _countLocalSCBB);
	}

	return energy;
}

/*
  Calcuate pairwise energy between to rotamers.  These must define unique atom sets (don't give _position1 = _position2 && _rotamer1 == _rotamer2, that would be self.
 */
double TwoBodyDistanceDependentPotentialTable::calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2, bool _countLocalSCBB){

	// Set active rotamer
	_sys.getPosition(_position1).setActiveRotamer(_rotamer1);

	// Get atoms of rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position1).getAtomPointers();

	// Set active rotamer
	_sys.getPosition(_position2).setActiveRotamer(_rotamer2);

	// Get atoms of rotamer
	AtomPointerVector &atoms2 = _sys.getPosition(_position2).getAtomPointers();

	//cout << "Pairwise "<<_position1<<" "<<_position2<<endl;
	double energy = calculatePairwiseNonBondedEnergy(_sys, atoms1, atoms2, _countLocalSCBB);

	return energy;
}

/*
  Calculate energy of rotamer with positions that are fixed.  Much like calculateTemplateEnergy, but use only non-bonded energies
 */
double TwoBodyDistanceDependentPotentialTable::calculateBackgroundEnergy(System &_sys, int _position, int _rotamer, bool _countLocalSCBB) {
	//Set active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomPointerVector &atoms1 = _sys.getPosition(_position).getAtomPointers();

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

		energy += calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtomPointers(), _countLocalSCBB);
	}

	return energy;
}

double TwoBodyDistanceDependentPotentialTable::calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, vector< vector< vector< vector<double> > > > & rotamerInteractions, vector<uint> & currentAllRotamers, bool _countLocalSCBB) {
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
			double energies = calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtomPointers(), _countLocalSCBB);
			energy += energies;
			rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]] = energies;
			rotamerInteractions[i][currentAllRotamers[i]][_position][_rotamer] = energies;
		}
		else {
			energy += rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]];
		}
	}
	
	return energy;
}

double TwoBodyDistanceDependentPotentialTable::calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet, bool _countLocalSCBB){

	double energies = 0.;
	for (int i = (_a.size() - 1); i >= 0; i--){
		// Adjust starting point for inner loop depending if _a == _b or not.
		int startJ = 0;
		if (_sameSet){
			startJ = i+1;
		}
		if (_a(i).getName().substr(0,1) == "H") { continue; }
		for (int j = (_b.size() - 1); j >= startJ; j--){
			if (_b(j).getName().substr(0,1) == "H") { continue; }
			else if ((_a(i).getChainId() == _b(j).getChainId()) && isBackbone(_a(i).getName()) && isBackbone(_b(j).getName()) && (abs(_a(i).getResidueNumber() - _b(j).getResidueNumber()) <= 7)) {}
			else if ((!_countLocalSCBB) && (_a(i).getChainId() == _b(j).getChainId()) && (abs(_a(i).getResidueNumber() - _b(j).getResidueNumber()) <= 7)) {}
			else if ((_a(i).getChainId() == _b(j).getChainId()) && (abs(_a(i).getResidueNumber() == _b(j).getResidueNumber()))) {}
			//else if ((_a(i).getChainId() == _b(j).getChainId()) && (isBackbone(_a(i).getName()) || isBackbone(_b(j).getName())) && (abs(_a(i).getResidueNumber() - _b(j).getResidueNumber()) <= 7)) {}
			//else if ((_a(i).getChainId() == _b(j).getChainId()) && (abs(_a(i).getResidueNumber() - _b(j).getResidueNumber()) <= 7)) {}
			else {
				double kb = 0;
				double dist = _a(i).distance(_b(j));
				if (dist < getMinDistCutoff()) {
					kb = getValueBelowCutoff();
				}
				else if (dist > getMaxDistCutoff()) { kb = getValueAboveCutoff(); }
				else {
					stringstream name1;
					name1 << _a(i).getResidueName() << "_" << _a(i).getName();
					stringstream name2;
					name2 << _b(j).getResidueName() << "_" << _b(j).getName();

					double distBin = getBin(dist);

					kb = getPotential(name1.str(),name2.str(),int(distBin));
				}
				energies += kb;
			}
		}
	}

	return energies;
}

bool TwoBodyDistanceDependentPotentialTable::isBackbone(string atomname) {
	if ((atomname == "N") || (atomname == "C") || (atomname == "CA") || (atomname == "O")) { return true; }
        else { return false; }
}
