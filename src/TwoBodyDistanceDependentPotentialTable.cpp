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

#include "TwoBodyDistanceDependentPotentialTable.h"
#include "MslExceptions.h"


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

double TwoBodyDistanceDependentPotentialTable::calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, vector< vector< vector< vector<double> > > > & rotamerInteractions, vector<uint> & currentAllRotamers) {
	//Set active rotamer
	_sys.getPosition(_position).setActiveRotamer(_rotamer);

	// Get atoms of active rotamer
	AtomVector &atoms1 = _sys.getPosition(_position).getAtoms();

	// Loop over all positions in the system.
	double energy  = 0.0;
	for (uint i = 0; i < _sys.positionSize();i++){

		Position &pos = _sys.getPosition(i);

		// Skip self..
		if (i == _position) {
			continue;
		}

		if (rotamerInteractions[_position][_rotamer][i][currentAllRotamers[i]] == MslTools::doubleMax) {
			double energies = calculatePairwiseNonBondedEnergy(_sys, atoms1, pos.getAtoms());
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

double TwoBodyDistanceDependentPotentialTable::calculatePairwiseNonBondedEnergy(System &_sys, AtomVector &_a, AtomVector &_b, bool _sameSet){

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

					kb = getPotential(name1.str(),name2.str(),distBin);
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
