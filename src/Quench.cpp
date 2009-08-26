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

#include "Quench.h"

Quench::Quench()
: ape("/library/charmmTopPar/par_all22_prot.inp"), pec("/library/charmmTopPar/par_all22_prot.inp"), currentRotamers(0), currentAllRotamers(0), monomericSurroundEnergies(0)
{
	topfile = "/library/charmmTopPar/top_all22_prot.inp";
	parfile = "/library/charmmTopPar/par_all22_prot.inp";
	rotlib = "/library/rotlib/balanced/rotlib-balanced-200.txt";
	numberLargeRotamers = -1;
	numberSmallRotamers = -1;
}

Quench::Quench(string _topfile, string _parfile, string _rotlib)
: ape(_parfile), pec(_parfile), currentRotamers(0), currentAllRotamers(0), monomericSurroundEnergies(0)
{
	topfile = _topfile;
	parfile = _parfile;
	rotlib = _rotlib;
	numberLargeRotamers = -1;
	numberSmallRotamers = -1;
}

Quench::~Quench() {
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem) {
	setUpSystem(_initialSystem, _outputSystem, 10);
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers) {

	uint maxRotamer = _numRotamers - 1;
	
	PolymerSequence pseq(_initialSystem);

	// Build a new system from polymer sequence and create energySet from energy terms
	CharmmSystemBuilder CSB(topfile,parfile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(_outputSystem,pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	int numAssignedAtoms = _outputSystem.assignCoordinates(_initialSystem.getAtoms(),false);

	// Build the all atoms without coordinates (not in initial PDB)
	_outputSystem.buildAllAtoms();
	
	// Build rotamers
	SystemRotamerLoader sysRot(_outputSystem, rotlib);

	for (uint i = 0; i < _outputSystem.positionSize();i++){
		Position &pos = _outputSystem.getPosition(i);
		string posName = pos.getResidueName();

		int numRots = maxRotamer;

		// Override maxRot if numberLargeRotamers is set.
		if (numberLargeRotamers != -1){
			numRots = numberLargeRotamers;
		}

		if (numberSmallRotamers != -1 && (posName == "SER" || posName == "THR" || posName == "CYS" || posName == "VAL" || posName == "ASN" || posName == "ASP")){
			numRots = numberSmallRotamers;
		}
		if ((posName != "GLY") && (posName != "ALA") && (posName != "PRO")) {
			sysRot.loadRotamers(&pos, "BALANCED-200", pos.getResidueName(), 0, numRots);
		}
	}

	for (uint i = 0; i < _outputSystem.positionSize(); i++) {
		Position & posVar = _outputSystem.getPosition(i);

		if (posVar.getTotalNumberOfRotamers() == 1) { currentAllRotamers.push_back(0); continue; }

		stringstream ss;
		ss << _outputSystem.getPosition(i).getChainId() << " " << _outputSystem.getPosition(i).getResidueNumber();

		selfEnergies[ss.str()] = vector<double>(0);

		double minSelf = MslTools::doubleMax;
		uint minSelfPos = 0;

		for (uint j = 0; j < posVar.getTotalNumberOfRotamers(); j++) {
			double self = ape.calculateSelfEnergy(_outputSystem,i,j);
			selfEnergies[ss.str()].push_back(self);

			if (self < minSelf) {
				minSelf = self;
				minSelfPos = j;
			}
		}

		currentRotamers.push_back(minSelfPos);
		currentAllRotamers.push_back(minSelfPos);

		posVar.setActiveRotamer(minSelfPos);
	}

}

void Quench::runPreSetUpQuench(System & _mySystem){
	runPreSetUpQuench(_mySystem, 10);
}

void Quench::runPreSetUpQuench(System & _mySystem, uint _numIterations){
	vector < vector < vector < vector<double> > > > surroundEnergies(_mySystem.positionSize());
	for (uint i = 0; i < surroundEnergies.size(); i++) {
		surroundEnergies[i].resize(_mySystem.getPosition(i).getTotalNumberOfRotamers());
		for (uint j = 0; j < surroundEnergies[i].size(); j++) {
			surroundEnergies[i][j].resize(_mySystem.positionSize());
			for (uint k = 0; k < surroundEnergies[i][j].size(); k++) {
				surroundEnergies[i][j][k].resize(_mySystem.getPosition(k).getTotalNumberOfRotamers(),MslTools::doubleMax);
			}
		}
	}

	uint numLoops = 0;
	uint changes = MslTools::intMax;
	// So that we only index variable positions in currentRotamers
	uint variablePosCounter = 0;
	PDBWriter writer;
	while ((numLoops < _numIterations) && (changes > 0)) {

		changes = 0;
		variablePosCounter = 0;

		for (uint i = 0; i < _mySystem.positionSize(); i++) {
			Position & posVar = _mySystem.getPosition(i);

			if (posVar.getTotalNumberOfRotamers() > 1) {
				stringstream ss;
				ss << _mySystem.getPosition(i).getChainId() << " " << _mySystem.getPosition(i).getResidueNumber();

				double minTotal = MslTools::doubleMax;
				uint minTotalPos = 0;

				for (uint j = 0; j < posVar.getTotalNumberOfRotamers(); j++) {
					double self = selfEnergies[ss.str()][j];
					double surround = ape.calculateSurroundingEnergy(_mySystem,i,j,surroundEnergies, currentAllRotamers);
					double total = self + surround;

					if (total < minTotal) {
						minTotal = total;
						minTotalPos = j;
					}
				}

				if (minTotal == MslTools::doubleMax) { cerr << "Warning: Bad clash!" << endl; }

				if (currentRotamers[variablePosCounter] != minTotalPos) {
					changes++;
					currentRotamers[variablePosCounter] = minTotalPos;
					currentAllRotamers[i] = minTotalPos;
				}
				posVar.setActiveRotamer(minTotalPos);
				variablePosCounter++;
			}
		}

		numLoops++;
	}
}

System Quench::runQuench(System & _initialSystem){
	System mySys;
	setUpSystem(_initialSystem, mySys);

	runPreSetUpQuench(mySys, 10);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, uint _numIterations){
	System mySys;
	setUpSystem(_initialSystem, mySys);

	runPreSetUpQuench(mySys, _numIterations);
	return mySys;
}

void Quench::setUpMonomericSurroundEnergies(System & _mySystem) {
	monomericSurroundEnergies.resize(_mySystem.positionSize());
	for (uint i = 0; i < monomericSurroundEnergies.size(); i++) {
		monomericSurroundEnergies[i].resize(_mySystem.getPosition(i).getTotalNumberOfRotamers());
		for (uint j = 0; j < monomericSurroundEnergies[i].size(); j++) {
			monomericSurroundEnergies[i][j].resize(_mySystem.positionSize());
			for (uint k = 0; k < monomericSurroundEnergies[i][j].size(); k++) {
				monomericSurroundEnergies[i][j][k].resize(_mySystem.getPosition(k).getTotalNumberOfRotamers(),MslTools::doubleMax);
			}
		}
	}
}

double Quench::runPreSetUpQuenchOnDimer(System & _mySystem){
	return runPreSetUpQuenchOnDimer(_mySystem, 10);
}

double Quench::runPreSetUpQuenchOnDimer(System & _mySystem, uint _numIterations){
	vector < vector < vector < vector<double> > > > surroundEnergies(_mySystem.positionSize());
	for (uint i = 0; i < surroundEnergies.size(); i++) {
		surroundEnergies[i].resize(_mySystem.getPosition(i).getTotalNumberOfRotamers());
		for (uint j = 0; j < surroundEnergies[i].size(); j++) {
			surroundEnergies[i][j].resize(_mySystem.positionSize());
			for (uint k = 0; k < surroundEnergies[i][j].size(); k++) {
				if (_mySystem.getPosition(i).getChainId() == _mySystem.getPosition(k).getChainId()) { surroundEnergies[i][j][k] = monomericSurroundEnergies[i][j][k]; }
				else { surroundEnergies[i][j][k].resize(_mySystem.getPosition(k).getTotalNumberOfRotamers(),MslTools::doubleMax); }
			}
		}
	}

	uint numLoops = 0;
	uint changes = MslTools::intMax;
	// So that we only index variable positions in currentRotamers
	uint variablePosCounter = 0;
	PDBWriter writer;
	while ((numLoops < _numIterations) && (changes > 0)) {

		changes = 0;
		variablePosCounter = 0;

		for (uint i = 0; i < _mySystem.positionSize(); i++) {
			Position & posVar = _mySystem.getPosition(i);

			if (posVar.getTotalNumberOfRotamers() > 1) {
				stringstream ss;
				ss << _mySystem.getPosition(i).getChainId() << " " << _mySystem.getPosition(i).getResidueNumber();

				double minTotal = MslTools::doubleMax;
				uint minTotalPos = 0;

				for (uint j = 0; j < posVar.getTotalNumberOfRotamers(); j++) {
					double self = selfEnergies[ss.str()][j];
					double surround = ape.calculateSurroundingEnergy(_mySystem,i,j,surroundEnergies, currentAllRotamers);
					double total = self + surround;

					if (total < minTotal) {
						minTotal = total;
						minTotalPos = j;
					}
				}

				if (minTotal == MslTools::doubleMax) { cerr << "Warning: Bad clash!" << endl; }

				if (currentRotamers[variablePosCounter] != minTotalPos) {
					changes++;
					currentRotamers[variablePosCounter] = minTotalPos;
					currentAllRotamers[i] = minTotalPos;
				}
				posVar.setActiveRotamer(minTotalPos);
				variablePosCounter++;
			}
		}

		numLoops++;
	}

	for (uint i = 0; i < surroundEnergies.size(); i++) {
		for (uint j = 0; j < surroundEnergies[i].size(); j++) {
			for (uint k = 0; k < surroundEnergies[i][j].size(); k++) {
				if (_mySystem.getPosition(i).getChainId() == _mySystem.getPosition(k).getChainId()) { monomericSurroundEnergies[i][j][k] = surroundEnergies[i][j][k]; }
			}
		}
	}
	return pec.calculateStateEnergy(_mySystem,currentRotamers);

}

