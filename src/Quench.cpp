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

#include "Quench.h"

using namespace MSL;
using namespace std;

#include "SysEnv.h"
static SysEnv SYSENV;

Quench::Quench()
	: calculator(SYSENV.getEnv("MSL_CHARMM_PAR")), currentRotamers(0), currentAllRotamers(0), monomericSurroundEnergies(0)
{
	topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
	rotlib =  SYSENV.getEnv("MSL_ROTLIB");
	numberLargeRotamers = -1;
	numberSmallRotamers = -1;
	rotLevel = "";

	calculator.setNonBondedCutoffs(8,12); // 0->8 is full 8-12 is switched, >12 is 0.
}

Quench::Quench(string _topfile, string _parfile, string _rotlib)
: calculator(_parfile), currentRotamers(0), currentAllRotamers(0), monomericSurroundEnergies(0)
{
	topfile = _topfile;
	parfile = _parfile;
	rotlib = _rotlib;
	numberLargeRotamers = -1;
	numberSmallRotamers = -1;
	rotLevel = "";

	calculator.setNonBondedCutoffs(8,12); // 0->8 is full 8-12 is switched, >12 is 0.
}

Quench::~Quench() {
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem) {
	setUpSystem(_initialSystem, _outputSystem, 10);
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, vector<int> & variablePositions) {
	setUpSystem(_initialSystem, _outputSystem, 10, variablePositions);
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, vector<int> & variablePositions) {

	uint maxRotamer = _numRotamers - 1;
	
	PolymerSequence pseq(_initialSystem);

	// Build a new system from polymer sequence and create energySet from energy terms
	CharmmSystemBuilder CSB(_outputSystem,topfile,parfile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	_outputSystem.assignCoordinates(_initialSystem.getAtomPointers(),false);

	// Build the all atoms with coordinates (in initial PDB)
	_outputSystem.buildAllAtoms();


	// Build rotamers
	SystemRotamerLoader sysRot(_outputSystem, rotlib);

	for (uint i = 0; i < variablePositions.size();i++){
		Position &pos = _outputSystem.getPosition(variablePositions[i]);
		string posName = pos.getResidueName();

		// Override maxRot if rotLevel is set.
		if (rotLevel != ""){
			if (sysRot.getRotamerLibrary()->getLevel(rotLevel, posName) > 1) {
				sysRot.loadRotamers(&pos, pos.getResidueName(), 0, rotLevel, "");
				if (!sysRot.loadRotamers(&pos, posName, rotLevel,"")) {
					cerr << "Cannot load rotamers " << pos.getPositionId() << "," << posName << endl;
					exit(0);
				}
			}
		} else {
			int numRots = maxRotamer;

			numRots = numberLargeRotamers;
			// Override maxRot if numberLargeRotamers is set.
			if (numberLargeRotamers != -1){
				numRots = numberLargeRotamers;
			}

			if (numberSmallRotamers != -1 && (posName == "SER" || posName == "THR" || posName == "CYS" || posName == "VAL" || posName == "ASN" || posName == "ASP")){
				numRots = numberSmallRotamers;
			}
			if ((posName != "GLY") && (posName != "ALA") && (posName != "PRO")) {
				sysRot.loadRotamers(&pos, posName, 0, numRots, "");
			}
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
			double self = calculator.calculateSelfEnergy(_outputSystem,i,j);
			self += calculator.calculateBackgroundEnergy(_outputSystem,i,j);
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

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers) {

	uint maxRotamer = _numRotamers - 1;
	
	PolymerSequence pseq(_initialSystem);

	// Build a new system from polymer sequence and create energySet from energy terms
	CharmmSystemBuilder CSB(_outputSystem,topfile,parfile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	_outputSystem.assignCoordinates(_initialSystem.getAtomPointers(),false);
	
	// Build the all atoms with coordinates (in initial PDB)
	_outputSystem.buildAllAtoms();
	
	// Build rotamers
	SystemRotamerLoader sysRot(_outputSystem, rotlib);


	for (uint i = 0; i < _outputSystem.positionSize();i++){
		Position &pos = _outputSystem.getPosition(i);
		string posName = pos.getResidueName();

		// Override maxRot if rotLevel is set.
		if (rotLevel != ""){
			if (sysRot.getRotamerLibrary()->getLevel(rotLevel, posName) > 1) {
				sysRot.loadRotamers(&pos, pos.getResidueName(), 0, rotLevel, "");
				if (!sysRot.loadRotamers(&pos, posName, rotLevel,"")) {
					cerr << "Cannot load rotamers " << pos.getPositionId() << "," << posName << endl;
					exit(0);
				}
			}
		} else {
			int numRots = maxRotamer;

			numRots = numberLargeRotamers;
			// Override maxRot if numberLargeRotamers is set.
			if (numberLargeRotamers != -1){
				numRots = numberLargeRotamers;
			}

			if (numberSmallRotamers != -1 && (posName == "SER" || posName == "THR" || posName == "CYS" || posName == "VAL" || posName == "ASN" || posName == "ASP")){
				numRots = numberSmallRotamers;
			}
			if ((posName != "GLY") && (posName != "ALA") && (posName != "PRO")) {
				sysRot.loadRotamers(&pos, posName, 0, numRots, "");
			}
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
			double self = calculator.calculateSelfEnergy(_outputSystem,i,j);
			self += calculator.calculateBackgroundEnergy(_outputSystem,i,j);
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

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, TwoBodyDistanceDependentPotentialTable & tbd) {
	setUpSystem(_initialSystem, _outputSystem, 10, tbd);
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd) {
	setUpSystem(_initialSystem, _outputSystem, 10, variablePositions, tbd);
}

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd) {

	uint maxRotamer = _numRotamers - 1;
	
	PolymerSequence pseq(_initialSystem);

	// Build a new system from polymer sequence and create energySet from energy terms
	CharmmSystemBuilder CSB(_outputSystem,topfile,parfile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	_outputSystem.assignCoordinates(_initialSystem.getAtomPointers(),false);

	// Build the all atoms with coordinates (in initial PDB)
	_outputSystem.buildAllAtoms();
	
	// Build rotamers
	SystemRotamerLoader sysRot(_outputSystem, rotlib);

	for (uint i = 0; i < variablePositions.size();i++){
		Position &pos = _outputSystem.getPosition(variablePositions[i]);
		string posName = pos.getResidueName();

		// Override maxRot if rotLevel is set.
		if (rotLevel != ""){
			if (sysRot.getRotamerLibrary()->getLevel(rotLevel, posName) > 1) {
				sysRot.loadRotamers(&pos, pos.getResidueName(), 0, rotLevel, "");
				if (!sysRot.loadRotamers(&pos, posName, rotLevel,"")) {
					cerr << "Cannot load rotamers " << pos.getPositionId() << "," << posName << endl;
					exit(0);
				}
			}
		} else {
			int numRots = maxRotamer;

			numRots = numberLargeRotamers;
			// Override maxRot if numberLargeRotamers is set.
			if (numberLargeRotamers != -1){
				numRots = numberLargeRotamers;
			}

			if (numberSmallRotamers != -1 && (posName == "SER" || posName == "THR" || posName == "CYS" || posName == "VAL" || posName == "ASN" || posName == "ASP")){
				numRots = numberSmallRotamers;
			}
			if ((posName != "GLY") && (posName != "ALA") && (posName != "PRO")) {
				sysRot.loadRotamers(&pos, posName, 0, numRots, "");
			}
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
			double self = tbd.calculateSelfEnergy(_outputSystem,i,j);
			self += tbd.calculateBackgroundEnergy(_outputSystem,i,j,true);
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

void Quench::setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, TwoBodyDistanceDependentPotentialTable & tbd) {

	uint maxRotamer = _numRotamers - 1;
	
	PolymerSequence pseq(_initialSystem);

	// Build a new system from polymer sequence and create energySet from energy terms
	CharmmSystemBuilder CSB(_outputSystem,topfile,parfile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Apply coordinates from structure from PDB
	//  Variable positions w/o WT identity don't get built properly.
	_outputSystem.assignCoordinates(_initialSystem.getAtomPointers(),false);
	
	// Build the all atoms with coordinates (in initial PDB)
	_outputSystem.buildAllAtoms();
	
	// Build rotamers
	SystemRotamerLoader sysRot(_outputSystem, rotlib);

	for (uint i = 0; i < _outputSystem.positionSize();i++){
		Position &pos = _outputSystem.getPosition(i);
		string posName = pos.getResidueName();

		// Override maxRot if rotLevel is set.
		if (rotLevel != ""){
			if (sysRot.getRotamerLibrary()->getLevel(rotLevel, posName) > 1) {
				sysRot.loadRotamers(&pos, pos.getResidueName(), 0, rotLevel, "");
				if (!sysRot.loadRotamers(&pos, posName, rotLevel,"")) {
					cerr << "Cannot load rotamers " << pos.getPositionId() << "," << posName << endl;
					exit(0);
				}
			}
		} else {
			int numRots = maxRotamer;

			numRots = numberLargeRotamers;
			// Override maxRot if numberLargeRotamers is set.
			if (numberLargeRotamers != -1){
				numRots = numberLargeRotamers;
			}

			if (numberSmallRotamers != -1 && (posName == "SER" || posName == "THR" || posName == "CYS" || posName == "VAL" || posName == "ASN" || posName == "ASP")){
				numRots = numberSmallRotamers;
			}
			if ((posName != "GLY") && (posName != "ALA") && (posName != "PRO")) {
				sysRot.loadRotamers(&pos, posName, 0, numRots, "");
			}
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
			double self = tbd.calculateSelfEnergy(_outputSystem,i,j);
			self += tbd.calculateBackgroundEnergy(_outputSystem,i,j,true);
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

void Quench::runPreSetUpQuench(System & _mySystem, TwoBodyDistanceDependentPotentialTable & tbd){
	runPreSetUpQuench(_mySystem, 10, tbd);
}

void Quench::runPreSetUpQuench(System & _mySystem, uint _numIterations, TwoBodyDistanceDependentPotentialTable & tbd){
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
	PDBWriter writer;
	RandomNumberGenerator rng(false);
	//rng.setRNGTimeBasedSeed();
	vector<int> shuffledOrder;
	vector<int> variablePosOrder;
	for (uint i = 0; i < _mySystem.positionSize(); i++) {
		Position & posVar = _mySystem.getPosition(i);
		variablePosOrder.push_back(shuffledOrder.size());
		if (posVar.getTotalNumberOfRotamers() > 1) {
			shuffledOrder.push_back(i);
		}
	}
	int totalNum = shuffledOrder.size();
	random_shuffle(shuffledOrder.begin(),shuffledOrder.end(),rng);
	double overallEnergy = 0.;
	while ((numLoops < _numIterations) && (changes > 0)) {

		overallEnergy = 0.;
		changes = 0;

		for (uint i = 0; i < totalNum; i++) {
			int thisPos = shuffledOrder[i];
			Position & posVar = _mySystem.getPosition(thisPos);

			stringstream ss;
			ss << _mySystem.getPosition(thisPos).getChainId() << " " << _mySystem.getPosition(thisPos).getResidueNumber();

			double minTotal = MslTools::doubleMax;
			uint minTotalPos = 0;

			for (uint j = 0; j < posVar.getTotalNumberOfRotamers(); j++) {
				double self = selfEnergies[ss.str()][j];
				double surround = tbd.calculateSurroundingEnergy(_mySystem,thisPos,j,surroundEnergies, currentAllRotamers, true);
				double total = self + surround;

				if (total < minTotal) {
					minTotal = total;
					minTotalPos = j;
				}
			}

			if (minTotal == MslTools::doubleMax) { cerr << "Warning: Bad clash!" << endl; }
			overallEnergy += minTotal;

			if (currentRotamers[variablePosOrder[thisPos]] != minTotalPos) {
				changes++;
				currentRotamers[variablePosOrder[thisPos]] = minTotalPos;
				currentAllRotamers[thisPos] = minTotalPos;
			}
			posVar.setActiveRotamer(minTotalPos);
			cout << thisPos << " " << minTotal << endl;
		}

		numLoops++;
		cout << "New loop: " << numLoops << " " << overallEnergy << endl;
	}
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
	PDBWriter writer;
	RandomNumberGenerator rng(false);
	//rng.setRNGTimeBasedSeed();
	vector<int> shuffledOrder;
	vector<int> variablePosOrder;
	for (uint i = 0; i < _mySystem.positionSize(); i++) {
		Position & posVar = _mySystem.getPosition(i);
		variablePosOrder.push_back(shuffledOrder.size());
		if (posVar.getTotalNumberOfRotamers() > 1) {
			shuffledOrder.push_back(i);
		}
	}
	int totalNum = shuffledOrder.size();
	random_shuffle(shuffledOrder.begin(),shuffledOrder.end(),rng);
	double overallEnergy = 0.;
	while ((numLoops < _numIterations) && (changes > 0)) {

		overallEnergy = 0.;
		changes = 0;

		for (uint i = 0; i < totalNum; i++) {
			int thisPos = shuffledOrder[i];
			Position & posVar = _mySystem.getPosition(thisPos);

			stringstream ss;
			ss << _mySystem.getPosition(thisPos).getChainId() << " " << _mySystem.getPosition(thisPos).getResidueNumber();

			double minTotal = MslTools::doubleMax;
			uint minTotalPos = 0;

			for (uint j = 0; j < posVar.getTotalNumberOfRotamers(); j++) {
				double self = selfEnergies[ss.str()][j];
				double surround = calculator.calculateSurroundingEnergy(_mySystem,thisPos,j,surroundEnergies, currentAllRotamers);
				double total = self + surround;

				if (total < minTotal) {
					minTotal = total;
					minTotalPos = j;
				}
			}

			if (minTotal == MslTools::doubleMax) { cerr << "Warning: Bad clash!" << endl; }
			overallEnergy += minTotal;

			if (currentRotamers[variablePosOrder[thisPos]] != minTotalPos) {
				changes++;
				currentRotamers[variablePosOrder[thisPos]] = minTotalPos;
				currentAllRotamers[thisPos] = minTotalPos;
			}
			posVar.setActiveRotamer(minTotalPos);
			cout << thisPos << " " << minTotal << endl;
		}

		numLoops++;
		cout << "New loop: " << numLoops << " " << overallEnergy << endl;
	}
}

System Quench::runQuench(System & _initialSystem){
	System mySys;
	setUpSystem(_initialSystem, mySys);

	runPreSetUpQuench(mySys, 10);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, TwoBodyDistanceDependentPotentialTable & tbd){
	System mySys;
	setUpSystem(_initialSystem, mySys, tbd);

	runPreSetUpQuench(mySys, 10, tbd);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, vector<int> & variablePositions){
	System mySys;
	setUpSystem(_initialSystem, mySys, variablePositions);

	runPreSetUpQuench(mySys, 10);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd){
	System mySys;

	setUpSystem(_initialSystem, mySys, variablePositions, tbd);

	runPreSetUpQuench(mySys, 10, tbd);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, uint _numIterations){
	System mySys;
	setUpSystem(_initialSystem, mySys);

	runPreSetUpQuench(mySys, _numIterations);
	return mySys;
}

System Quench::runQuench(System & _initialSystem, uint _numIterations, TwoBodyDistanceDependentPotentialTable & tbd){
	System mySys;
	setUpSystem(_initialSystem, mySys, tbd);

	runPreSetUpQuench(mySys, _numIterations, tbd);
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
					double surround = calculator.calculateSurroundingEnergy(_mySystem,i,j,surroundEnergies, currentAllRotamers);
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
        OnTheFlyManager otfmanager(&_mySystem, parfile);
	return otfmanager.calculateStateEnergy(_mySystem,currentRotamers);

}

