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


#include "SelfConsistentMeanField.h"

SelfConsistentMeanField::SelfConsistentMeanField() {
	pFixed = NULL;
	pSelfE = NULL;
	pPairE = NULL;
	pBaseLines = NULL;
	setup();
}

SelfConsistentMeanField::SelfConsistentMeanField(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setup();
	setEnergyTables(NULL, &_selfEnergies, &_pairEnergies, NULL);
}

SelfConsistentMeanField::SelfConsistentMeanField(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setup();
	setEnergyTables(&_fixEnergy, &_selfEnergies, &_pairEnergies, NULL);
}

SelfConsistentMeanField::SelfConsistentMeanField(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setup();
	setEnergyTables(NULL, &_selfEnergies, &_pairEnergies, &_baselines);
}

SelfConsistentMeanField::SelfConsistentMeanField(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setup();
	setEnergyTables(&_fixEnergy, &_selfEnergies, &_pairEnergies, &_baselines);
}

SelfConsistentMeanField::~SelfConsistentMeanField() {
	deletePointers();
}

void SelfConsistentMeanField::setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setEnergyTables(NULL, &_selfEnergies, &_pairEnergies, NULL);
}

void SelfConsistentMeanField::setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setEnergyTables(NULL, &_selfEnergies, &_pairEnergies, &_baselines);
}

void SelfConsistentMeanField::setEnergyTables(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setEnergyTables(&_fixEnergy, &_selfEnergies, &_pairEnergies, NULL);
}

void SelfConsistentMeanField::setEnergyTables(double & _fixEnergy, vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setEnergyTables(&_fixEnergy, &_selfEnergies, &_pairEnergies, &_baselines);
}

void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines) {
	
	/*******************************************************
	 *  CHECK THE DIMENSIONS OF THE TABLES FIRST
	 *******************************************************/
	
	// self energy table should not be emtpy
	if (_pSelfEnergies->size() == 0) {
		cerr << "ERROR 7146: the self energy table has zero size void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(7146);
	}

	// pair and self sizes should match
	if (_pSelfEnergies->size() != _pPairEnergies->size()) {
		cerr << "ERROR 7149: the self energy table (" << _pSelfEnergies->size() << ") has different size ththan the pair energy table (" << _pPairEnergies->size() << ") at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(7149);
		for (int i=0; i<_pPairEnergies->size(); i++) {
			// pair and self rotamers sizes should match
			if ((*_pSelfEnergies)[i].size() != (*_pPairEnergies)[i].size()) {
				cerr << "ERROR 7152: at position " << i << " the pair energy table (" << (*_pPairEnergies)[i].size() << ") has different size than the self energy table (" << (*_pSelfEnergies)[i].size() << ") at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
				exit(7152);
			}
			for (int ir=0; ir<(*_pPairEnergies)[i].size(); ir++) {
				// pos i rot ir should have an entry for every prior pos
				if ((*_pPairEnergies)[i][ir].size() != i) {
					cerr << "ERROR 7155: at position " << i << ", rotamer " << ir << ", unexpected size (" << (*_pPairEnergies)[i][ir].size() << " != " << i << ")  at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
					exit(7155);
				}
				for (int j=0; j<(*_pPairEnergies)[i][ir].size(); j++) {
					// at pos j is got to be the same number or rotamers than the self energies at j
					if ((*_pSelfEnergies)[j].size() != (*_pPairEnergies)[i][ir][j].size()) {
						cerr << "ERROR 7158: at position " << i << ", rotamer " << ir << ", second position " << j << ", the pair energy table (" << (*_pPairEnergies)[i][ir][j].size() << ") has different size than the self energy table (" << (*_pSelfEnergies)[j].size() << ") at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
						exit(7158);
					}
				}
			}
		}
	}
	// baseline and self sizes should match
	if (_pBaselines != NULL) {
	       	if (_pBaselines->size() != _pSelfEnergies->size()) {
			cerr << "ERROR 7161: the baseline table (" << _pBaselines->size() << ") has different size than the selfEnergy table (" << _pSelfEnergies->size() << ") at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
			exit(7161);
		}
		for (int i=0; i<_pSelfEnergies->size(); i++) {
			if ((*_pSelfEnergies)[i].size() != (*_pBaselines)[i].size()) {
				cerr << "ERROR 7164: at position " << i << " the baseline table (" << (*_pBaselines)[i].size() << ") has different size than the selfEnergy table (" << (*_pSelfEnergies)[i].size() << ") at void SelfConsistentMeanField::setEnergyTables(double * _pFixEnergy, vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines)" << endl;
				exit(7164);
			}
		}
	}

	/*******************************************************
	 *  ASSIGN THE POINTERS
	 *******************************************************/
	pSelfE = _pSelfEnergies;
	pPairE = _pPairEnergies;
	pFixed = _pFixEnergy;
	pBaseLines = _pBaselines;
	initialize();
}

void SelfConsistentMeanField::setup() {
	T = 298.0; 
	RT = MslTools::R * T;
	lambda = 0.9;
	sumSquare = 0.0;
	cycleCounter = 0;

	deleteRng = true;
	verbose = false;
	pRng = new RandomNumberGenerator;
	
}

void SelfConsistentMeanField::deletePointers() {
	if (deleteRng == true) {
		delete pRng;
	}
}

void SelfConsistentMeanField::setT(double _T) {
	T = _T;
	RT = MslTools::R * T;
}

double SelfConsistentMeanField::getT() const {
	return T;
}

void SelfConsistentMeanField::initialize() {
	initializeMask();
	initializeProbs();
}

void SelfConsistentMeanField::initializeProbs() {
	p.clear();
	selfConsE.clear();
	cycleCounter = 0;
	currentState.clear();
	//cout << "UUUQ pSelfE size " << pSelfE->size() << endl;
	for (int i=0; i<pSelfE->size(); i++) {
		//cout << "UUUQ pSelfE[" << i << "] size " << (*pSelfE)[i].size() << endl;
		// count the number of alive rotamers
		int aliveRots = 0;
		for (unsigned int j=0; j<mask[i].size(); j++) {
			if (mask[i][j]) {
				aliveRots++;
			}
		}
		if (aliveRots == 0) {
			cerr << "ERROR 54912: no available state in mask in void SelfConsistentMeanField::initialize()" << endl;
			exit(54912);
		}
		//int aliveRots = (*pSelfE)[i].size();
		// assign a default 1/num_of_rots probability to everything (including the dead ones)
		p.push_back(vector<double>((*pSelfE)[i].size(), (double)1/(double)aliveRots));
		// zero the dead rots
		for (unsigned int j=0; j<mask[i].size(); j++) {
			//cout << "UUUQ " << i << endl;
			if (!mask[i][j]) {
				p[i][j] = 0.0;
			}
			//cout << "  UUUQ " << i << " " << j << " " << p[i][j] << endl;
		}
		selfConsE.push_back(vector<double>((*pSelfE)[i].size(), 0.0));
		currentState.push_back(-1);
	}
	//for (unsigned int i=0; i<p.size(); i++) {
	//	for (unsigned int j=0; j<p[i].size(); j++) {
	//	}
	//}
	//initializeMask();
}

void SelfConsistentMeanField::initializeMask() {
	mask.clear();
	for (int i=0; i<pSelfE->size(); i++) {
		int rots = (*pSelfE)[i].size();
		mask.push_back(vector<bool>(rots, true));
	}
}

void SelfConsistentMeanField::setMask(vector<vector<bool> > _mask) {
	if (pSelfE->size() != _mask.size()) {
		// ERROR
	}
	for (unsigned int i=0; i<pSelfE->size(); i++) {
		if ((*pSelfE)[i].size() != _mask[i].size()) {
			// ERROR
		}
	}
	mask = _mask;
	initializeProbs();
}

vector<vector<bool> > SelfConsistentMeanField::getMask() const {
	return mask;
}


void SelfConsistentMeanField::readMaskFromFile(string _filename) {

	ifstream maskFile_fs(_filename.c_str(), ios::in | ios::binary);
	int dimension = 0;
	maskFile_fs.read((char*)&dimension, sizeof(int));
	if (maskFile_fs.fail()) {
		cerr << "ERROR 2842: unable to read mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
		exit(2842);
	}
	int expectedSize = sizeof(int);
	if (dimension != 2) {
		cerr << "ERROR 2841: malformed mask table in file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
		exit(2841);
	}
	int size = 0;
	maskFile_fs.read((char*)&size, sizeof(int));
	if (maskFile_fs.fail()) {
		cerr << "ERROR 2846: unable to read mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
		exit(2846);
	}
	expectedSize += sizeof(int);
	if (size != pSelfE->size()) {
		cerr << "ERROR 2849: unmatching size of mask table in file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
		exit(2849);
	}
	vector<vector<bool> > newMask;
	for (unsigned int i=0; i<size; i++) {
		int size2 = 0;
		maskFile_fs.read((char*)&size2, sizeof(int));
		if (maskFile_fs.fail()) {
			cerr << "ERROR 2853: unable to read mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
			exit(2853);
		}
		expectedSize += sizeof(int);
		if (size2 != (*pSelfE)[i].size()) {
			cerr << "ERROR 2857: unmatching size of mask table in file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
			exit(2857);
		}
		newMask.push_back(vector<bool>(size2, true));
		expectedSize += sizeof(int) * size2;
	}

	/*************************************************
	 *  Error checking: stat the file and check if the 
	 *  file size is what expected
	 **************************************************/ 
	struct stat results;
	if (stat(_filename.c_str(), &results) == 0) {
		if (expectedSize != results.st_size) {
			cerr << "ERROR 2861: unmatching size of mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
			exit(2861);
		}
	} else {
		cerr << "ERROR 2865: cannot stat mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
		exit(2865);
	}

	for (unsigned int i=0; i<newMask.size(); i++) {
		for (unsigned int j=0; j<newMask[i].size(); j++) {
			int flag = 0;
			maskFile_fs.read((char*)&flag, sizeof(int));
			if (maskFile_fs.fail()) {
				cerr << "ERROR 2868: unable to read mask table file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
				exit(2868);
			}
			if (flag == 0) {
				newMask[i][j] = false;
			} else if (flag != 1) {
				newMask[i][j] = true;
			} else {
				cerr << "ERROR 2872: malformed mask table in file " << _filename << " in void SelfConsistentMeanField::readMaskFromFile(string _filename)" << endl;
				exit(2872);
			}
		}
	}
	
	initialize();
	mask = newMask;
	
}

void SelfConsistentMeanField::cycle() {
	cycleCounter++;

	for (unsigned int i=0; i<p.size(); i++) {
		for (unsigned int j=0; j<p[i].size(); j++) {
		}
	}
	/************************************************
	 *  Calculate the new average energy of the rotames, based
	 *  on the current probabilities
	 ************************************************/
	for (int i=0; i<pSelfE->size(); i++) {
		for (int ir=0; ir<(*pSelfE)[i].size(); ir++) {
			if (mask[i][ir]) {
				if (pBaseLines != NULL) {
					selfConsE[i][ir] = (*pBaseLines)[i][ir] + (*pSelfE)[i][ir];
				} else {
					selfConsE[i][ir] = (*pSelfE)[i][ir];
				}
			} else {
				selfConsE[i][ir] = 1e+100;
			}
		}
	}
	for (int i=0; i<pPairE->size(); i++) {
		for (int j=0; j<i; j++) {
			for (int ir=0; ir<(*pPairE)[i].size(); ir++) {
				if (mask[i][ir]) {
					for (int jr=0; jr<(*pPairE)[j].size(); jr++) {
						if (mask[j][jr]) {
							selfConsE[i][ir] += p[j][jr] * (*pPairE)[i][ir][j][jr];
							selfConsE[j][jr] += p[i][ir] * (*pPairE)[i][ir][j][jr];
						}
					}
				}
			}
		}
	}

	/************************************************
	 *  Recalculate the probabilities based on the
	 *  new energies
	 ************************************************/
	pPrev = p;
	vector<vector<double> > subtractedE(selfConsE.size(), vector<double>());
	vector<double> norm(selfConsE.size(), 0.0);
	for (int i=0; i<selfConsE.size(); i++) {
		double min = 0.0;
		bool foundFirst = false;
		for (int ir=0; ir<selfConsE[i].size(); ir++) {
			if (mask[i][ir]) {
				if (!foundFirst || selfConsE[i][ir] < min) {
					min = selfConsE[i][ir];
					foundFirst = true;
				}
			}
		}
		for (int ir=0; ir<selfConsE[i].size(); ir++) {
			if (mask[i][ir]) {
				subtractedE[i].push_back(selfConsE[i][ir]-min);
				norm[i] += exp(-(subtractedE[i][ir]/RT));
			} else {
				subtractedE[i].push_back(1e+100);
			}
		}
	}
	
	sumSquare = 0.0;
	for (int i=0; i<selfConsE.size(); i++) {
		for (int ir=0; ir<selfConsE[i].size(); ir++) {
			if (mask[i][ir]) {
				p[i][ir] = (exp(-(subtractedE[i][ir]/RT)))/norm[i];
				// lambda introduces some memory, this is needed to 
				// avoid cyclic obsillations
				if (lambda < 1.0) {
					p[i][ir] = lambda * p[i][ir] + (1.0 - lambda) * pPrev[i][ir];
					sumSquare += pow((p[i][ir] - pPrev[i][ir]),2);
				}
			} else {
				p[i][ir] = 0;
			}
		}
	}
}

void SelfConsistentMeanField::setLambda(double theLambda) {
	if (theLambda < 0.0 || theLambda > 1.0) {
		cerr << "ERROR 7183: lambda must be between 0.0 and 1.0 in void SelfConsistentMeanField::setLambda(double theLambda)" << endl;
		exit(7183);
	}
	lambda = theLambda;
}


double SelfConsistentMeanField::getLambda() const {
	return lambda;
}

vector<vector<double> > & SelfConsistentMeanField::getP() {
	return p;
}

vector<vector<double> > & SelfConsistentMeanField::getE() {
	return selfConsE;
}

double SelfConsistentMeanField::getStateP(vector<unsigned int> states) const {
	double out = 1.0;
	if (states.size() != p.size()) {
		cerr << "ERROR 7194: the size of the vector of the states is not equal to the number of positions in double SelfConsistentMeanField::getStateP(vector<int> states) const" << endl;
		exit(7194);
	}

	for (int i=0; i<states.size(); i++) {
		if (states[i] < 0 || states[i] >= p[i].size()) {
			cerr << "ERROR 7197: the state (" << states[i] << ") at the position index " << i << " is out of the range of the number of rotamers (0-" << p[i].size() << ") in double SelfConsistentMeanField::getStateP(vector<int> states) const" << endl;
			exit(7197);
		}
		out *= p[i][states[i]];
	}
	return out;
}

double SelfConsistentMeanField::getPvariation() const {
	return pow(sumSquare, 0.5);
}

double SelfConsistentMeanField::getNumberOfCycles() const {
	return cycleCounter;
}

vector<unsigned int> SelfConsistentMeanField::getMostProbableState() {
	currentState.clear();
	for (int i=0; i<p.size(); i++) {
		int rot = -1;
		double max = 0.0;
		for (int j=0; j<p[i].size(); j++) {
			if (j==0 || p[i][j] > max) {
				rot = j;
				max = p[i][j];
			}
		}
		currentState.push_back(rot);
	}
	return currentState;
			
}

vector<unsigned int> SelfConsistentMeanField::getRandomState() {
	currentState.clear();
	for (int i=0; i<p.size(); i++) {
		currentState.push_back(selectRandomStateAtPosition(i));
	}
	return currentState;
}

vector<unsigned int> SelfConsistentMeanField::moveRandomState()  {
	return moveRandomState(1);
}

vector<unsigned int> SelfConsistentMeanField::moveRandomState(unsigned int _numOfMoves)  {
	/****************************************************
	 *  this function selects one or more moves proportionally
	 *  to the SCMF probabilities
	 ****************************************************/

	if (_numOfMoves > currentState.size()) {
		cerr << "WARNING 3205: the number of step is larger than the number of positions in vector<int> SelfConsistentMeanField::moveRandomState(unsigned int _numOfMoves)" << endl;
		_numOfMoves = currentState.size();
	}
	
	/*********************************************
	 *  Pick a position, with a bias for those that have
	 *  more probable alternative states than the current
	 *********************************************/
	vector<bool> alreadySelected(currentState.size(), false);
	for (unsigned int i=0; i<_numOfMoves; i++) {
		vector<double> residualP;
		double sumP = 0.0;
		for (unsigned int j=0; j<currentState.size(); j++) {
			if (alreadySelected[j]) {
				// if the position was already choosen, no prob to it
				residualP.push_back(0);
			} else {
				if (currentState[j] != -1) {
					// if there is a state, the p is 1 minus the p of the state
					residualP.push_back(1 - p[j][currentState[j]]);
				} else {
					// if there is no state each position has the same prob
					residualP.push_back(1);
				}
			}
			sumP += residualP.back();
		}


		if (sumP == 0.0) {
			// there is no available move, try again
			continue;
		}
	
		pRng->setDiscreteProb(residualP);
		unsigned int randomPos = pRng->getRandomDiscreteIndex();
		/*
		for (unsigned int i=0; i<residualP.size(); i++) {
			cout << " UUU " << i << " (" << currentState[i] << ") " << residualP[i] << endl;
		}
		cout << "UUU selected " << randomPos << endl;
		*/
		currentState[randomPos] = selectRandomStateAtPosition(randomPos);
		alreadySelected[randomPos] = true;
	}

	return currentState;
}

int SelfConsistentMeanField::selectRandomStateAtPosition(int _position) const {
	if (_position >= p.size()) {
		cerr << "ERROR 7210: position " << _position << " out of range in int SelfConsistentMeanField::selectRandomStateAtPosition(int _position) const" << endl;
		exit(7210);
	}

	if (p[_position].size() == 1) {
		return 0;
	}

	vector<double> residualP;
	double sumP = 0.0;
	for (int i=0; i<p[_position].size(); i++) {
		if (i == currentState[_position] || !mask[_position][i]) {
			residualP.push_back(0.0);
		} else {
			residualP.push_back(p[_position][i]);
			sumP = p[_position][i];
		}
	}

	if (sumP == 0.0) {
		// no move avalable, stay there
		return currentState[_position];
	}

	pRng->setDiscreteProb(residualP);
	return pRng->getRandomDiscreteIndex();
}


void SelfConsistentMeanField::setCurrentState(vector<unsigned int> _state) {
	if (_state.size() != p.size()) {
		cerr << "ERROR 7223: the size of the vector of the states is not equal to the number of positions in void SelfConsistentMeanField::setCurrentState(vector<int> _state)" << endl;
		exit(7223);
	}

	for (int i=0; i<_state.size(); i++) {
		if (_state[i] < 0 || _state[i] >= p[i].size()) {
			cerr << "ERROR 7227: the state (" << _state[i] << ") at the position index " << i << " is out of the range of the number of rotamers (0-" << p[i].size() << ") in void SelfConsistentMeanField::setCurrentState(vector<int> _state)" << endl;
			exit(7227);
		}
	}

	currentState = _state;
}

double SelfConsistentMeanField::getStateEnergy(vector<unsigned int> _state) const {
	double energy = 0.0;
	if (pFixed != NULL) {
		energy += *pFixed;
	}
	
	for (int i=0; i<pSelfE->size(); i++) {
		if (pBaseLines != NULL) {
			energy += (*pBaseLines)[i][_state[i]] + (*pSelfE)[i][_state[i]];
		} else {
			energy += (*pSelfE)[i][_state[i]];
		}
	}
			
	for (int i=0; i<pPairE->size(); i++) {
		for (int j=0; j<i; j++) {
			energy += (*pPairE)[i][_state[i]][j][_state[j]];
		}
	}
	return energy;
}

double SelfConsistentMeanField::getAverageEnergy() {
	double energy = 0.0;
	if (pFixed != NULL) {
		energy += *pFixed;
	}
	
	for (int i=0; i<pSelfE->size(); i++) {
		for (int ir=0; ir<(*pSelfE)[i].size(); ir++) {
			if (pBaseLines != NULL) {
				energy += p[i][ir] * ((*pBaseLines)[i][ir] + (*pSelfE)[i][ir]);
			} else {
				energy += p[i][ir] * (*pSelfE)[i][ir];
			}
		}
	}
			
	for (int i=0; i<pPairE->size(); i++) {
		for (int ir=0; ir<(*pPairE)[i].size(); ir++) {
			for (int j=0; j<i; j++) {
				for (int jr=0; jr<(*pPairE)[j].size(); jr++) {
					energy += p[i][ir] * p[j][jr] * (*pPairE)[i][ir][j][jr];
				}
			}
		}
	}
	return energy;
}

vector<unsigned int> SelfConsistentMeanField::runMC(double _startingTemperature, double _endingTemperature, int _scheduleCycles, int _scheduleShape, int _maxRejectionsNumber, int _convergedSteps, double _convergedE) {
	/******************************************************************************
	 *              === SCMF Biased MONTE CARLO OPTIMIZATION ===
	 *   Based on the paper "Using aself-consistent fields to bias Monte Carlo methods
	 *   with applications to designing and sampling protein sequences" - Jinming Zou
	 *   and Jeffery Saven, Journal of chemical physics 2003.
	 ******************************************************************************/
	time_t startMCOtime, endMCOtime;
	double MCOTime;

	time (&startMCOtime);

	if (verbose) {
		cout << "===================================" << endl;
		cout << "Run SCMF Biased Monte Carlo Optimization" << endl;
		cout << endl;
	}
	MonteCarloManager MCMngr(_startingTemperature, _endingTemperature,_scheduleCycles, _scheduleShape, _maxRejectionsNumber, _convergedSteps, _convergedE);
	MCMngr.setRandomNumberGenerator(pRng);
	
	unsigned int cycleCounter = 0;
	unsigned int moveCounter = 0;

	// We wish to start from the most probable state
	vector<unsigned int> bestState = getMostProbableState();
	double bestEnergy = getStateEnergy(bestState);
	double prevStateP = getStateP(bestState);
	double stateP = 0.0;

	vector<unsigned int> prevStateVec = bestState;
	vector<unsigned int> stateVec;

	currentState = bestState; // should already be true .. BUT just in case
	MCMngr.setEner(bestEnergy);

	while (!MCMngr.getComplete()) {
		// make atleast 1 move and update the current state
		stateVec = moveRandomState(pRng->getRandomInt(pSelfE->size() - 1) + 1); 
		//stateVec = moveRandomState(); 
		stateP = getStateP(stateVec);
		if (verbose) {
			for (int i = 0; i < stateVec.size(); i++){
				cout << stateVec[i] << ",";
			}
			cout << endl;
		}
		double oligomerEnergy = getStateEnergy(stateVec);

		if (verbose) {
			cout << "SCMFBMC [" << cycleCounter << "]: ";
		}

		if(oligomerEnergy < bestEnergy) {
			bestEnergy = oligomerEnergy;
			bestState = stateVec;
		}
		if (!MCMngr.accept(oligomerEnergy, stateP/prevStateP)) {
			setCurrentState(prevStateVec);
			if (verbose) {
				cout << "SCMFBMC: State not accepted, E=" << oligomerEnergy << "\n";
			}
		} else {
			prevStateVec = stateVec;
			prevStateP = stateP;
			if (verbose) {
				cout << "SCMFBMC: State accepted, E=" << oligomerEnergy << "\n";
			}
		}
		cycleCounter++;
		moveCounter++;
	}

	time (&endMCOtime);
	MCOTime = difftime (endMCOtime, startMCOtime);
	if (verbose) {
		cout << endl;
		cout << "Best SCMF Biased MC Energy: " << bestEnergy << endl;
		cout << "SCMF Biased MCO Time: " << MCOTime << " seconds" << endl;
		cout << "===================================" << endl;
	}
	return bestState;
}

