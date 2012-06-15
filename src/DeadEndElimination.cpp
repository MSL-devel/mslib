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

#include "DeadEndElimination.h"

using namespace MSL;
using namespace std;


DeadEndElimination::DeadEndElimination() {
	setInitialVariables();
}

DeadEndElimination::DeadEndElimination(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setInitialVariables();
	setEnergyTables(&_selfEnergies, &_pairEnergies, NULL);
}

DeadEndElimination::DeadEndElimination(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setInitialVariables();
	setEnergyTables(&_selfEnergies, &_pairEnergies, &_baselines);
}

void DeadEndElimination::setInitialVariables() {
	selfEnergy = NULL;
	pairEnergy = NULL;
	pBaseLines = NULL;
	verboseLevel = 1;
	setVerbose(true, 1); // default low level verbose mode
	cycles = 0;
	eliminatedCounter = 0;
	flaggedCounter = 0;
	enerOffset = 0.0;
	afterPair_flag = false;
	responsibleForEnergyTableMemory = false;
	totalNumPositions = 0;
	totalNumRotamers  = 0;
}

DeadEndElimination::~DeadEndElimination() {
}

void DeadEndElimination::setVerbose(bool _flag) {
	setVerbose(_flag, verboseLevel);
}
	
void DeadEndElimination::setVerbose(bool _flag, unsigned int _level) {

	/*************************************
	 *  Verbose levels:
	 *
	 *  0 not verbose
	 *  1 very little: only print summary at the end of whore run (like Simple Goldstein Single)
	 *  2 some: print summary at every cycle
	 *  3 very verbose: print every elimination
	 *************************************/
	verboseLevel = _level;
	if (_level == 0 || _flag == false) {
		verboseLevel1_flag = false;
		verboseLevel2_flag = false;
		verboseLevel3_flag = false;
	} else if (_level == 1) {
		verboseLevel1_flag = true;
		verboseLevel2_flag = false;
		verboseLevel3_flag = false;
	} else if (_level == 2) {
		verboseLevel1_flag = true;
		verboseLevel2_flag = true;
		verboseLevel3_flag = false;
	} else if (_level == 3) {
		verboseLevel1_flag = true;
		verboseLevel2_flag = true;
		verboseLevel3_flag = true;
	} else {
		cerr << "ERROR 3808: invalid verbose level in void DeadEndElimination::setVerbose(bool _flag, unsigned int _level)" << endl;
		exit(3808);
	}}

void DeadEndElimination::setEnergyOffset(double _offset) {
	if (_offset < 0.0) {
		cerr << "ERROR 3810: offset (" << _offset << ") cannot be negative in void DeadEndElimination::setEnergyOffset(double _offset)" << endl;
		exit(3810);
	}
	enerOffset = _offset;
}

double DeadEndElimination::getEnergyOffset() const {
	return enerOffset;
}

bool DeadEndElimination::getVerbose() const {
	return verboseLevel1_flag;
}

unsigned int DeadEndElimination::getVerboseLevel() const {
	return verboseLevel;
}

unsigned int DeadEndElimination::getEliminatedCounter() const {
	return eliminatedCounter;
}

unsigned int DeadEndElimination::getFlaggedCounter() const {
	return flaggedCounter;
}

unsigned int DeadEndElimination::getRunCycles() const {
	return cycles;
}

double DeadEndElimination::getTotalCombinations() const {
	double total = 1;
	for (unsigned int i=0; i<alive.size(); i++) {
		unsigned int count = 0;
		for (unsigned int j=0; j<alive[i].size(); j++) {
			if (alive[i][j]) {
				count++;
			}
		}
		total *= count;
	}
	return total;
}

bool DeadEndElimination::isAlive(vector<int> _states) const {
	if (_states.size() != alive.size()) {
		cerr << "ERROR 3812: the number of positions is incorrect in bool DeadEndElimination::isAlive(vector<int> _states) const" << endl;
		exit(3812);
	}

	bool out = true;
	for (unsigned int i=0; i<_states.size(); i++) {
		if (_states[i] < 0 || _states[i] > alive[i].size()) {
			cerr << "ERROR 3816: state out of range at position " << i << " in bool DeadEndElimination::isAlive(vector<int> _states) const" << endl;
			exit(3816);
		}
		if (!alive[i][_states[i]]) {
			return false;
		}
	}
	return out;
}

void DeadEndElimination::setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies) {
	setEnergyTables(&_selfEnergies, &_pairEnergies, NULL);
}

void DeadEndElimination::setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines) {
	setEnergyTables(&_selfEnergies, &_pairEnergies, &_baselines);
}

void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies) {
	setEnergyTables(_selfEnergynergies, _pairEnergynergies, NULL);
}

void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines) {
	

	if (_selfEnergynergies == NULL) {
		cerr << "ERROR 3818: null pointer for self energy table in void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(3818);
	}
	
	if (_pairEnergynergies == NULL) {
		cerr << "ERROR 3819: null pointer for pair energy table in void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(3819);
	}
	
	/*******************************************************
	 *  CHECK THE DIMENSIONS OF THE TABLES FIRST
	 *******************************************************/
	
	// self energy table should not be emtpy
	if (_selfEnergynergies->size() == 0) {
		cerr << "ERROR 3820: the self energy table has zero size in void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(3820);
	}

	// pair and self sizes should match
	if (_selfEnergynergies->size() != _pairEnergynergies->size()) {
		cerr << "ERROR 3824: the self energy table (" << _selfEnergynergies->size() << ") has different size ththan the pair energy table (" << _pairEnergynergies->size() << ") at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
		exit(3824);


		/*	
			THIS CODE IS NEVER REACHED...A BUG?
		 */
		for (int i=0; i<_pairEnergynergies->size(); i++) {
			// pair and self rotamers sizes should match
			if ((*_selfEnergynergies)[i].size() != (*_pairEnergynergies)[i].size()) {
				cerr << "ERROR 3828: at position " << i << " the pair energy table (" << (*_pairEnergynergies)[i].size() << ") has different size than the self energy table (" << (*_selfEnergynergies)[i].size() << ") at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
				exit(3828);
			}
			for (int ir=0; ir<(*_pairEnergynergies)[i].size(); ir++) {
				// pos i rot ir should have an entry for every prior pos
				if ((*_pairEnergynergies)[i][ir].size() != i) {
					cerr << "ERROR 3832: at position " << i << ", rotamer " << ir << ", unexpected size (" << (*_pairEnergynergies)[i][ir].size() << " != " << i << ")  at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
					exit(3832);
				}
				for (int j=0; j<(*_pairEnergynergies)[i][ir].size(); j++) {
					// at pos j is got to be the same number or rotamers than the self energies at j
					if ((*_selfEnergynergies)[j].size() != (*_pairEnergynergies)[i][ir][j].size()) {
						cerr << "ERROR 3836: at position " << i << ", rotamer " << ir << ", second position " << j << ", the pair energy table (" << (*_pairEnergynergies)[i][ir][j].size() << ") has different size than the self energy table (" << (*_selfEnergynergies)[j].size() << ") at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
						exit(3836);
					}
				}
			}
		}
	}
	/*
	// baseline and self sizes should match
	if (_pBaselines != NULL) {
	       	if (_pBaselines->size() != _selfEnergynergies->size()) {
			cerr << "ERROR 3840: the baseline table (" << _pBaselines->size() << ") has different size than the selfEnergy table (" << _selfEnergynergies->size() << ") at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
			exit(3840);
		}
		for (int i=0; i<_selfEnergynergies->size(); i++) {
			if ((*_selfEnergynergies)[i].size() != (*_pBaselines)[i].size()) {
				cerr << "ERROR 3844: at position " << i << " the baseline table (" << (*_pBaselines)[i].size() << ") has different size than the selfEnergy table (" << (*_selfEnergynergies)[i].size() << ") at void DeadEndElimination::setEnergyTables(vector<vector<double> > * _selfEnergynergies, vector<vector<vector<vector<double> > > > * _pairEnergynergies, vector<vector<double> > * _pBaselines)" << endl;
				exit(3844);
			}
		}
	}
	*/

	/*******************************************************
	 *  ASSIGN THE POINTERS
	 *******************************************************/
	selfEnergy = _selfEnergynergies;
	pairEnergy = _pairEnergynergies;
	if (_pBaselines != NULL) {
		setBaselines(_pBaselines);
	}
	initializeMask();
}

void DeadEndElimination::setBaselines(vector<vector<double> > & _baselines) {
	setBaselines(&_baselines);
}

void DeadEndElimination::setBaselines(vector<vector<double> > * _pBaselines) {
	if (_pBaselines == NULL) {
		cerr << "ERROR 3848: null pointer for baselines in void DeadEndElimination::setBaselines(vector<vector<double> > * _pBaselines)" << endl;
		exit(3848);
	}

	// baseline and self sizes should match
	if (_pBaselines->size() != selfEnergy->size()) {
		cerr << "ERROR 3840: the baseline table (" << _pBaselines->size() << ") has different size than the selfEnergy table (" << selfEnergy->size() << ") at void DeadEndElimination::setBaselines(vector<vector<double> > * _pBaselines)" << endl;
		exit(3840);
	}
	for (int i=0; i<selfEnergy->size(); i++) {
		if ((*selfEnergy)[i].size() != (*_pBaselines)[i].size()) {
			cerr << "ERROR 3844: at position " << i << " the baseline table (" << (*_pBaselines)[i].size() << ") has different size than the selfEnergy table (" << (*selfEnergy)[i].size() << ") at void DeadEndElimination::setBaselines(vector<vector<double> > * _pBaselines)" << endl;
			exit(3844);
		}
	}

	/*******************************************************
	 *  ASSIGN THE POINTER
	 *******************************************************/
	pBaseLines = _pBaselines;
}

void DeadEndElimination::initializeMask() {
	alive.clear();
	for (unsigned int ip=0; ip<selfEnergy->size(); ip++) {
		unsigned int ir = (*selfEnergy)[ip].size();
		alive.push_back(vector<bool>(ir, true));
	}
	flaggedPair.clear();
	for (unsigned int ip=0; ip<pairEnergy->size(); ip++) {
		flaggedPair.push_back(vector<vector<vector<bool> > >());
		for (unsigned int ir=0; ir<(*pairEnergy)[ip].size(); ir++) {
			flaggedPair[ip].push_back(vector<vector<bool> >());
			for (unsigned int jp=0; jp<(*pairEnergy)[ip][ir].size(); jp++) {
				flaggedPair[ip][ir].push_back(vector<bool>());
				for (unsigned int jr=0; jr<(*pairEnergy)[ip][ir][jp].size(); jr++) {
					flaggedPair[ip][ir][jp].push_back(false);
				}
			}
		}
	}
}

void DeadEndElimination::setMask(vector<vector<bool> > theMask) {
	alive = theMask;
}

vector<vector<bool> > DeadEndElimination::getMask() const {
	return alive;
}

bool DeadEndElimination::runSimpleGoldsteinSingles() {
	if (verboseLevel1_flag) {
		cout << "Starting Simple Goldstein Singles (SGS)" << endl;
	}
	cycles = 0;
	eliminatedCounter = 0;
	unsigned int prevCounter = 0;
	while(true) {
		bool eliminatedSomething = false;
		// for each position...
		for (unsigned int posI=0; posI<pairEnergy->size(); posI++) {
			// ... for each rotamer at the position...
			for (unsigned int rotR=0; rotR<(*pairEnergy)[posI].size(); rotR++) {
				if (!alive[posI][rotR]) {
					// eliminated already
					continue;
				}
				for (unsigned int rotT=rotR+1; rotT<(*pairEnergy)[posI].size(); rotT++) {
					// ... for each other rotamer at the position...
					if (!alive[posI][rotT]) {
						continue;
					}
					// ... run the simple Goldstein single.
					// if eliminateSomething was false, 
					// set it to true if the function returns true
					if (simpleGoldsteinSinglesIteration(posI, rotR, rotT)) {
						eliminatedSomething = true;
						eliminatedCounter++;
					}
				}
			}
		}
		if (verboseLevel2_flag) {
			cout << "DEE SGS " << cycles << ": eliminated " << eliminatedCounter - prevCounter << endl;
			prevCounter = eliminatedCounter;
		}
		if (!eliminatedSomething) {
			if (verboseLevel1_flag) {
				cout << "Ended Simple Goldstein Singles (SGS), eliminated " << eliminatedCounter << " in " << cycles+1 << " cycles" << endl;
			}
			// end if we did not eliminate anything
			if (cycles == 0) {
				// the Goldstein Single run did not eliminate anything
				return false;
			} else {
				// the Goldstein Single run eliminated something
				return true;
			}
		}
		cycles++;
	}
		
}

bool DeadEndElimination::simpleGoldsteinSinglesIteration(unsigned int _posI, unsigned int _rotR, unsigned int _rotT) {

	/********************************************
	 *  At position I, rotamer It eliminates rotamer Ir
	 *  if 
	 *  E(Ir) - E(It) + SUM(  min[E(IrJu)-E(ItJu)] ) > 0
	 *                 J,J!=I  u
	 *  
	 ********************************************/
	//cout << "UUU DEE 0 " << alive.size() << " " << _posI << endl;
	//cout << "UUU DEE 0 " << alive[_posI].size() << " " << _rotR << " " << _rotT << endl;
	if (!alive[_posI][_rotR] || !alive[_posI][_rotT]) {
		//cout << "UUU DEE 1" << endl;
		return false;
	}
	//cout << "UUU DEE 2" << endl;
	// add the self energies
	double Er = (*selfEnergy)[_posI][_rotR] - (*selfEnergy)[_posI][_rotT];
	//cout << "UUU DEE 3" << endl;
	// possibly add the baselines
	//cout << "UUU DEE 4" << endl;
	if (pBaseLines != NULL) {
		//cout << "UUU DEE 5" << endl;
		//cout << "UUU DEE (*pBaseLines).size() = " << (*pBaseLines).size() << " vs " << _posI << endl;
		//cout << "UUU DEE (*pBaseLines)[_posI].size() = " << (*pBaseLines)[_posI].size() << " vs " << _rotR << " " << _rotT << endl;
		Er += (*pBaseLines)[_posI][_rotR] - (*pBaseLines)[_posI][_rotT];
		//cout << "UUU DEE 6" << endl;
	}
	//cout << "UUU DEE 7" << endl;
	double Et = -Er;
	//cout << "UUU DEE 8" << endl;
	
	if (afterPair_flag) {	
	
		// for all other positions
		for (unsigned int posJ=0; posJ<pairEnergy->size(); posJ++) {
			double min = 0.0;
			double max = 0.0;
			//bool found = false;
			if (posJ==_posI) {
				continue;
			}
			bool eliminateR = true;
			bool eliminateT = true;
			minDiffIrItJuSingleAfterPair(_posI, _rotR, _rotT, posJ, min, max, eliminateR, eliminateT);
			if (eliminateR || eliminateT) {
				if (eliminateR) {
					alive[_posI][_rotR] = false;
					if (verboseLevel3_flag) {
						cout << "DEE SGS: pos/rot " << _posI << "/" << _rotR << " eliminated because there are no unflagged pairs with pos " << posJ << endl;
					}
				}
				if (eliminateT) {
					alive[_posI][_rotT] = false;
					if (verboseLevel3_flag) {
						cout << "DEE SGS: pos/rot " << _posI << "/" << _rotT << " eliminated because there are no unflagged pairs with pos " << posJ << endl;
					}
				}
				return true;
			}
			Er += min;
			Et -= max;
		}
	} else {

		// for all other positions
		for (unsigned int posJ=0; posJ<pairEnergy->size(); posJ++) {
			double min = 0.0;
			double max = 0.0;
			//bool found = false;
			if (posJ==_posI) {
				continue;
			}
			minDiffIrItJuSingle(_posI, _rotR, _rotT, posJ, min, max);
			Er += min;
			Et -= max;
		}
	}

	// the normal criteria is 0.  if enerOffset is above, more solutions are
	// saved, those that are not worse than the global minimum by more than
	// enerOffset
	if (Er > enerOffset) {
		// eliminate rotamer R
		alive[_posI][_rotR] = false;
		if (verboseLevel3_flag) {
			cout << "DEE SGS: pos/rot " << _posI << "/" << _rotR << " eliminated by " << _posI << "/" << _rotT << endl;
		}
		return true;  // we eliminated rot _posI/_rotR
	} else if (Et > enerOffset) {
		// eliminate rotamer T
		alive[_posI][_rotT] = false;
		if (verboseLevel3_flag) {
			cout << "DEE SGS: pos/rot " << _posI << "/" << _rotT << " eliminated by " << _posI << "/" << _rotR << endl;
		}
		return true;  // we eliminated rot _posI/_rotT
	}
	return false; // no elimination
	
}

void DeadEndElimination::minDiffIrItJuSingle(unsigned int _posI, unsigned int _rotR, unsigned int _rotT, unsigned int _posJ, double & _min, double & _max) {

	/*****************************************
	 *  Returns
	 * 
	 *   min[E(IrJu)-E(ItJu)]
	 *    u
	 *****************************************/

	// reset min and max in case nothing is found
	_min = 0.0;
	_max = 0.0;
	
	bool found = false;
	// the table of pair energies (with N positions) is the half of the matrix with I = 0 -> N and J = 0 -> (I-1)
	if (_posI > _posJ) {
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posI][_rotR][_posJ].size(); rotJ++) {
			if (alive[_posJ][rotJ]) {
				double diff = (*pairEnergy)[_posI][_rotR][_posJ][rotJ] - (*pairEnergy)[_posI][_rotT][_posJ][rotJ];
				if (!found) {
					found = true;
					_min = diff;
					_max = diff;
				} else {
					if (diff < _min) {
						_min = diff;
					}
					if (diff > _max) {
						_max = diff;
					}
				}
			}
		}
	} else {
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posJ].size(); rotJ++) {
			if (alive[_posJ][rotJ]) {
				double diff = (*pairEnergy)[_posJ][rotJ][_posI][_rotR] - (*pairEnergy)[_posJ][rotJ][_posI][_rotT];
				if (!found) {
					found = true;
					_min = diff;
					_max = diff;
				} else {
					if (diff < _min) {
						_min = diff;
					}
					if (diff > _max) {
						_max = diff;
					}
				}
			}
		}
	}
}

void DeadEndElimination::minDiffIrItJuSingleAfterPair(unsigned int _posI, unsigned int _rotR, unsigned int _rotT, unsigned int _posJ, double & _min, double & _max, bool & _eliminateR, bool & _eliminateT) {

	/*****************************************
	 *  This version applies when we have done 
	 *  SGP first
	 * 
	 *  Returns
	 * 
	 *   min[E(IrJu)-E(ItJu)]
	 *    u
	 *****************************************/

	// reset min and max in case nothing is found
	_min = 0.0;
	_max = 0.0;
	
	bool found = false;
	_eliminateR = true;
	_eliminateT = true;
	// the table of pair energies (with N positions) is the half of the matrix with I = 0 -> N and J = 0 -> (I-1)
	if (_posI > _posJ) {
		// if all R (or T) pairs with J rotamers are flagged, eliminate R (or T)
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posI][_rotR][_posJ].size(); rotJ++) {
			if (!alive[_posJ][rotJ]) {
				continue;
			}
			if (!flaggedPair[_posI][_rotR][_posJ][rotJ]) {
				_eliminateR = false;
				if (!_eliminateT) {
					break;
				}
			}
			if (!flaggedPair[_posI][_rotT][_posJ][rotJ]) {
				_eliminateT = false;
				if (!_eliminateR) {
					break;
				}
			}
		}
		if (_eliminateR || _eliminateT) {
			return;
		}
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posI][_rotR][_posJ].size(); rotJ++) {
			if (alive[_posJ][rotJ]) {
				double diff = (*pairEnergy)[_posI][_rotR][_posJ][rotJ] - (*pairEnergy)[_posI][_rotT][_posJ][rotJ];
				if (!found) {
					found = true;
					_min = diff;
					_max = diff;
				} else {
					if (diff < _min) {
						_min = diff;
					}
					if (diff > _max) {
						_max = diff;
					}
				}
			}
		}
	} else {
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posJ].size(); rotJ++) {
			if (!alive[_posJ][rotJ]) {
				continue;
			}
			if (!flaggedPair[_posJ][rotJ][_posI][_rotR]) {
				_eliminateR = false;
				if (!_eliminateT) {
					break;
				}
			}
			if (!flaggedPair[_posJ][rotJ][_posI][_rotT]) {
				_eliminateT = false;
				if (!_eliminateR) {
					break;
				}
			}
		}
		if (_eliminateR || _eliminateT) {
			return;
		}
		for (unsigned int rotJ=0; rotJ<(*pairEnergy)[_posJ].size(); rotJ++) {
			if (alive[_posJ][rotJ]) {
				double diff = (*pairEnergy)[_posJ][rotJ][_posI][_rotR] - (*pairEnergy)[_posJ][rotJ][_posI][_rotT];
				if (!found) {
					found = true;
					_min = diff;
					_max = diff;
				} else {
					if (diff < _min) {
						_min = diff;
					}
					if (diff > _max) {
						_max = diff;
					}
				}
			}
		}
	}
}

bool DeadEndElimination::runSimpleGoldsteinPairs() {

	unsigned int totalFlagged = 0;
	
	if (verboseLevel1_flag) {
		cout << "Starting Simple Goldstein Pairs (SGP)" << endl;
	}
	cycles = 0;
	while(true) {
		unsigned int flagged = runSimpleGoldsteinPairsOnce();
		totalFlagged += flagged;
		if (!flagged) {
			if (verboseLevel1_flag) {
				cout << "Ended Simple Goldstein Pairs (SGS), flagged " << totalFlagged << " in " << cycles+1 << " cycles" << endl;
			}
			// end if we did not eliminate anything
			if (cycles == 0) {
				// the Goldstein Single run did not flag anything
				return false;
			} else {
				// the Goldstein Single run flagged something
				return true;
			}
		}
		cycles++;
	}

}
	
unsigned int DeadEndElimination::runSimpleGoldsteinPairsOnce() {
	afterPair_flag = true;
	// for each position pair

	unsigned int flagged = flaggedCounter; // subtract from flaggedCounter at the end to get number of pairs flagged in this iteration

	for (unsigned int posI1=0; posI1<pairEnergy->size(); posI1++) {
		for (unsigned int posI2=0; posI2<posI1; posI2++) {
			// ... for each pair of rotamers at the position...
			for (unsigned int rotR1=0; rotR1<(*pairEnergy)[posI1].size(); rotR1++) {
				if (!alive[posI1][rotR1]) {
					// eliminated already
					continue;
				}
				for (unsigned int rotR2=0; rotR2<(*pairEnergy)[posI1][rotR1][posI2].size(); rotR2++) {
					if (flaggedPair[posI1][rotR1][posI2][rotR2] || !alive[posI2][rotR2]) {
						// eliminated already
						continue;
					}
					// ... for each other pair of rotamers at the position...
					for (unsigned int rotT1=0; rotT1<(*pairEnergy)[posI1].size(); rotT1++) {
						if (!alive[posI1][rotT1]) {
							continue;
						}
						for (unsigned int rotT2=0; rotT2<(*pairEnergy)[posI1][rotR1][posI2].size(); rotT2++) {
							if (flaggedPair[posI1][rotR1][posI2][rotR2] || flaggedPair[posI1][rotT1][posI2][rotT2] || !alive[posI2][rotT2] || (rotR1 == rotT1 && rotR2 == rotT2)) {
								continue;
							}
							// ... run the simple Goldstein pair.
							// if eliminateSomething was false, 
							// set it to true if the function returns true
							if (simpleGoldsteinPairIteration(posI1, rotR1, rotT1, posI2, rotR2, rotT2)) {
								flaggedCounter++;
							}
						}
					}
				}
			}
		}
	}
	flagged = flaggedCounter - flagged;
	if (verboseLevel2_flag) {
		cout << "DEE SGP " <<  ": flagged " << flagged << endl;
	}
	return flagged;
		
}

bool DeadEndElimination::simpleGoldsteinPairIteration(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2) {
	/*
	if (_posI1 < _posI2) {
		cerr << "GAZ!" << endl;
		exit(1);
		unsigned int tmp = _posI1;
		_posI1 = _posI2;
		_posI2 = tmp;
		tmp = _rotR1;
		_rotR1 = _rotR2;
		_rotR2 = tmp;
		tmp = _rotT1;
		_rotT1 = _rotT2;
		_rotT2 = tmp;
	}
	*/
	
	/********************************************
	 *  At position I1 and I2, super-rotamer It eliminates rotamer Ir
	 *  if 
	 *  E(Ir) - E(It) + E(Ir12) - E(It12) + SUM(  min[E(IrJu)-E(ItJu)] ) > 0
	 *                           J,J!=I  u
	 *  
	 *  Ir and It are pair super-rotamers
	 *
	 *  E(Ir)   = E(I1r1) + E(I2r2)
	 *  E(It)   = E(I1t1) + E(I2t2)
	 *  E(Ir12) = E(I1r1I2r2)
	 *  E(It12) = E(I1t1I2t2)
	 *  
	 ********************************************/
	// add th self energies
	double Er = (*selfEnergy)[_posI1][_rotR1] + (*selfEnergy)[_posI2][_rotR2] - (*selfEnergy)[_posI1][_rotT1] - (*selfEnergy)[_posI2][_rotT2];
	//cout << (*selfEnergy)[_posI1][_rotR1] << " " << (*selfEnergy)[_posI2][_rotR2] << " " << (*selfEnergy)[_posI1][_rotT1] << " " << (*selfEnergy)[_posI2][_rotT2] << endl;

	// possibly add the baselines
	if (pBaseLines != NULL) {
		Er += (*pBaseLines)[_posI1][_rotR1] - (*pBaseLines)[_posI1][_rotT1] + (*pBaseLines)[_posI2][_rotR2] - (*pBaseLines)[_posI2][_rotT2];
	}

	// add the pair interactions of the R and T pairs
	Er += (*pairEnergy)[_posI1][_rotR1][_posI2][_rotR2] - (*pairEnergy)[_posI1][_rotT1][_posI2][_rotT2];
	//cout << (*pairEnergy)[_posI1][_rotR1][_posI2][_rotR2] << " " << (*pairEnergy)[_posI1][_rotT1][_posI2][_rotT2] << endl;
	double Et = -Er;

//DEE SGP: pair pos/rot-pos/rot 10/2-6/1 flagged for elimination by 10/6-6/30
	
	// for all other positions
	for (unsigned int posJ=0; posJ<pairEnergy->size(); posJ++) {
		double min = 0.0;
		double max = 0.0;
		if (posJ==_posI1 || posJ==_posI2) {
			continue;
		}

		minDiffIrItJuDouble(_posI1, _rotR1, _rotT1, _posI2, _rotR2, _rotT2, posJ, min, max);
		//cout << Er << " " << _posI1 << " " <<  _rotR1 << " " <<  _rotT1 << " " << _posI2 << " " << _rotR2 << " " <<  _rotT2 << " " << posJ << " " << min << " " << max << endl;

		Er += min;
		Et -= max;

	}

	if (Er > enerOffset) {
		// flag rotamer pair R
		flaggedPair[_posI1][_rotR1][_posI2][_rotR2] = true;
		if (verboseLevel3_flag) {
			cout << "DEE SGP: pair pos/rot-pos/rot " << _posI1 << "/" << _rotR1 << "-" << _posI2 << "/" << _rotR2 << " flagged for elimination by " << _posI1 << "/" << _rotT1 << "-" << _posI2 << "/" << _rotT2 << endl;
		}
		return true; 
	} else if (Et > enerOffset) {
		// flag rotamer pair T
		flaggedPair[_posI1][_rotT1][_posI2][_rotT2] = true;
		if (verboseLevel3_flag) {
			cout << "DEE SGP: pair pos/rot-pos/rot " << _posI1 << "/" << _rotT1 << "-" << _posI2 << "/" << _rotT2 << " flagged for elimination by " << _posI1 << "/" << _rotR1 << "-" << _posI2 << "/" << _rotR2 << endl;
		}
		return true;
	}
	return false; // nothing flagged

	
}

void DeadEndElimination::minDiffIrItJuDouble(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2, unsigned int _posJ, double & _min, double & _max) {
	
	/*****************************************
	 *  Returns
	 * 
	 *   min[E(IrJu)-E(ItJu)]
	 *    u
	 *
	 *  for the super-rotamer I made by pos I1 and I2
	 *****************************************/
	if (_posI1 < _posI2) {
		cerr << "ERROR 3852: out of order positions at void DeadEndElimination::minDiffIrItJuDouble(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2, unsigned int _posJ, double & _min, double & _max)" << endl;
		exit(3852);
	}
	
	// reset min and max in case nothing is found
	_min = 0.0;
	_max = 0.0;
	
	bool found = false;
	// the table of pair energies (with N positions) is the half of the matrix with I = 0 -> N and J = 0 -> (I-1)
	if (_posI2 > _posJ) {
		for (unsigned int jr=0; jr<(*pairEnergy)[_posI1][_rotR1][_posJ].size(); jr++) {
			double diff = (*pairEnergy)[_posI1][_rotR1][_posJ][jr] + (*pairEnergy)[_posI2][_rotR2][_posJ][jr] - (*pairEnergy)[_posI1][_rotT1][_posJ][jr] - (*pairEnergy)[_posI2][_rotT2][_posJ][jr];

			if (!found) {
				found = true;
				_min = diff;
				_max = diff;
			} else {
				if (diff < _min) {
					_min = diff;
				}
				if (diff > _max) {
					_max = diff;
				}
			}
		}
	} else if (_posI1 > _posJ) {
		for (unsigned int jr=0; jr<(*pairEnergy)[_posI1][_rotR1][_posJ].size(); jr++) {
			double diff = (*pairEnergy)[_posI1][_rotR1][_posJ][jr] + (*pairEnergy)[_posJ][jr][_posI2][_rotR2] - (*pairEnergy)[_posI1][_rotT1][_posJ][jr] - (*pairEnergy)[_posJ][jr][_posI2][_rotT2];
			if (!found) {
				found = true;
				_min = diff;
				_max = diff;
			} else {
				if (diff < _min) {
					_min = diff;
				}
				if (diff > _max) {
					_max = diff;
				}
			}
		}
	} else {
		for (unsigned int jr=0; jr<(*pairEnergy)[_posJ].size(); jr++) {
			double diff = (*pairEnergy)[_posJ][jr][_posI1][_rotR1] + (*pairEnergy)[_posJ][jr][_posI2][_rotR2] - (*pairEnergy)[_posJ][jr][_posI1][_rotT1] - (*pairEnergy)[_posJ][jr][_posI2][_rotT2];
			if (!found) {
				found = true;
				_min = diff;
				_max = diff;
			} else {
				if (diff < _min) {
					_min = diff;
				}
				if (diff > _max) {
					_max = diff;
				}
			}
		}
	}
	/*
	if (!found) {
		cerr << "ERROR 3856: everything eliminated at position " << _posJ << " at void DeadEndElimination::minDiffIrItJuDouble(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2, unsigned int _posJ, double & _min, double & _max)" << endl;
		exit(3856);
	}
	*/
}


vector<vector<unsigned int> > DeadEndElimination::getAliveStates() const {

	vector<vector<unsigned int> > aliveStates;
	
	for (unsigned int i=0; i<alive.size(); i++) {
		aliveStates.push_back(vector<unsigned int>());
		for (unsigned int j=0; j<alive[i].size(); j++) {
			if (alive[i][j]) {
				aliveStates[i].push_back(j);
			}
		}	
	}	
	return aliveStates;
}
	
bool DeadEndElimination::writeMaskAsBinaryFile(const string _filename) const {

	ofstream dataFile_fs(_filename.c_str(), ios::out | ios::binary);
	if (dataFile_fs.fail()) {
		cerr << "WARNING 4244: cannot write mask table in data file " << _filename << " in bool DeadEndElimination::writeMaskAsBinaryFile(const string _filename) const" << endl;
		return false;
	}

	int expectedSize = 0;
	
	/*************************************************
	 *  Write the table
	 **************************************************/ 
	// dimension and sizes
	int dimension = 2;
	dataFile_fs.write((char*)&dimension, sizeof(int));
	expectedSize += sizeof(int);
	int size = alive.size();
	dataFile_fs.write((char*)&size, sizeof(int));
	expectedSize += sizeof(int);
	for (unsigned int i=0; i<alive.size(); i++) {
		size = alive[i].size();
		dataFile_fs.write((char*)&size, sizeof(int));
		expectedSize += sizeof(int);
	}

	// data
	for (unsigned int i=0; i<alive.size(); i++) {
		for (unsigned int j=0; j<alive[i].size(); j++) {
			// waste of bytes converting the bool to int but
			// I don't know how to write bool to binary file
			// (each bool doesn't take a number of bytes but a bit
			int data = 0;
			if (alive[i][j]) {
				data = 1;
			}
			dataFile_fs.write((char*)&data, sizeof(int));
		}
		expectedSize += sizeof(int) * alive[i].size();
	}


	/*************************************************
	 *  Done writing file: error checking and close
	 **************************************************/ 
	if (dataFile_fs.fail()) {
		cerr << "WARNING 4248: error while writing the mask table in data file " << _filename << " in bool DeadEndElimination::writeMaskAsBinaryFile(const string _filename) const" << endl;
		return false;
	}

	dataFile_fs.close();




	/*************************************************
	 *  Error checking: stat the file and check the file size
	 **************************************************/ 
	struct stat results;
	if (stat(_filename.c_str(), &results) == 0) {
		if (expectedSize != results.st_size) {
			cerr << "WARNING 5499: unexpected file size after writing the mask table in data file " << _filename << " in bool DeadEndElimination::writeMaskAsBinaryFile(const string _filename) const" << endl;
			return false;
		}
	} else {
		cerr << "WARNING 5504: canoot stat the mask table data file " << _filename << " in bool DeadEndElimination::writeMaskAsBinaryFile(const string _filename) const" << endl;
		return false;
	}

	return true;


}



void DeadEndElimination::readEnergyTable(string _filename){

	/*
	  ENERGY TABLE FILE FORMAT:
	  
	  # COMMENT LINES
	  # SELF TERMS
	  #  POSITION_INDEX ROTAMER_INDEX SELF_ENERGY
	  0 0 -0.5
	  0 1 -0.3
	  1 0  0.1
	  ..

	  # PAIR TERMS
	  #  POSITION_INDEX ROTAMER_INDEX POSITION_INDEX ROTAMER_INDEX PAIR_ENERGY
	  0 0 1 0 -0.2
	  0 0 1 1  0.5
	  ...
	  
	*/

	// This object is now responsible for the energy table memory.
	responsibleForEnergyTableMemory = true;

	selfEnergy = new vector<vector<double> >();
	pairEnergy = new vector<vector<vector<vector<double> > > >();

	ifstream fin;
	string line;
	try {
		fin.open(_filename.c_str());
		if (fin.fail()){
			cerr << "ERROR 8904 in DeadEndElimination::readEnergyTable opening file "<<_filename<<endl;
			exit(8904);
		}
		
		bool firstPair = true;
		while (!fin.fail()){
			getline(fin,line);
			if (line[0] == '#') continue;
			if (line.length() <= 1) continue;



			vector<string> toks = MslTools::tokenize(line);

			
			// self Energy Line has 3 numbers on it
			if (toks.size() == 3){
				int posIndex = MslTools::toInt(toks[0]);
				if (selfEnergy->size() < posIndex+1){
					vector<double> tmp;
					selfEnergy->push_back(tmp);
				}

				(*selfEnergy)[posIndex].push_back(MslTools::toDouble(toks[2]));
				
			}



			// pair Energy Line has 5 numbers on it
			/*
			  TODO:
			     Make sure pairEnergy is not a full sized table.
			      it should only be a triangle - lower triangular.
			      resize(i+1) for pairEnergy[i]?
			 */
			if (toks.size() == 5){

				if (firstPair){
					pairEnergy->resize(selfEnergy->size());
					// For each position resize to
					// num rotamers
					for (uint i = 0; i < selfEnergy->size();i++){
						(*pairEnergy)[i].resize((*selfEnergy)[i].size());

						// For each rotamer,
						// resize to num positions
						for (uint j = 0; j < (*selfEnergy)[i].size();j++){
							(*pairEnergy)[i][j].resize((*selfEnergy).size());

							// For each position resize to num rotamers
							for (uint k = 0; k < (*selfEnergy).size();k++){
								(*pairEnergy)[i][j][k].resize((*selfEnergy)[k].size());

							}
						}

					}
				}
				firstPair = false;
				
				// Add to pairEnergy object. Indices pairEnergy[POS1][ROT1][POS2][ROT2] = ENERGY;
				(*pairEnergy)[MslTools::toInt(toks[0])][MslTools::toInt(toks[1])][MslTools::toInt(toks[2])][MslTools::toInt(toks[3])] = MslTools::toDouble(toks[4]);

				// Add symmetric entries into the pairEnergy table ****NO NEED****
				// pairEnergy[MslTools::toInt(toks[2])][MslTools::toInt(toks[3])][MslTools::toInt(toks[0])][MslTools::toInt(toks[1])] = MslTools::toDouble(toks[4]);
				
			}
		}
	} catch (...){
		cerr << "ERROR 7809 in DeadEndElimination::readEnergyTable in file: "<<_filename<<endl;
		exit(7809);
	}


	totalNumPositions = (*selfEnergy).size();
	alive.clear();
	alive.resize(totalNumPositions);
	for (uint i = 0; i < alive.size();i++){
		alive[i].resize((*selfEnergy)[i].size());
		totalNumRotamers += alive[i].size();
		for (uint j = 0; j < alive[i].size();j++){
			alive[i][j] = true;
		}
	}

}


void DeadEndElimination::printMe(bool _selfOnly){
	fprintf(stdout,"Self terms:\n");
	for (uint i = 0; i < (*selfEnergy).size();i++){
		for (uint j = 0; j < (*selfEnergy)[i].size();j++){

			fprintf(stdout, "    %4d %4d %8.3f", i, j, (*selfEnergy)[i][j]);
			if (alive[i][j]) {
				fprintf(stdout, " **** ");
			}
			fprintf(stdout,"\n");


		}

	}

	if (_selfOnly){
		return;
	}

	fprintf(stdout,"Pair terms:\n");
	for (uint i = 0; i < (*pairEnergy).size();i++){
		for (uint j = 0; j < (*pairEnergy)[i].size();j++){
			//for (uint k = i+1 ; k < (*pairEnergy)[i][j].size();k++){	
			for (uint k = 0 ; k < (*pairEnergy)[i][j].size();k++){	
				if (i == k){
					continue;
				}
				for (uint l = 0 ; l < (*pairEnergy)[i][j][k].size();l++){	
					fprintf(stdout, "    %4d %4d %4d %4d %8.3f", i, j, k, l, (*pairEnergy)[i][j][k][l]);

					if (alive[i][j] && alive[k][l]) {
						fprintf(stdout, " **** ");
					}
					fprintf(stdout,"\n");
				}
			}
		}

	}
}
