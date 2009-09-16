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

#ifndef QUENCH_H
#define QUENCH_H

// MSL Includes
#include <string>
#include <sstream>
#include <vector>
#include <ostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"
#include "SystemRotamerLoader.h"
#include "PairwiseEnergyCalculator.h"
#include "AtomicPairwiseEnergy.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "RandomNumberGenerator.h"

// Namespaces
using namespace std;


class Quench {

	public:
		Quench();
		Quench(string _topfile, string _parfile, string _rotlib);
		~Quench();

		System runQuench(System & _initialSystem);
		System runQuench(System & _initialSystem, vector<int> variablePositions);
		System runQuench(System & _initialSystem, uint _numIterations);

		void runPreSetUpQuench(System & _mySystem);
		void runPreSetUpQuench(System & _mySystem, uint _numIterations);

		void setUpSystem(System & _initialSystem, System & _outputSystem);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers);
		void setUpSystem(System & _initialSystem, System & _outputSystem, vector<int> variablePositions);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, vector<int> variablePositions);

		void setUpMonomericSurroundEnergies(System & _mySystem);
		double runPreSetUpQuenchOnDimer(System & _mySystem); // Returns CHARMM energy
		double runPreSetUpQuenchOnDimer(System & _mySystem, uint _numIterations); // Returns CHARMM energy


		void setVariableNumberRotamers(int _largeSideChainsNumRot, int _smallSideChainsNumRot);

	protected:

		string topfile;
		string parfile;
		string rotlib;

		int numberLargeRotamers;
		int numberSmallRotamers;

		AtomicPairwiseEnergy ape;
		PairwiseEnergyCalculator pec;
		vector<uint> currentRotamers;
		vector<uint> currentAllRotamers;
		map<string, vector<double> > selfEnergies;
		vector < vector < vector < vector<double> > > > monomericSurroundEnergies;

};
inline void Quench::setVariableNumberRotamers(int _largeSideChainsNumRot, int _smallSideChainsNumRot) { numberLargeRotamers = _largeSideChainsNumRot; numberSmallRotamers = _smallSideChainsNumRot; }
#endif
