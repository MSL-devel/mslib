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
#include "CharmmEnergyCalculator.h"
#include "OnTheFlyManager.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "RandomNumberGenerator.h"
#include "TwoBodyDistanceDependentPotentialTable.h"

// Namespaces


namespace MSL { 
class Quench {

	public:
		Quench();
		Quench(std::string _topfile, std::string _parfile, std::string _rotlib);
		~Quench();

		System runQuench(System & _initialSystem);
		System runQuench(System & _initialSystem, uint _numIterations);
		System runQuench(System & _initialSystem, std::vector<int> & variablePositions);

		void runPreSetUpQuench(System & _mySystem);
		void runPreSetUpQuench(System & _mySystem, uint _numIterations);

		void setUpSystem(System & _initialSystem, System & _outputSystem);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers);
		void setUpSystem(System & _initialSystem, System & _outputSystem, std::vector<int> & variablePositions);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, std::vector<int> & variablePositions);

		void setUpMonomericSurroundEnergies(System & _mySystem);
		double runPreSetUpQuenchOnDimer(System & _mySystem); // Returns CHARMM energy
		double runPreSetUpQuenchOnDimer(System & _mySystem, uint _numIterations); // Returns CHARMM energy

		void setRotamerLevel(std::string _rotLevel);
		void setVariableNumberRotamers(int _largeSideChainsNumRot, int _smallSideChainsNumRot);

		// For knowledge based potentials
		System runQuench(System & _initialSystem, TwoBodyDistanceDependentPotentialTable & tbd);
		System runQuench(System & _initialSystem, uint _numIterations, TwoBodyDistanceDependentPotentialTable & tbd);
		System runQuench(System & _initialSystem, std::vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd);
		void runPreSetUpQuench(System & _mySystem, TwoBodyDistanceDependentPotentialTable & tbd);
		void runPreSetUpQuench(System & _mySystem, uint _numIterations, TwoBodyDistanceDependentPotentialTable & tbd);
		void setUpSystem(System & _initialSystem, System & _outputSystem, TwoBodyDistanceDependentPotentialTable & tbd);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, TwoBodyDistanceDependentPotentialTable & tbd);
		void setUpSystem(System & _initialSystem, System & _outputSystem, std::vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd);
		void setUpSystem(System & _initialSystem, System & _outputSystem, uint _numRotamers, std::vector<int> & variablePositions, TwoBodyDistanceDependentPotentialTable & tbd);

	protected:

		std::string topfile;
		std::string parfile;
		std::string rotlib;

		std::string rotLevel;
		int numberLargeRotamers;
		int numberSmallRotamers;

		CharmmEnergyCalculator calculator;
		std::vector<uint> currentRotamers;
		std::vector<uint> currentAllRotamers;
		std::map<std::string, std::vector<double> > selfEnergies;
		std::vector < std::vector < std::vector < std::vector<double> > > > monomericSurroundEnergies;

};
inline void Quench::setVariableNumberRotamers(int _largeSideChainsNumRot, int _smallSideChainsNumRot) { numberLargeRotamers = _largeSideChainsNumRot; numberSmallRotamers = _smallSideChainsNumRot;}
inline void Quench::setRotamerLevel(std::string _rotLevel) { rotLevel = _rotLevel; }

}

#endif
