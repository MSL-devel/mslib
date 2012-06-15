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

#ifndef DEADENDELIM_H
#define DEADENDELIM_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
//#include <math.h>
#include "MslTools.h"

/*! \brief Dead End Elimination class
 */

namespace MSL { 
class DeadEndElimination {
	public:
		DeadEndElimination();
		DeadEndElimination(std::vector<std::vector<double> > & _selfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > & _pairEnergies);
		DeadEndElimination(std::vector<std::vector<double> > & _selfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > & _pairEnergies, std::vector<std::vector<double> > & _baselines);
		~DeadEndElimination();


		void readEnergyTable(std::string _filename);
		void setEnergyTables(std::vector<std::vector<double> > * _pSelfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > * _pPairEnergies);
		void setEnergyTables(std::vector<std::vector<double> > & _selfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > & _pairEnergies);
		void setEnergyTables(std::vector<std::vector<double> > & _selfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > & _pairEnergies, std::vector<std::vector<double> > & _baselines);
		void setEnergyTables(std::vector<std::vector<double> > * _pSelfEnergies, std::vector<std::vector<std::vector<std::vector<double> > > > * _pPairEnergies, std::vector<std::vector<double> > * _pBaselines);

		void setBaselines(std::vector<std::vector<double> > & _baselines);
		void setBaselines(std::vector<std::vector<double> > * _pBaselines);

		void setEnergyOffset(double _offset);
		double getEnergyOffset() const;

		void setVerbose(bool _flag);
		void setVerbose(bool _flag, unsigned int _level);
		bool getVerbose() const;
		unsigned int getVerboseLevel() const;
		unsigned int getEliminatedCounter() const;
		unsigned int getFlaggedCounter() const;
		unsigned int getRunCycles() const;
		
		void setMask(std::vector<std::vector<bool> > theMask);
		std::vector<std::vector<bool> > getMask() const;

		bool runSimpleGoldsteinSingles();
		bool runSimpleGoldsteinPairs();  // might take too long
		unsigned int runSimpleGoldsteinPairsOnce();  // run one iteration of Pairs, may be used for large optimization problems

		double getTotalCombinations() const;

		bool isAlive(std::vector<int> _states) const;

		std::vector<std::vector<unsigned int> > getAliveStates() const;

		bool writeMaskAsBinaryFile(const std::string _filename) const;


	        void printMe(bool _selfOnly="true");

		int getTotalNumberRotamers();
	private:
		
		void setInitialVariables();
		void initializeMask();
		//void setInitialVariables();
		bool simpleGoldsteinSinglesIteration(unsigned int _posI, unsigned int _rotR, unsigned int _rotT);
		bool simpleGoldsteinPairIteration(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2);
		void minDiffIrItJuSingle(unsigned int _posI, unsigned int _rotR, unsigned int _rotT, unsigned int _posJ, double & _min, double & _max);
		void minDiffIrItJuSingleAfterPair(unsigned int _posI, unsigned int _rotR, unsigned int _rotT, unsigned int _posJ, double & _min, double & _max, bool & _eliminateR, bool & _eliminateT);
		void minDiffIrItJuDouble(unsigned int _posI1, unsigned int _rotR1, unsigned int _rotT1, unsigned int _posI2, unsigned int _rotR2, unsigned int _rotT2, unsigned int _posJ, double & _min, double & _max);

		std::vector<std::vector<double> > * selfEnergy;
		std::vector<std::vector<std::vector<std::vector<double> > > > * pairEnergy;
		std::vector<std::vector<double> > * pBaseLines;
		

		// the mask specifying the living (true) vs eliminated (false) rotamers
		std::vector<std::vector<bool> > alive;
		std::vector<std::vector<std::vector<std::vector<bool> > > > flaggedPair;
		unsigned int verboseLevel; // there are 4 levels, 0 (none), 1 (little), 2 (some), 3 (very)
		bool verboseLevel1_flag;
		bool verboseLevel2_flag;
		bool verboseLevel3_flag;
		unsigned int cycles;
		unsigned int eliminatedCounter;
		unsigned int flaggedCounter;
		double enerOffset;
		bool afterPair_flag;

		bool responsibleForEnergyTableMemory;

		int totalNumPositions;
		int totalNumRotamers;

};

inline int DeadEndElimination::getTotalNumberRotamers() { return totalNumRotamers; }
}

#endif
