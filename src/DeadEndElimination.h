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

#ifndef DEADENDELIM_H
#define DEADENDELIM_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
//#include <math.h>
#include "MslTools.h"
using namespace std;

/*! \brief Dead End Elimination class
 */

class DeadEndElimination {
	public:
		DeadEndElimination();
		DeadEndElimination(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies);
		DeadEndElimination(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		~DeadEndElimination();


		void readEnergyTable(string _filename);
		void setEnergyTables(vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies);
		void setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies);
		void setEnergyTables(vector<vector<double> > & _selfEnergies, vector<vector<vector<vector<double> > > > & _pairEnergies, vector<vector<double> > & _baselines);
		void setEnergyTables(vector<vector<double> > * _pSelfEnergies, vector<vector<vector<vector<double> > > > * _pPairEnergies, vector<vector<double> > * _pBaselines);

		void setBaselines(vector<vector<double> > & _baselines);
		void setBaselines(vector<vector<double> > * _pBaselines);

		void setEnergyOffset(double _offset);
		double getEnergyOffset() const;

		void setVerbose(bool _flag);
		void setVerbose(bool _flag, unsigned int _level);
		bool getVerbose() const;
		unsigned int getVerboseLevel() const;
		unsigned int getEliminatedCounter() const;
		unsigned int getFlaggedCounter() const;
		unsigned int getRunCycles() const;
		
		void setMask(vector<vector<bool> > theMask);
		vector<vector<bool> > getMask() const;

		bool runSimpleGoldsteinSingles();
		bool runSimpleGoldsteinPairs(); // buggy, the code will exit if it is called
		
		double getTotalCombinations() const;

		bool isAlive(vector<int> _states) const;

		vector<vector<int> > getAliveStates() const;

		bool writeMaskAsBinaryFile(const string _filename) const;


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

		vector<vector<double> > * selfEnergy;
		vector<vector<vector<vector<double> > > > * pairEnergy;
		vector<vector<double> > * pBaseLines;
		

		// the mask specifying the living (true) vs eliminated (false) rotamers
		vector<vector<bool> > alive;
		vector<vector<vector<vector<bool> > > > flaggedPair;
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
#endif
