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

#ifndef LINEARPROGRAMMINGOPTIMIZATION_H
#define LINEARPROGRAMMINGOPTIMIZATION_H

#include <vector>
#include <string>

#include "glpk.h"

using namespace std;

class LinearProgrammingOptimization {


	public:
		LinearProgrammingOptimization();
		~LinearProgrammingOptimization();

		void readEnergyTable(string _filename);
		void addEnergyTable(vector<vector<double> > &_selfEnergy, vector<vector<vector<vector<double> > > > &_pairEnergy);

		void analyzeEnergyTable();
	        void createLP();
		void solveLP();
		void printMe(bool _selfOnly=true);

		void setVerbose(bool _flag){ verbose = _flag;}
		bool getVerbose() { return verbose; }

		vector<vector<bool> > getMask();
		double getTotalEnergy();
		
		void setInputRotamerMasks(vector<vector<bool> > &_inputMasks);
		string getRotString();

		int getNumPositions();
		int getNumRotamers(int _index);

		vector<int>& getRotamerSelection();
		
	private:	    
		
		vector<vector<double> > *selfEnergy;
		vector<vector<vector<vector<double > > > > *pairEnergy;
		vector<vector<bool> > pairType;
		vector<vector<bool> > masks;
		vector<int> rotamerSelection;
		vector<vector<bool> > inputMasks;

		// Decision Variables
		vector<vector<int> > rotamerState;
		vector<vector<int> > rotamerPair;

		int totalNumRotamers;
		int totalNumPositions;
		int numConstraints;
		int numDecisionVariables;
		

		bool responsibleForEnergyTableMemory;
		glp_prob *lp;


		bool verbose;
};

inline void LinearProgrammingOptimization::setInputRotamerMasks(vector<vector<bool> > &_inputMasks) { inputMasks = _inputMasks; }
inline int LinearProgrammingOptimization::getNumPositions() { if (selfEnergy == NULL) { return 0; } return (*selfEnergy).size();}
inline int LinearProgrammingOptimization::getNumRotamers(int _index) { 
	if (selfEnergy == NULL || _index >= (*selfEnergy).size()) { 
		return 0; 
	} 

	return (*selfEnergy)[_index].size();
}

inline vector<int> & LinearProgrammingOptimization::getRotamerSelection() { return rotamerSelection; }

#endif
