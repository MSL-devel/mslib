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


namespace MSL { 
class LinearProgrammingOptimization {


	public:
		LinearProgrammingOptimization();
		~LinearProgrammingOptimization();

		void readEnergyTable(std::string _filename);
		void addEnergyTable(std::vector<std::vector<double> > &_selfEnergy, std::vector<std::vector<std::vector<std::vector<double> > > > &_pairEnergy);

		void analyzeEnergyTable();
	        void createLP();
		void solveLP();
		void printMe(bool _selfOnly=true);

		void setVerbose(bool _flag){ verbose = _flag;}
		bool getVerbose() { return verbose; }

		std::vector<std::vector<bool> > getMask();
		double getTotalEnergy();
		
		void setInputRotamerMasks(std::vector<std::vector<bool> > &_inputMasks);
		std::string getRotString();

		int getNumPositions();
		int getNumRotamers(int _index);

		std::vector<int>& getRotamerSelection();
		
	private:	    
		
		std::vector<std::vector<double> > *selfEnergy;
		std::vector<std::vector<std::vector<std::vector<double > > > > *pairEnergy;
		std::vector<std::vector<bool> > pairType;
		std::vector<std::vector<bool> > masks;
		std::vector<int> rotamerSelection;
		std::vector<std::vector<bool> > inputMasks;

		// Decision Variables
		std::vector<std::vector<int> > rotamerState;
		std::vector<std::vector<int> > rotamerPair;

		int totalNumRotamers;
		int totalNumPositions;
		int numConstraints;
		int numDecisionVariables;
		

		bool responsibleForEnergyTableMemory;
		glp_prob *lp;


		bool verbose;
};

inline void LinearProgrammingOptimization::setInputRotamerMasks(std::vector<std::vector<bool> > &_inputMasks) { inputMasks = _inputMasks; }
inline int LinearProgrammingOptimization::getNumPositions() { if (selfEnergy == NULL) { return 0; } return (*selfEnergy).size();}
inline int LinearProgrammingOptimization::getNumRotamers(int _index) { 
	if (selfEnergy == NULL || _index >= (*selfEnergy).size()) { 
		return 0; 
	} 

	return (*selfEnergy)[_index].size();
}

inline std::vector<int> & LinearProgrammingOptimization::getRotamerSelection() { return rotamerSelection; }

}

#endif
