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

#ifndef LINEARPROGRAMMINGOPTIMIZATION_H
#define LINEARPROGRAMMINGOPTIMIZATION_H

#include <vector>
#include <string>
#include <ctime>

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
		void solveLP(); // uses the rotamer with the most weight for each position and updates rotamerSelection 
		void solveMIP();
		void printMe(bool _selfOnly=true);
		int writeCPLEXFile(std::string _filename); // write the LP in CPLEX format - returns 0 on success

		void setVerbose(bool _flag){ verbose = _flag;}
		bool getVerbose() { return verbose; }

		std::vector<std::vector<bool> > getMask(); // masks out everything except the solution rotamers
		double getTotalEnergy();
		
		void setInputRotamerMasks(std::vector<std::vector<bool> > &_inputMasks);
		std::string getRotString();

		int getNumPositions();
		int getNumRotamers(int _index);

		std::vector<unsigned int>& getRotamerSelection();
		// call after setting up the energytables - returns the final rotamerSelection
		std::vector<unsigned int>& getSolution(bool _runMIP);

		
	private:	    
		
		void deletePointers();
		void deleteEnergyTables();

		std::vector<std::vector<double> > *selfEnergy;
		std::vector<std::vector<std::vector<std::vector<double > > > > *pairEnergy;
		std::vector<std::vector<bool> > pairType;
		std::vector<unsigned int> rotamerSelection;
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

inline std::vector<unsigned int> & LinearProgrammingOptimization::getRotamerSelection() { return rotamerSelection; }

}

#endif
