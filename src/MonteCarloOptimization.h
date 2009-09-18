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

#ifndef MONTECARLOOPTIMIZATION_H
#define MONTECARLOOPTIMIZATION_H

#include <vector>
#include <string>
#include <map>
#include <queue>

#include "RandomNumberGenerator.h"
#include "MslTools.h"

using namespace std;
using namespace MslTools;

class MonteCarloOptimization {


	public:

	        // Constructors
		MonteCarloOptimization();
		~MonteCarloOptimization();


		// Enum. Types
		enum ANNEALTYPES { NO_ANNEAL=0, LIN_TEMP_ANNEAL=1, EXP_TEMP_ANNEAL=2, SAWTOOTH_TEMP_ANNEAL=3,EXPCYCLE_TEMP_ANNEAL=4 } ;
		enum INITTYPE    { RANDOM=0, LOWESTSELF=1, QUICKSCAN=2, USERDEF=3};


		// Read or Add Energy information
		void readEnergyTable(string _filename);
		void addEnergyTable(vector<vector<double> > &_selfEnergy, vector<vector<vector<vector<double> > > > &_pairEnergy);


		// Initialize the system (choose rotamers for each position)
		void initialize();

		// Run the MonteCarlo
		void runMC();

		// Print
		void printMe(bool _selfOnly="true");
		void printSampledConfigurations();
		
		// Getter/Setters...

		// Number of Cycles of MC to do
		void setNumberOfCycles(int _cycles);
		int getNumberOfCycles();

		void setInitializationState(int _state, string _userDef="");
		int getInitializationState();
		
		void setAnnealSchedule(int _annealType, double _startTemp, double _endTemp,int _numCycles=1);

		void setVerbose(bool _flag){ verbose = _flag;}
		bool getVerbose() { return verbose; }

		void setNumberOfStoredConfigurations(int _numConfs);
		int getNumberOfStoredConfigurations();
		
		void setRandomSeed(int _seed);
		int getRandomSeed();

		double getEnergy(int _pos, int _rot);
		double getTotalEnergy();
		
		int getCurrentStep(); // can't set this one.

		vector<vector<bool> > getMask();

		int getNumPositions();
		int getNumRotamers(int _index);

		void setInputRotamerMasks(vector<vector<bool> > &_inputMasks);
		void linkPositions(int _pos1, int _pos2);


		priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > & getSampledConformations();
		
		
	private:	    
		
		void selectRotamer(int _pos,int _rot);
		vector<int> getRandomRotamer();

		void setCurrentTemp(double _temp);
		double getCurrentTemp();

		string getRotString(int _pos, int _rot);
		string getRotString();
		void annealTemperature(double initialTemp, double finalTemp, int step, int totalsteps);


		// Member Variables
		vector<vector<double> > *selfEnergy;
		vector<vector<vector<vector<double > > > > *pairEnergy;
		vector<vector<bool> > masks;
		vector<vector<bool> > inputMasks;
		vector<int> rotamerSelection;
		map<string,double> configurationMap;

		/*
		  Functions tha linkedPositions need to be added:
		   selectRotamer
		   getEnergy
		   initialize
		 */

		map<int,int> linkedPositions;
		int numStoredConfigurations;

		// MC Parameters
		int cycles;
		int annealType;
		double startTemp;
		double endTemp;
		double temp;
		int numAnnealCycles;
		int initType;
		string initConf;
		int currentStep;


		// Energy table parameters
		int totalNumRotamers;
		int totalNumPositions;
		bool responsibleForEnergyTableMemory;

		// Utility variables
		bool verbose;
		RandomNumberGenerator rng;
		int randomSeed;

		priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > sampledConfigurations;

};

inline void MonteCarloOptimization::setNumberOfCycles(int _cycles) { cycles = _cycles; }
inline int MonteCarloOptimization::getNumberOfCycles() { return cycles; }


inline void MonteCarloOptimization::setCurrentTemp(double _temp) { temp  = _temp; }
inline double MonteCarloOptimization::getCurrentTemp() { return temp; }

inline void MonteCarloOptimization::setAnnealSchedule(int _annealType, double _startTemp, double _endTemp, int _cycles) { 
	annealType = _annealType;
	startTemp  = _startTemp;
	endTemp    = _endTemp;
	numAnnealCycles  = _cycles; 
	
}


inline void MonteCarloOptimization::setInitializationState(int _state, string _userDef){
	initType = _state;
	initConf = _userDef;
}
inline int MonteCarloOptimization::getInitializationState() { return initType; }


inline void MonteCarloOptimization::setNumberOfStoredConfigurations(int _numConfs){ numStoredConfigurations = _numConfs; }
inline int MonteCarloOptimization::getNumberOfStoredConfigurations() { return numStoredConfigurations; }


inline void MonteCarloOptimization::setRandomSeed(int _seed){ randomSeed = _seed;}
inline int MonteCarloOptimization::getRandomSeed(){ return randomSeed; }

inline int MonteCarloOptimization::getCurrentStep(){ return currentStep; }

inline int MonteCarloOptimization::getNumPositions() { if (selfEnergy == NULL) { return 0; } return (*selfEnergy).size();}
inline int MonteCarloOptimization::getNumRotamers(int _index) { 
	if (selfEnergy == NULL || _index >= (*selfEnergy).size()) { 
		return 0; 
	} 

	return (*selfEnergy)[_index].size();
}

inline void MonteCarloOptimization::setInputRotamerMasks(vector<vector<bool> > &_inputMasks) { inputMasks = _inputMasks; }


//inline void MonteCarloOptimization::linkedPositions(int _pos1, int _pos2){
//
//	linkedPositions[_pos1].push_back(_pos2);
//	linkedPositions[_pos2].push_back(_pos1);
//}


inline priority_queue< pair<double,string>, vector< pair<double,string> >, less<pair<double,string> > > & MonteCarloOptimization::getSampledConformations() { return sampledConfigurations; }

#endif
