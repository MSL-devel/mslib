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

#ifndef MONTECARLOOPTIMIZATION_H
#define MONTECARLOOPTIMIZATION_H

#include <vector>
#include <string>
#include <map>
#include <queue>
#include <algorithm>
#include <math.h>

#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "SelfPairManager.h"
#include "MslTools.h"


namespace MSL { 
class SelfPairManager;
class MonteCarloOptimization {


	public:

	        // Constructors
		MonteCarloOptimization();
		MonteCarloOptimization(SelfPairManager *_pSpm); // must be set if in onTheFlyMode
		~MonteCarloOptimization();


		// Enum. Types
		enum INITTYPE    { RANDOM=0, LOWESTSELF=1, QUICKSCAN=2, USERDEF=3};


		// Read or Add Energy information
		void readEnergyTable(std::string _filename);
		// The pairTable has to be lower triangular.  
		void addEnergyTable(std::vector<std::vector<double> > &_selfEnergy, std::vector<std::vector<std::vector<std::vector<double> > > > &_pairEnergy); 
		void setSelfPairManager(SelfPairManager* _pSpm);//must be set if in onTheFlyMode


		// Run the MonteCarlo and get the best State
		std::vector<unsigned int> runMC(double _startingTemperature, double _endingTemperature, int _scheduleCycles, int _scheduleShape, int _maxRejectionsNumber, int _convergedSteps, double _convergedE);

		// Print
		void printMe(bool _selfOnly="true");
		void printSampledConfigurations();
		
		// Getter/Setters...

		// sets initType and initState appropriately
		void setInitializationState(INITTYPE _type, std::string _state="");
		void setInitializationState(std::vector<unsigned int>& _state); // set the starting state to _state and type to USERDEF		
		int getInitialization();
		
		std::vector<unsigned int>& getInitState(); /*! return the starting state */

		void setCurrentState(std::vector<unsigned int> _state);
		std::vector<unsigned int>& getCurrentState(); /*! return the current state */
		
		std::vector<unsigned int> getRandomState(); /*! return a state selected uniformly at random  and update current state*/
		std::vector<unsigned int> moveRandomState(); /*! change randomly just one position of the current state */
		std::vector<unsigned int> moveRandomState(unsigned int _numOfSteps); /*! change randomly the specified number of positions of the current state */

		int getNumPositions();
		int getNumRotamers(int _index); 

		void setNumberOfStoredConfigurations(int _numConfs);
		int getNumberOfStoredConfigurations();
		
		void seed(unsigned int _seed);
		unsigned int getSeed() const;

		double getEnergy(int _pos, int _rot);
		double getStateEnergy(std::vector<unsigned int> _states);
		double getStateEnergy();

		std::vector<std::vector<bool> > getMask(); // everything except the best state will be masked out

		void setInputRotamerMasks(std::vector<std::vector<bool> > &_inputMasks); // true if rotamer is alive
		//void linkPositions(int _pos1, int _pos2); // not implemented, what's for?


		std::priority_queue< std::pair<double,std::string>, std::vector< std::pair<double,std::string> >, std::less<std::pair<double,std::string> > > & getSampledConformations();
		
		void setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG);
		RandomNumberGenerator * getRandomNumberGenerator() const;


	private:	    
		
		// Initialize the system - update current state and initState appropriately
		void setup();
		void initialize();
		void deletePointers();
		void deleteEnergyTables();
		void selectRotamer(int _pos,int _rot); // update currentState
		int selectRandomStateAtPosition(int _position) const;
		std::string getRotString(int _pos, int _rot);
		std::string getRotString();


		// Member Variables
		std::vector<std::vector<double> > *selfEnergy;
		std::vector<std::vector<std::vector<std::vector<double > > > > *pairEnergy;
		std::vector<std::vector<bool> > inputMasks; 
		std::map<std::string,double> configurationMap;

		/*
		  Functions tha linkedPositions need to be added:
		   selectRotamer
		   getEnergy
		   initialize
		 */

		//std::map<int,int> linkedPositions; // unused, remove?
		int numStoredConfigurations;

		int initType;
		std::vector<unsigned int> initState;
		std::vector<unsigned int> currentState;
		std::vector<unsigned int> bestState;


		// Energy table parameters
		int totalNumPositions;
		bool responsibleForEnergyTableMemory;

		// Utility variables
		RandomNumberGenerator * pRng;
		SelfPairManager * pSpm;
		bool deleteRng;

		std::priority_queue< std::pair<double,std::string>, std::vector< std::pair<double,std::string> >, std::less<std::pair<double,std::string> > > sampledConfigurations;

};


inline void MonteCarloOptimization::setCurrentState(std::vector<unsigned int> _state) { currentState = _state;}
inline std::vector<unsigned int>& MonteCarloOptimization::getCurrentState() { return currentState;}
inline int MonteCarloOptimization::getInitialization() { return initType; }
inline std::vector<unsigned int>& MonteCarloOptimization::getInitState() { return initState; }


inline void MonteCarloOptimization::setNumberOfStoredConfigurations(int _numConfs){ numStoredConfigurations = _numConfs; }
inline int MonteCarloOptimization::getNumberOfStoredConfigurations() { return numStoredConfigurations; }

inline int MonteCarloOptimization::getNumPositions() { if (selfEnergy == NULL) { return 0; } return (*selfEnergy).size();}
inline int MonteCarloOptimization::getNumRotamers(int _index) { 
	if (selfEnergy == NULL || _index >= (*selfEnergy).size()) { 
		return 0; 
	} 

	return (*selfEnergy)[_index].size();
}

inline void MonteCarloOptimization::setInputRotamerMasks(std::vector<std::vector<bool> > &_inputMasks) { inputMasks = _inputMasks; }


//inline void MonteCarloOptimization::linkedPositions(int _pos1, int _pos2){
//
//	linkedPositions[_pos1].push_back(_pos2);
//	linkedPositions[_pos2].push_back(_pos1);
//}


inline std::priority_queue< std::pair<double,std::string>, std::vector< std::pair<double,std::string> >, std::less<std::pair<double,std::string> > > & MonteCarloOptimization::getSampledConformations() { return sampledConfigurations; }

inline void MonteCarloOptimization::setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}
inline RandomNumberGenerator * MonteCarloOptimization::getRandomNumberGenerator() const {
	return pRng;
}
inline void MonteCarloOptimization::seed(unsigned int _seed){
	if (_seed == 0) {
		pRng->setTimeBasedSeed();
	}
	else {
		pRng->setSeed(_seed);
	}
}


inline unsigned int MonteCarloOptimization::getSeed() const {
	return pRng->getSeed();
}
}

#endif
