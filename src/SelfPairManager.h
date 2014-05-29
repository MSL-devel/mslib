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

#ifndef SELFPARIMANAGER_H
#define SELFPARIMANAGER_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "DeadEndElimination.h"
#include "Enumerator.h"
#include "MonteCarloManager.h"
#include "MonteCarloOptimization.h"
#ifdef __GLPK__
	#include "LinearProgrammingOptimization.h"
#endif
#include "SelfConsistentMeanField.h"
#include "System.h"
#include "EnergySet.h"
#include "MslTools.h"

/***************************************************************
   TODO LINK POSITIONS
 ***************************************************************/

namespace MSL { 
class SelfPairManager {
	public:
		SelfPairManager();
		SelfPairManager(System * _pSystem);
		virtual ~SelfPairManager();

		void calculateEnergies(); // calls calculateFixedEnergies(),calculateSelfEnergies(), calculatePairEnergies();
		void recalculateNonSavedEnergies(std::vector<std::vector<std::vector<std::vector<bool> > > > savedPairEnergies);  // calls calculateFixedEnergies(),calculateSelfEnergies(), recalculateNonSavedPairEnergies();

		void setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG);
		RandomNumberGenerator * getRandomNumberGenerator() const;

		void seed(unsigned int _seed); // use 0 for time based seed
		unsigned int getSeed() const; // get back the seed (useful to get the time based one)
		//  GET THE ENERGIES
		double getStateEnergy(std::vector<unsigned int> _overallRotamerStates, std::string _term="");
		double getStateEnergy(std::vector<unsigned int> _residueStates, std::vector<unsigned int> _rotamerStates);
		double getStateEnergy(std::vector<std::string> _residueNames, std::vector<unsigned int> _rotamerStates);
		// may be called in onTheFly mode
		double computeSelfE(unsigned pos, unsigned rot);
		double computePairE(unsigned pos1, unsigned rot1, unsigned pos2, unsigned rot2, string _term="");

		double computeSelfE(string _posId, string _resName, string _term="");
		double computeBestPairE(string _posId, string _resName, string _posId2, string _resName2, string _term);
		unsigned int getStateInteractionCount(std::vector<unsigned int> _overallRotamerStates, std::string _term="");

		void setSystem(System * _pSystem);
		System * getSystem() const;

		unsigned int getNumberOfVariablePositions() const;
		std::vector<unsigned int> getNumberOfRotamers() const;
		double getAliveCombinations() const;

		std::string getSummary(std::vector<unsigned int> _overallRotamerStates, unsigned int _precision=6);
		std::vector<std::string> getStateDescriptors(std::vector<unsigned int> _overallRotamerStates) const;
		std::vector<std::vector<unsigned int> > getStatePositionIdentityRotamerIndeces(std::vector<unsigned int> _overallRotamerStates) const;
		
		void saveEnergiesByTerm(bool _save); // false by default
		void saveInteractionCounts(bool _save); // false by default

		double getFixedEnergy() const;

		// gets the reduced self energy table if cutoff is applied
		std::vector<std::vector<double> > & getSelfEnergy(); 		

		// gets the reduced pair energy table if cutoff is applied
		std::vector<std::vector<std::vector<std::vector<double> > > > & getPairEnergy(); 

		// Side Chain Optimization Functions
		void setRunDEE(bool _singles, bool _pairs = false); 
		void setRunSCMFBiasedMC(bool _toogle);
		void setRunUnbiasedMC(bool _toogle);
		void setRunSCMF(bool _toogle);
		void setRunEnum(bool _toogle);

		void setMCOptions(double _startT, double _endT, int _nCycles, int _shape, int _maxReject, int _deltaSteps, double _minDeltaE);

		void setOnTheFly(bool _onTheFly);
		
		void setEnumerationLimit(int _enumLimit);

		std::vector<std::vector<bool> > getSelfEnergyMask() const;

		void setVerbose(bool _toggle);

		void runOptimizer();
		//SGFC runGreedyOptimizer can accept a mask to exclude particular rotamers
		void runGreedyOptimizer(int _cycles, std::vector< std::vector<bool> > _mask);
		void runGreedyOptimizer(int _cycles) ;
		std::vector<unsigned int> runLP(bool _runMIP = false); // Run the LP/MIP formulation

		std::vector<double> getMinBound();
		std::vector<vector<unsigned int> > getMinStates();

		std::vector<std::vector<unsigned int> > getDEEAliveRotamers();  // after applying cutoff
		std::vector<std::vector<bool> > getDEEAliveMask(); // after applying cutoff

		std::vector<unsigned int> getSCMFstate();
		std::vector<unsigned int> getBestSCMFBiasedMCState();
		std::vector<unsigned int> getBestUnbiasedMCState();

		void updateWeights();

	private:
		void setup();
		void copy(const SelfPairManager & _sysBuild);
		void deletePointers();
		void findVariablePositions();
		void subdivideInteractions();

		//double computePairEbyTerm(unsigned pos1, unsigned rot1, unsigned pos2, unsigned rot2, std::string _term);
		int getNumPositions();
		int getNumRotamers(int _index);
		double getInteractionEnergy(int _pos, int _rot, vector<unsigned int>& _currentState);

		void calculateFixedEnergies();
		void calculateSelfEnergies();
		void calculatePairEnergies();
		void recalculateNonSavedPairEnergies(std::vector<std::vector<std::vector<std::vector<bool> > > > savedPairEnergies);

		double runDeadEndElimination(); // returns the finalCombinations
		void runEnumeration();
		void runSelfConsistentMeanField();
		void runUnbiasedMonteCarlo();

		bool deleteRng;

		std::map<std::string, double> weights; // weights for the individual terms (obtained from the energy set)
		

		System * pSys;
		EnergySet * pESet;
		RandomNumberGenerator * pRng;

		std::map<std::string, std::vector<Interaction*> > * pEnergyTerms;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, std::vector<Interaction*> > > > > > subdividedInteractions;
		std::vector<Position*> variablePositions;
		//std::vector<vector<Position*> > slavePositions; // these positions are linked by symmetry to the variable
		std::vector<std::vector<Residue*> > variableIdentities;
		std::vector<std::vector<std::vector<Residue*> > > slaveIdentities;

		double fixE;
		std::vector<std::vector<double> > selfE;
		std::vector<std::vector<std::vector<std::vector<double> > > > pairE;
		std::vector<std::vector<std::vector<std::vector<bool> > > > pairEFlag; // true if the energy is computed already and available

		bool saveEbyTerm;
		std::map<std::string, double> fixEbyTerm;
		std::vector<std::vector<std::map<std::string, double> > > selfEbyTerm;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, double> > > > > pairEbyTerm;

		bool saveInteractionCount;
		unsigned int fixCount;
		std::vector<std::vector<unsigned int> > selfCount;
		std::vector<std::vector<std::vector<std::vector<unsigned int> > > > pairCount;

		std::map<std::string, unsigned int> fixCountByTerm;
		std::vector<std::vector<std::map<std::string, unsigned int> > > selfCountByTerm;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, unsigned int> > > > > pairCountByTerm;

		std::vector<unsigned int> variableCount;
		//std::map<Position*, bool> variablePosMap;
		std::map<Position*, unsigned int> variablePosIndex;
		//std::map<Position*, unsigned int> slavePosIndex;
		//std::map<std::vector<unsigned int>, unsigned int> IdRotAbsIndex;

		std::vector<std::vector<std::string> > rotamerDescriptors;
		std::vector<std::vector<std::vector<unsigned int> > > rotamerPos_Id_Rot;

		// Side Chain Optimization Functions
		bool runDEE;
		bool runUnbiasedMC;
		bool runSCMFBiasedMC;
		bool runSCMF;
		bool runEnum;
		bool verbose;

		bool onTheFly; // if true, pair energies are not precomputed

		std::vector<std::vector<unsigned int> > aliveRotamers;
		std::vector<std::vector<bool> > aliveMask;
		std::vector<unsigned int> mostProbableSCMFstate;
		std::vector<unsigned int> bestSCMFBiasedMCstate;
		std::vector<unsigned int> bestUnbiasedMCstate;

		vector<double> minBound;
		vector<vector<unsigned int> > minStates;

		void saveMin(double _boundE, vector<unsigned int> _stateVec, int _maxSaved); 
		int selectRandomStateAtPosition(int _position,vector<unsigned int>& _currentState) const ;
		void getRandomState(vector<unsigned int>& _currentState ) ;

		// DEE Options
		double DEEenergyOffset;
		bool DEEdoSimpleGoldsteinSingle;
		bool DEEdoSimpleGoldsteinPair;

		// Enumeration Options
		int enumerationLimit;

		// SCMF Options
		int maxSavedResults;
		int SCMFtemperature;
		int SCMFcycles;

		// MCO Options
		double mcStartT;
		double mcEndT;
		int mcCycles;
		int mcShape;
		int mcMaxReject;
		int mcDeltaSteps;
		double mcMinDeltaE;
		

};

inline void SelfPairManager::setSystem(System * _pSystem) {
	pSys = _pSystem;
	pESet = _pSystem->getEnergySet();
	pEnergyTerms = pESet->getEnergyTerms();
	_pSystem->updateVariablePositions();
	findVariablePositions();
	subdivideInteractions();
	updateWeights();
}
inline void SelfPairManager::updateWeights() {
	weights = pESet->getWeightMap();
}
inline System * SelfPairManager::getSystem() const { return pSys;}

inline unsigned int SelfPairManager::getNumberOfVariablePositions() const {return variableCount.size();}
inline std::vector<unsigned int> SelfPairManager::getNumberOfRotamers() const {return variableCount;}


inline std::vector<std::string> SelfPairManager::getStateDescriptors(std::vector<unsigned int> _overallRotamerStates) const {
	std::vector<std::string> out;
	for (unsigned int i=0; i<_overallRotamerStates.size(); i++) {
		out.push_back(rotamerDescriptors[i][_overallRotamerStates[i]]);
	}
	return out;
}
inline std::vector<std::vector<unsigned int> > SelfPairManager::getStatePositionIdentityRotamerIndeces(std::vector<unsigned int> _overallRotamerStates) const {
	std::vector<std::vector<unsigned int> > out;
	for (unsigned int i=0; i<_overallRotamerStates.size(); i++) {
		out.push_back(rotamerPos_Id_Rot[i][_overallRotamerStates[i]]);
	}
	return out;
}

inline void SelfPairManager::setMCOptions(double _startT, double _endT, int _nCycles, int _shape, int _maxReject, int _deltaSteps, double _minDeltaE) {
	mcStartT = _startT;
	mcEndT = _endT;
	mcCycles = _nCycles;
	mcShape  = _shape;
	mcMaxReject = _maxReject;
	mcDeltaSteps = _deltaSteps;
	mcMinDeltaE = _minDeltaE;

}
inline void SelfPairManager::setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}
inline RandomNumberGenerator * SelfPairManager::getRandomNumberGenerator() const {
	return pRng;
}

inline void SelfPairManager::seed(unsigned int _seed){
	if (_seed == 0) {
		pRng->setTimeBasedSeed();
	}
	else {
		pRng->setSeed(_seed);
	}
}


inline unsigned int SelfPairManager::getSeed() const {
	return pRng->getSeed();
}

inline void SelfPairManager::saveEnergiesByTerm(bool _save) {
	saveEbyTerm = _save;
}
inline void SelfPairManager::saveInteractionCounts(bool _save) {
	saveInteractionCount = _save;
}
inline void SelfPairManager::setOnTheFly(bool _onTheFly) {
	onTheFly = _onTheFly;
}
inline void SelfPairManager::setEnumerationLimit(int _enumLimit) {
	enumerationLimit = _enumLimit;
}

inline int SelfPairManager::getNumPositions() { 
	return selfE.size();
}

inline int SelfPairManager::getNumRotamers(int _position) {

	if ( _position >= selfE.size()) { 
		return 0; 
	} 

	return selfE[_position].size();
}


}

#endif
