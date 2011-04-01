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

#ifndef SELFPARIMANAGER_H
#define SELFPARIMANAGER_H

#include <iostream>
#include <vector>

#include "DeadEndElimination.h"
#include "Enumerator.h"
#include "MonteCarloManager.h"
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
		~SelfPairManager();

		void calculateEnergies();

		void setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG);
		RandomNumberGenerator * getRandomNumberGenerator() const;

		void seed(unsigned int _seed); // use 0 for time based seed
		unsigned int getSeed() const; // get back the seed (useful to get the time based one)
		//  GET THE ENERGIES
		double getStateEnergy(std::vector<unsigned int> _overallRotamerStates, std::string _term="");
		double getStateEnergy(std::vector<unsigned int> _residueStates, std::vector<unsigned int> _rotamerStates);
		double getStateEnergy(std::vector<std::string> _residueNames, std::vector<unsigned int> _rotamerStates);

		unsigned int getStateInteractionCount(std::vector<unsigned int> _overallRotamerStates, std::string _term="");

		void setSystem(System * _pSystem);
		System * getSystem() const;

		unsigned int getNumberOfVariablePositions() const;
		std::vector<unsigned int> getNumberOfRotamers() const;

		std::string getSummary(std::vector<unsigned int> _overallRotamerStates);
		std::vector<std::string> getStateDescriptors(std::vector<unsigned int> _overallRotamerStates) const;
		std::vector<std::vector<unsigned int> > getStatePositionIdentityRotamerIndeces(std::vector<unsigned int> _overallRotamerStates) const;
		
		void saveEnergiesByTerm(bool _save); // false by default

		double getFixedEnergy() const;

		// gets the reduced self energy table if cutoff is applied
		std::vector<std::vector<double> > & getSelfEnergy(); 		

		// gets the reduced pair energy table if cutoff is applied
		std::vector<std::vector<std::vector<std::vector<double> > > > & getPairEnergy(); 

		// Side Chain Optimization Functions
		void setRunDEE(bool _toogle);
		void setRunMC(bool _toogle);
		void setRunSCMF(bool _toogle);
		void setRunEnum(bool _toogle);

		// reject rotamers that have a high selfEnergy with respect to the best rotamer in that position
		void setSelfECutoff(double _cutoff);
		void setUseSelfECutoff(bool _use);
		std::vector<std::vector<bool> > getSelfEnergyMask() const;

		void setVerbose(bool _toggle);

		void runOptimizer();

		std::vector<double> getMinBound();
		std::vector<vector<unsigned int> > getMinStates();

		std::vector<std::vector<unsigned int> > getDEEAliveRotamers();  // after applying cutoff
		std::vector<std::vector<bool> > getDEEAliveMask(); // after applying cutoff

		std::vector<unsigned int> getSCMFstate();
		std::vector<unsigned int> getMCOstate();
		std::vector<unsigned int> getMCfinalState();

	private:
		void setup();
		void copy(const SelfPairManager & _sysBuild);
		void deletePointers();
		void findVariablePositions();
		void subdivideInteractions();
		bool deleteRng;

		System * pSys;
		EnergySet * pESet;
		RandomNumberGenerator * pRng;

		std::map<std::string, std::vector<Interaction*> > * pEnergyTerms;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, std::vector<Interaction*> > > > > > subdividedInteractions;
		std::vector<Position*> variablePositions;
		std::vector<std::vector<Residue*> > variableIdentities;

		double fixE;
		std::vector<std::vector<double> > selfE;
		std::vector<std::vector<std::vector<std::vector<double> > > > pairE;

		bool saveEbyTerm;
		std::map<std::string, double> fixEbyTerm;
		std::vector<std::vector<std::map<std::string, double> > > selfEbyTerm;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, double> > > > > pairEbyTerm;

		unsigned int fixCount;
		std::vector<std::vector<unsigned int> > selfCount;
		std::vector<std::vector<std::vector<std::vector<unsigned int> > > > pairCount;

		std::map<std::string, unsigned int> fixCountByTerm;
		std::vector<std::vector<std::map<std::string, unsigned int> > > selfCountByTerm;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, unsigned int> > > > > pairCountByTerm;

		std::vector<unsigned int> variableCount;
		//std::map<Position*, bool> variablePosMap;
		std::map<Position*, unsigned int> variablePosIndex;
		//std::map<std::vector<unsigned int>, unsigned int> IdRotAbsIndex;

		std::vector<std::vector<std::string> > rotamerDescriptors;
		std::vector<std::vector<std::vector<unsigned int> > > rotamerPos_Id_Rot;

		// Side Chain Optimization Functions
		bool runDEE;
		bool runMC;
		bool runSCMF;
		bool runEnum;
		bool verbose;

		bool useSelfECutoff; // false by default
		double selfECutoff; // 15.0 by default
		std::vector<std::vector<bool> > selfEMask;
		std::vector<unsigned int> getOriginalRotamerState(vector<unsigned int>& _reducedState );
		std::vector<unsigned int> getReducedRotamerState(vector<unsigned int>& _origState );
		double getInternalStateEnergy(std::vector<unsigned int> _overallRotamerStates, std::string _term="");
		unsigned int getInternalStateInteractionCount(std::vector<unsigned int> _overallRotamerStates, std::string _term="");

		std::vector<std::vector<unsigned int> > aliveRotamers;
		std::vector<std::vector<bool> > aliveMask;
		std::vector<unsigned int> mostProbableSCMFstate;
		std::vector<unsigned int> finalMCOstate;

		vector<double> minBound;
		vector<vector<unsigned int> > minStates;

		void saveMin(double _boundE, vector<unsigned int> _stateVec, int _maxSaved); 

		// DEE Options
		double DEEenergyOffset;
		bool DEEdoSimpleGoldsteinSingle;

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

inline void SelfPairManager::setRandomNumberGenerator(RandomNumberGenerator * _pExternalRNG) {
	if (deleteRng == true) {
		delete pRng;
	}
	pRng = _pExternalRNG;
	deleteRng = false;
}

inline std::vector<std::vector<bool> > SelfPairManager::getSelfEnergyMask() const {
	return selfEMask;
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

inline void SelfPairManager::setUseSelfECutoff(bool _use) {
	useSelfECutoff = _use;
}

inline void SelfPairManager::setSelfECutoff(double _cutoff) {
	selfECutoff = _cutoff;
}
}

#endif
