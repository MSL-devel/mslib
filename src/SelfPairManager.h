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

#include "System.h"
#include "EnergySet.h"

using namespace std;

class SelfPairManager {
	public:
		SelfPairManager();
		SelfPairManager(System * _pSystem);
		~SelfPairManager();

		void calculateEnergies();

		//  GET THE ENERGIES
		double getStateEnergy(vector<unsigned int> _overallRotamerStates, string _term="");
		double getStateEnergy(vector<unsigned int> _residueStates, vector<unsigned int> _rotamerStates);
		double getStateEnergy(vector<string> _residueNames, vector<unsigned int> _rotamerStates);

		unsigned int getStateInteractionCount(vector<unsigned int> _overallRotamerStates, string _term="");

		void setSystem(System * _pSystem);
		System * getSystem() const;

		unsigned int getNumberOfVariablePositions() const;
		vector<unsigned int> getNumberOfRotamers() const;

		string getSummary(vector<unsigned int> _overallRotamerStates);
		vector<string> getStateDescriptors(vector<unsigned int> _overallRotamerStates) const;
		vector<vector<unsigned int> > getStatePositionIdentityRotamerIndeces(vector<unsigned int> _overallRotamerStates) const;
		

	private:
		void setup();
		void copy(const SelfPairManager & _sysBuild);
		void deletePointers();
		void findVariablePositions();
		void subdivideInteractions();

		System * pSys;
		EnergySet * pESet;

		map<string, vector<Interaction*> > * pEnergyTerms;
		vector<vector<vector<vector<map<string, vector<Interaction*> > > > > > subdividedInteractions;
		vector<Position*> variablePositions;
		vector<vector<Residue*> > variableIdentities;

		double fixE;
		vector<vector<double> > selfE;
		vector<vector<vector<vector<double> > > > pairE;

		map<string, double> fixEbyTerm;
		vector<vector<map<string, double> > > selfEbyTerm;
		vector<vector<vector<vector<map<string, double> > > > > pairEbyTerm;

		unsigned int fixCount;
		vector<vector<unsigned int> > selfCount;
		vector<vector<vector<vector<unsigned int> > > > pairCount;

		map<string, unsigned int> fixCountByTerm;
		vector<vector<map<string, unsigned int> > > selfCountByTerm;
		vector<vector<vector<vector<map<string, unsigned int> > > > > pairCountByTerm;

		vector<unsigned int> variableCount;
		//map<Position*, bool> variablePosMap;
		map<Position*, unsigned int> variablePosIndex;
		map<vector<unsigned int>, unsigned int> IdRotAbsIndex;

		vector<vector<string> > rotamerDescriptors;
		vector<vector<vector<unsigned int> > > rotamerPos_Id_Rot;
};

inline void SelfPairManager::setSystem(System * _pSystem) {
	pSys = _pSystem;
	pESet = _pSystem->getEnergySet();
	pEnergyTerms = pESet->getEnergyTerms();
	findVariablePositions();
	subdivideInteractions();
	_pSystem->updateVariablePositions();
}
inline System * SelfPairManager::getSystem() const { return pSys;}

inline unsigned int SelfPairManager::getNumberOfVariablePositions() const {return variableCount.size();}
inline vector<unsigned int> SelfPairManager::getNumberOfRotamers() const {return variableCount;}


inline vector<string> SelfPairManager::getStateDescriptors(vector<unsigned int> _overallRotamerStates) const {
	vector<string> out;
	for (unsigned int i=0; i<_overallRotamerStates.size(); i++) {
		out.push_back(rotamerDescriptors[i][_overallRotamerStates[i]]);
	}
	return out;
}
inline vector<vector<unsigned int> > SelfPairManager::getStatePositionIdentityRotamerIndeces(vector<unsigned int> _overallRotamerStates) const {
	vector<vector<unsigned int> > out;
	for (unsigned int i=0; i<_overallRotamerStates.size(); i++) {
		out.push_back(rotamerPos_Id_Rot[i][_overallRotamerStates[i]]);
	}
	return out;
}

#endif
