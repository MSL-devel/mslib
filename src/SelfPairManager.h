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


namespace MSL { 
class SelfPairManager {
	public:
		SelfPairManager();
		SelfPairManager(System * _pSystem);
		~SelfPairManager();

		void calculateEnergies();

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
		

	private:
		void setup();
		void copy(const SelfPairManager & _sysBuild);
		void deletePointers();
		void findVariablePositions();
		void subdivideInteractions();

		System * pSys;
		EnergySet * pESet;

		std::map<std::string, std::vector<Interaction*> > * pEnergyTerms;
		std::vector<std::vector<std::vector<std::vector<std::map<std::string, std::vector<Interaction*> > > > > > subdividedInteractions;
		std::vector<Position*> variablePositions;
		std::vector<std::vector<Residue*> > variableIdentities;

		double fixE;
		std::vector<std::vector<double> > selfE;
		std::vector<std::vector<std::vector<std::vector<double> > > > pairE;

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
		std::map<std::vector<unsigned int>, unsigned int> IdRotAbsIndex;

		std::vector<std::vector<std::string> > rotamerDescriptors;
		std::vector<std::vector<std::vector<unsigned int> > > rotamerPos_Id_Rot;
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

}

#endif
