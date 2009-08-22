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

#ifndef PAIRWISEENERGYCALCULATOR
#define PAIRWISEENERGYCALCULATOR


#include "System.h"
#include "AtomicPairwiseEnergy.h"
#include <vector>

class PairwiseEnergyCalculator {
	public:

		PairwiseEnergyCalculator(string _charmmParamterFile);
		~PairwiseEnergyCalculator();


		// Different ways to calculate energy
		double calculateTotalEnergy(System &_sys);
		double calculateStateEnergy(System &_sys, vector<unsigned int> &_stateVector);
		double calculateEnergyTable(System &_sys);


		double getStateEnergy(System &_sys, vector<unsigned int> &_stateVector);
		void printSummary();
		//void getSummaryByGroup();
		void printPairwiseTable();

		// After calling calculateEnergyTable, these will be meaningful.
		vector<vector<double> > & getSelfTable();
		vector<vector<double> > & getTemplateTable();
		vector<vector<vector<vector<double> > > > & getPairTable();

		// Only in debug mode, will get Interactions from AtomicPairwiseEnergy
		map<Interaction*,int> & getInteractions();

		AtomicPairwiseEnergy& getAtomicPairwiseEnergyObject();

	private:

		PairwiseEnergyCalculator();

		// Private object to do atom-atom energy calculations.
		AtomicPairwiseEnergy *pairwiseEnergy;


		/*
		  Indices are: 
		       position,conformation
		               OR
		       position,conformation,position,conformation
			       
		*/

		vector<vector<double> > selfEnergy;                
		vector<vector<double> > templateEnergy;
		vector<vector<vector<vector<double> > > >  pairEnergy;



		

};
inline map<Interaction*,int> & PairwiseEnergyCalculator::getInteractions(){return pairwiseEnergy->getInteractions();}
inline vector<vector<double> > & PairwiseEnergyCalculator::getSelfTable() { return selfEnergy; }
inline vector<vector<double> > & PairwiseEnergyCalculator::getTemplateTable() { return templateEnergy; }
inline vector<vector<vector<vector<double> > > >  & PairwiseEnergyCalculator::getPairTable() { return pairEnergy; }
inline AtomicPairwiseEnergy& PairwiseEnergyCalculator::getAtomicPairwiseEnergyObject() { return *pairwiseEnergy; }

#endif
