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

#ifndef PAIRWISEENERGYCALCULATOR
#define PAIRWISEENERGYCALCULATOR


#include "System.h"
#include "AtomicPairwiseEnergy.h"
#include "TwoBodyDistanceDependentPotentialTable.h"
#include <vector>

namespace MSL { 
class PairwiseEnergyCalculator {
	public:

		PairwiseEnergyCalculator(std::string _charmmParamterFile);
		~PairwiseEnergyCalculator();


		// Different ways to calculate energy
		double calculateTotalEnergy(System &_sys);
		double calculateStateEnergy(System &_sys, std::vector<unsigned int> &_stateVector);
		double calculateEnergyTable(System &_sys);

		double calculateTotalEnergy(System &_sys, TwoBodyDistanceDependentPotentialTable & tbd);

		double getStateEnergy(System &_sys, std::vector<unsigned int> &_stateVector);
		void printSummary();
		//void getSummaryByGroup();
		void printPairwiseTable();

		// After calling calculateEnergyTable, these will be meaningful.
		std::vector<std::vector<double> > & getSelfTable();
		std::vector<std::vector<double> > & getTemplateTable();
		std::vector<std::vector<std::vector<std::vector<double> > > > & getPairTable();

		// Only in debug mode, will get Interactions from AtomicPairwiseEnergy
		std::map<Interaction*,int> & getInteractions();

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

		std::vector<std::vector<double> > selfEnergy;                
		std::vector<std::vector<double> > templateEnergy;
		std::vector<std::vector<std::vector<std::vector<double> > > >  pairEnergy;



		

};
inline std::map<Interaction*,int> & PairwiseEnergyCalculator::getInteractions(){return pairwiseEnergy->getInteractions();}
inline std::vector<std::vector<double> > & PairwiseEnergyCalculator::getSelfTable() { return selfEnergy; }
inline std::vector<std::vector<double> > & PairwiseEnergyCalculator::getTemplateTable() { return templateEnergy; }
inline std::vector<std::vector<std::vector<std::vector<double> > > >  & PairwiseEnergyCalculator::getPairTable() { return pairEnergy; }
inline AtomicPairwiseEnergy& PairwiseEnergyCalculator::getAtomicPairwiseEnergyObject() { return *pairwiseEnergy; }

}

#endif
