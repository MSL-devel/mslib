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

#include "UserDefinedEnergySetBuilder.h"


UserDefinedEnergySetBuilder::UserDefinedEnergySetBuilder(){}

void UserDefinedEnergySetBuilder::setup(){
}

void UserDefinedEnergySetBuilder::addEnergyInteractions(System &_sys){
	
	
}

void UserDefinedEnergySetBuilder::addTwoBodyInteractions(System &_sys, string _filename){
	
	TwoBodyDistanceDependentPotentialTable *twoBody = new TwoBodyDistanceDependentPotentialTable();
	twoBody->readPotentialTable(_filename);


	UserDefinedEnergy::instance()->addTwoBodyPotential(twoBody);


	/*
	  For each position i , each position j where abs(i - j) >= skipping number
	               Each identity of i, each identity of j
		             All atoms pos i, id i  vs All atoms pos j, id j
	 */

	//map<string,TwoBodyDistanceDependentPotentialTable*> tbl = UserDefinedEnergy::instance()->getAllTwoBodyPotentialTables();
	for (uint posI = 0 ; posI < _sys.positionSize();posI++){
		
		for (uint id1 = 0; id1 < _sys.getPosition(posI).size();id1++){

			for (uint posJ = 0 ; posJ < _sys.positionSize();posJ++){


				// Skipping Number
				if ( _sys.getPosition(posI).getChainId() == _sys.getPosition(posJ).getChainId() && 
				     abs(_sys.getPosition(posI).getResidueNumber() - _sys.getPosition(posJ).getResidueNumber()) < twoBody->getResidueSkippingNumber()) {
					continue;
				}

				for (uint id2 = 0; id2 < _sys.getPosition(posJ).size();id2++){


					// Now loop over all atom-atom interactions...
					for (uint a1 = 0; a1 < _sys.getPosition(posI)(id1).size();a1++){
						for (uint a2 = 0; a2 < _sys.getPosition(posJ)(id2).size();a2++){


							UserDefinedInteraction *udInt = new UserDefinedInteraction();
							udInt->setAtoms(_sys.getPosition(posI)(id1).getAtom(a1),_sys.getPosition(posJ)(id2).getAtom(a2));
							udInt->setName(twoBody->getPotentialName());

							_sys.getEnergySet()->addInteraction(udInt);
							
						}
					}
				}
			}
		}

	}
}
