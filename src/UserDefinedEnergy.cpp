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

#include "UserDefinedEnergy.h"

using namespace MSL;
using namespace std;


void UserDefinedEnergy::operator= (const UserDefinedEnergy &_instance){

}

void UserDefinedEnergy::addTwoBodyPotential(TwoBodyDistanceDependentPotentialTable *_p){
	twoBodyPotentials[_p->getPotentialName()] = _p;
}

double UserDefinedEnergy::getTwoBodyPotentialValue(string _potentialName, string _name1, string _name2, double _dist)  {

	map<string, TwoBodyDistanceDependentPotentialTable *>::iterator it;
	it = twoBodyPotentials.find(_potentialName);

	// Check if we have the potential asked for..
	if (it == twoBodyPotentials.end()){
		cerr << "ERROR 2342 UserDefinedEnergy::getTwoBodyPotentialValue(..) potential not defined: "<<_potentialName<<endl;
		for (it = twoBodyPotentials.begin();it != twoBodyPotentials.end();it++){
			cerr << "\t Pot = "<<it->first<<endl;
		}
		exit(2342);
		
	}

	// Get Distance bin index
	int bin = (it->second)->getBin(_dist);

	// Make sure valid distance bin..
	if (bin == -1){
		cerr << "ERROR 2343 UserDefinedEnergy::getTwoBodyPotentialValue(..) dist not defined by a bin: "<<_dist<<endl; 
		exit(2343);
	}

	// Return the minimum if below cutoff
	if (_dist < (it->second)->getMinDistCutoff()){
		return (it->second)->getValueBelowCutoff();
	}

	// Return the maximum if above cutoff
	if (_dist > (it->second)->getMaxDistCutoff()){
		return (it->second)->getValueAboveCutoff();
	}


	
	return (it->second)->getPotential(_name1,_name2,bin);
	
}
