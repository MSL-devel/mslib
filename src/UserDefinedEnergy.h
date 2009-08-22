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

#ifndef USERDEFINEDENERGY_H
#define USERDEFINEDENERGY_H

#include "TwoBodyDistanceDependentPotentialTable.h"

using namespace std;

class UserDefinedEnergy {
	public:
	        static UserDefinedEnergy * instance();

		void addTwoBodyPotential(TwoBodyDistanceDependentPotentialTable *_p);
		double getTwoBodyPotentialValue(string _potentialName, string _name1, string _name2, double _dist) ;
	
		bool potentialExists(string _potentialName);

		TwoBodyDistanceDependentPotentialTable& getTwoBodyPotentialTable(string _potentialName);
		//map<string, TwoBodyDistanceDependentPotentialTable*> &getAllTwoBodyPotentialTables;

	protected:
		UserDefinedEnergy();
		UserDefinedEnergy(const UserDefinedEnergy & _instance);
		void operator=(const UserDefinedEnergy & _instance);


		map<string, TwoBodyDistanceDependentPotentialTable*> twoBodyPotentials;
};

/* INLINES */


inline UserDefinedEnergy * UserDefinedEnergy::instance(){
	static UserDefinedEnergy inst;
	return &inst;
}

inline UserDefinedEnergy::UserDefinedEnergy(){}
inline UserDefinedEnergy::UserDefinedEnergy(const UserDefinedEnergy &_instance){}
inline TwoBodyDistanceDependentPotentialTable& UserDefinedEnergy::getTwoBodyPotentialTable(string _potentialName){

	map<string,TwoBodyDistanceDependentPotentialTable*>::iterator it;
	it = twoBodyPotentials.find(_potentialName);
	if (it == twoBodyPotentials.end()){
		cerr << "ERROR 6254 UserDefinedEnergy::getTwoBodyPotentialTable(..) potential not found: "<<_potentialName<<endl;
		for (it = twoBodyPotentials.begin();it != twoBodyPotentials.end();it++){
			cerr << "\t Pot = "<<it->first<<endl;
		}
		exit(6254);
	}

	return (*(it->second));
	
}

inline bool UserDefinedEnergy::potentialExists(string _potentialName){

	map<string,TwoBodyDistanceDependentPotentialTable*>::iterator it;
	it = twoBodyPotentials.find(_potentialName);
	if (it == twoBodyPotentials.end()){
		cerr << "WARNING UserDefinedEnergy::getTwoBodyPotentialTable(..) potential not found: "<<_potentialName<<endl;
		for (it = twoBodyPotentials.begin();it != twoBodyPotentials.end();it++){
			cerr << "\t Pot = "<<it->first<<endl;
		}
		return false;
	}

	return true;
	
}
#endif
