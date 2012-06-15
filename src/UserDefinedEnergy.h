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

#ifndef USERDEFINEDENERGY_H
#define USERDEFINEDENERGY_H

#include "TwoBodyDistanceDependentPotentialTable.h"


namespace MSL { 
class UserDefinedEnergy {
	public:
	        static UserDefinedEnergy * instance();

		void addTwoBodyPotential(TwoBodyDistanceDependentPotentialTable *_p);
		double getTwoBodyPotentialValue(std::string _potentialName, std::string _name1, std::string _name2, double _dist) ;
	
		bool potentialExists(std::string _potentialName);

		TwoBodyDistanceDependentPotentialTable& getTwoBodyPotentialTable(std::string _potentialName);
		//std::map<std::string, TwoBodyDistanceDependentPotentialTable*> &getAllTwoBodyPotentialTables;

	protected:
		UserDefinedEnergy();
		UserDefinedEnergy(const UserDefinedEnergy & _instance);
		void operator=(const UserDefinedEnergy & _instance);


		std::map<std::string, TwoBodyDistanceDependentPotentialTable*> twoBodyPotentials;
};

/* INLINES */


inline UserDefinedEnergy * UserDefinedEnergy::instance(){
	static UserDefinedEnergy inst;
	return &inst;
}

inline UserDefinedEnergy::UserDefinedEnergy(){}
inline UserDefinedEnergy::UserDefinedEnergy(const UserDefinedEnergy &_instance){}
inline TwoBodyDistanceDependentPotentialTable& UserDefinedEnergy::getTwoBodyPotentialTable(std::string _potentialName){

	std::map<std::string,TwoBodyDistanceDependentPotentialTable*>::iterator it;
	it = twoBodyPotentials.find(_potentialName);
	if (it == twoBodyPotentials.end()){
		std::cerr << "ERROR 6254 UserDefinedEnergy::getTwoBodyPotentialTable(..) potential not found: "<<_potentialName<<std::endl;
		for (it = twoBodyPotentials.begin();it != twoBodyPotentials.end();it++){
			std::cerr << "\t Pot = "<<it->first<<std::endl;
		}
		exit(6254);
	}

	return (*(it->second));
	
}

inline bool UserDefinedEnergy::potentialExists(std::string _potentialName){

	std::map<std::string,TwoBodyDistanceDependentPotentialTable*>::iterator it;
	it = twoBodyPotentials.find(_potentialName);
	if (it == twoBodyPotentials.end()){
		std::cerr << "WARNING UserDefinedEnergy::getTwoBodyPotentialTable(..) potential not found: "<<_potentialName<<std::endl;
		for (it = twoBodyPotentials.begin();it != twoBodyPotentials.end();it++){
			std::cerr << "\t Pot = "<<it->first<<std::endl;
		}
		return false;
	}

	return true;
	
}
}

#endif
