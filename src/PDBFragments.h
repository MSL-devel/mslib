/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

#ifndef PDBFRAGMENTS_H
#define PDBFRAGMENTS_H

#include <string>

#include "AtomPointerVector.h"
#include "System.h"

namespace MSL { 
class PDBFragments{

	public:
		PDBFragments();
		PDBFragments(std::string _fragDbFile, std::string _BBQTableForBackboneAtoms="");
		~PDBFragments();


		void loadFragmentDatabase();

		
		int searchForMatchingFragments(System &_sys, std::vector<std::string> &_stemResidues,int _numResiduesInFragment=-1);

		System & getSystem();
		AtomPointerVector & getAtomPointers();
		enum dbAtoms { caOnly=0, allAtoms=1 };
		void printMe();

	private:
		std::string fragDbFile;
		dbAtoms fragType;
		std::string bbqTable;
		AtomPointerVector fragDB;

		System *lastResults;
};

inline PDBFragments::PDBFragments() { 	fragType   = caOnly; lastResults = NULL; }
inline PDBFragments::PDBFragments(std::string _fragDbFile,std::string _BBQTableForBackboneAtoms) {
	fragDbFile = _fragDbFile;
	lastResults = NULL;

	if (_BBQTableForBackboneAtoms == ""){
		fragType   = allAtoms;
	} else {
		fragType   = caOnly;
	}
	bbqTable = _BBQTableForBackboneAtoms;

}
inline PDBFragments::~PDBFragments() {
	if (lastResults != NULL){
		delete(lastResults);
	}
}
inline void PDBFragments::loadFragmentDatabase(){
	fragDB.load_checkpoint(fragDbFile);

}

inline System & PDBFragments::getSystem(){
	return (*lastResults);
}

inline AtomPointerVector & PDBFragments::getAtomPointers(){
	return (*lastResults).getAllAtomPointers();
}

inline void PDBFragments::printMe(){
	for (uint i = 0; i < fragDB.size();i++){
		std::cout << fragDB(i).toString()<<fragDB(i).getSegID()<<std::endl;
	}
}
}

#endif
