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

#ifndef PDBFRAGMENTS_H
#define PDBFRAGMENTS_H

#include <string>

#include "AtomVector.h"
#include "System.h"
#include "Chain.h"

class PDBFragments{

	public:
		PDBFragments();
		PDBFragments(string _fragDbFile);
		~PDBFragments();


		void loadFragmentDatabase();

		
		int searchForMatchingFragments(Chain &_ch, vector<int> &_stemResidueIndices,int _numResiduesInFragment=-1);

		System & getLastSearchResults();

		enum dbAtoms { caOnly=0, allAtoms=1 };
		void printMe();

	private:
		string fragDbFile;
		dbAtoms fragType;
		AtomVector fragDB;

		System *lastResults;
};

inline PDBFragments::PDBFragments() { 	fragType   = caOnly; lastResults = NULL;}
inline PDBFragments::PDBFragments(string _fragDbFile) {
	fragDbFile = _fragDbFile;
	fragType   = caOnly;
	lastResults = NULL;
}
inline PDBFragments::~PDBFragments() {
	if (lastResults != NULL){
		delete(lastResults);
	}
}
inline void PDBFragments::loadFragmentDatabase(){
	fragDB.load_checkpoint(fragDbFile);

}

inline System & PDBFragments::getLastSearchResults(){
	return (*lastResults);
}

inline void PDBFragments::printMe(){
	for (uint i = 0; i < fragDB.size();i++){
		cout << fragDB(i).toString()<<fragDB(i).getSegID()<<endl;
	}
}
#endif
