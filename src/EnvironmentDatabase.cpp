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


#include "EnvironmentDatabase.h"
#include "CartesianGeometry.h"
#include "PDBWriter.h"
#include "AtomPointerVector.h"
#include "AtomSelection.h"

using namespace MSL;
using namespace std;


EnvironmentDatabase::EnvironmentDatabase(){
	type = "";
	lastSearchedKey ="";
}

EnvironmentDatabase::EnvironmentDatabase(string _type){
	type = _type;
}

EnvironmentDatabase::EnvironmentDatabase(EnvironmentDatabase &_ed){
	copy(_ed);
}

EnvironmentDatabase::~EnvironmentDatabase(){
	
	for (uint i = 0; i < descriptors.size();i++){
		delete(descriptors[i]);
	}
	descriptors.clear();
	lookupTable.clear();
}

void EnvironmentDatabase::operator=(EnvironmentDatabase &_ed){
	copy(_ed);
}

void EnvironmentDatabase::createDatabase(System &_sys, string _systemName){


	for (uint i = 0; i < _sys.positionSize();i++){
		Residue res = _sys.getResidue(i);

		if (!res.atomExists("CA") || !res.atomExists("N") || !res.atomExists("C")) continue; // Skip non-amino acids

		if (type == "" || res.getResidueName() == type){
			cout << "Working on "<<res.getChainId()<<" "<<res.getResidueName()<<" "<<res.getResidueNumber()<<endl;


			// Add refFrame
 			Frame refFrame;
 			refFrame.setName("ref");
 			if (res.getResidueName() != "GLY"){
				if (!res.atomExists("CB")) continue;
 				refFrame.computeFrameFrom3Atoms(res("N"), res("CA"), res("CB"));
 			} else {
 				CartesianPoint CBcoor = CartesianGeometry::build(res("CA").getCoor(), res("N").getCoor(), res("C").getCoor(), 1.521, 110.5, -122.5);
 				Atom CB;
 				CB.setCoor(CBcoor);
 				refFrame.computeFrameFrom3Atoms(res("N"), res("CA"), CB);
 			}



			// Add environment..
			int startRes = res.getResidueNumber() - 7;
			if (startRes < 0){
				startRes = 0;
			}
			int endRes = res.getResidueNumber() + 7;
			if (endRes > _sys(res.getChainId()).positionSize()){
				endRes = _sys(res.getChainId()).positionSize();
			}

 			char a[200];
			sprintf(a,"%10d, name CA and not resi %d-%d WITHIN %d OF ((chain %s and resi %d) and name CA)",i,startRes,endRes,10,res.getChainId().c_str(),res.getResidueNumber());

			cout << "\tUsing environment selection: "<<(string)a<<endl;
 			AtomSelection sel(_sys.getAtomPointers());
 			AtomPointerVector env = sel.select(string(a));

			// Bail out and don't add if less than 3 residues in environment, most likely a surface residue.
			if (env.size() < 3) continue;

			// Add core atoms
			AtomPointerVector ats = res.getAtomPointers();
			descriptors.push_back(new EnvironmentDescriptor);
			descriptors.back()->setName(_systemName);
			descriptors.back()->setCore(ats);			
 			descriptors.back()->setReferenceFrame(refFrame);
 			descriptors.back()->setEnvironment("CA", env);


			// Add some meta data for easy lookup.
			// For instance the diagnol of the angle matrix.
			// Bin angle in 10 degree bins?


			/*
			  Key lookup: 
			    DistanceBin:EigenVector1Bin:EigenVector2Bin:EigenVector3Bin
			*/

			string key = descriptors.back()->generateLookupKey("CA");

			//vector<EnvironmentDescriptor*> edList;
			//map<string,vector<EnvironmentDescriptor*> >::iterator it;
			//it = lookupTable.find(key);
			//cout << "\tAdding Meta Data (Key): "<<key<<endl;
			lookupTable[key].push_back(descriptors.back());			
		} 
		
	}
	

}


bool EnvironmentDatabase::searchForEnvironment(EnvironmentDescriptor &_ed, string _envType){


	lastSearchedKey = _ed.generateLookupKey(_envType);
	map<string,vector<EnvironmentDescriptor*> >::iterator it;
	it = lookupTable.find(lastSearchedKey);

	
	if (it != lookupTable.end()){
		//cout << "Key Found is: "<<lastSearchedKey<<endl;
		return true;
	}


	return false;
}


vector<EnvironmentDescriptor*>& EnvironmentDatabase::getSearchResults(){

	return lookupTable[lastSearchedKey];
}

vector<EnvironmentDescriptor*>& EnvironmentDatabase::getAllDescriptors(){

	return descriptors;
}

void EnvironmentDatabase::copy(EnvironmentDatabase &_ed){

	vector<EnvironmentDescriptor *> allDescriptors = _ed.getAllDescriptors();

	for (uint i = 0; i < allDescriptors.size();i++){
		descriptors.push_back(new EnvironmentDescriptor(*allDescriptors[i]));

		string key = descriptors.back()->generateLookupKey("CA");
		lookupTable[key].push_back(descriptors.back());			
	}
	

}
