#include "VectorHashing.h"


#include "MslOut.h"
static MslOut MSLOUT("VectorHashing");

VectorHashing::VectorHashing(){
	archiveType = "binary";

	distanceGridSize = 1.0;
	angleGridSize    = 180;
	dihedralGridSize = 60;
	filterVectorPairs = false;
}

VectorHashing::VectorHashing(const VectorHashing &_copyThis){
	archiveType = _copyThis.archiveType;
	distanceGridSize = _copyThis.distanceGridSize;
	angleGridSize = _copyThis.angleGridSize;
	dihedralGridSize = _copyThis.dihedralGridSize;
	filterVectorPairs = _copyThis.filterVectorPairs;
}

VectorHashing::~VectorHashing(){

}


bool VectorHashing::addToVectorHash(System &_sys, string _id, bool _printIt){

	for (uint c =0; c < _sys.chainSize();c++){
		bool result = addToVectorHash(_sys.getChain(c),_id,_printIt);
		if (!result) return false;
		
	}

	return true;
}

bool VectorHashing::addToVectorHash(Chain &_ch,string _id,bool _printIt){

	stringstream chainSpecificId;
	chainSpecificId << _id << ":" ;

	vector<string> chainPositionIds;
	//MSLOUT.stream() << "Working on chain "<<chainSpecificId.str()<<endl;

	// For each position
	for (uint i = 0; i < _ch.positionSize();i++){

		Position &posI = _ch.getPosition(i);
		
		//MSLOUT.stream()<<"\tPosition: "<<posI.toString()<<endl;
		// For each rotamer
		for (uint rotI = 0; rotI < posI.getTotalNumberOfRotamers();rotI++){

			//MSLOUT.stream()<<"\t\tRotamer: "<<rotI<<endl;
			posI.setActiveRotamer(rotI);
			
			// Skip over residues that don't have CA, CB atoms..
			if (!posI.atomExists("CA") || (posI.getResidueName() != "GLY" && !posI.atomExists("CB") ) ) continue;


			Atom *caI  = &posI.getAtom("CA");
			Atom *cbI  = NULL;

			// Compute a CB for glycines
			if (posI.getResidueName() == "GLY"){
				cbI = PDBTopology::getPseudoCbeta(posI.getCurrentIdentity());
			} else {
				cbI = &posI.getAtom("CB");
			}

			stringstream posIid;
			posIid << chainSpecificId.str()<<","<<posI.getRotamerId();
		
			chainPositionIds.push_back(posIid.str());

			for (uint j = i+1; j < _ch.positionSize();j++){


				Position &posJ = _ch.getPosition(j);

				//MSLOUT.stream()<<"\t\t\tPosition: "<<posJ.toString()<<endl;

				
				// For each rotamer
				for (uint rotJ = 0; rotJ < posJ.getTotalNumberOfRotamers();rotJ++){

					//MSLOUT.stream()<<"\t\t\t\tRotamer: "<<rotJ<<endl;

					posJ.setActiveRotamer(rotJ);
			
					// Skip over residues that don't have CA, CB atoms..
					if (!posJ.atomExists("CA") || (posJ.getResidueName() != "GLY" && !posJ.atomExists("CB") ) ) continue;


					Atom *caJ  = &posJ.getAtom("CA");
					Atom *cbJ  = NULL;

					// Compute a CB for glycines
					if (posJ.getResidueName() == "GLY"){
						cbJ = PDBTopology::getPseudoCbeta(posJ.getCurrentIdentity());
					} else {
						cbJ = &posJ.getAtom("CB");
					}


					stringstream posJid;
					posJid << chainSpecificId.str()<<","<<posJ.getRotamerId();

					//MSLOUT.stream() <<"\t\t\t\t\t\t"<<caI->toString()<<endl;
					//MSLOUT.stream() <<"\t\t\t\t\t\t"<<cbI->toString()<<endl;
					//MSLOUT.stream() <<"\t\t\t\t\t\t"<<caJ->toString()<<endl;
					//MSLOUT.stream() <<"\t\t\t\t\t\t"<<cbJ->toString()<<endl;
					VectorPair vp(caI->getCoor(),cbI->getCoor(),caJ->getCoor(),cbJ->getCoor(),posIid.str(),posJid.str());
					vp.calcAll();

					// Filter VectorPair Option
					//					if (filterVectorPairs && filterVectorPair(vp,posI.getResidueName(),posJ.getResidueName())){
					//						continue;
					//					}
					
					pairPositionHash.insert(make_pair(vp.getVectorPairId(),vp));
					//MSLOUT.stream()<< "VP: "<<pairPositionHash[vp.getVectorPairId()].toString()<<endl;

					string geometricKey = getHashKey(vp);
					//if (_printIt) {MSLOUT.stream() << geometricKey<<endl;}
					geometricHash[geometricKey].push_back(&pairPositionHash[vp.getVectorPairId()]);			


					// Visualize the hash
					cout << getHashKey(vp)<<endl;
	

				} // POSJ
			} // J

			positionIds.push_back(chainPositionIds);
		} // POSI

	} // I
	return true;
}


/*
  Discover "complete systems" as defined by:
				  
  Number of total edges = _vh.positionIds.size() - 1;
  Number of "acceptable" edges = _numAcceptableEdges;

  So _vh has position ids A,34   A,56    A,72 , then "complete systems" with 3 edges
				  
  Each VectorPair is complete if:
  A,34 - A,56 matches
  A,34 - A,72 matches
*/
void VectorHashing::searchForVectorMatchAll(VectorHashing &_vh, int _numAcceptableEdges){

	
	for (uint chain1 = 0; chain1 < positionIds.size();chain1++){

		for (uint p1 = 0; p1 < positionIds[chain1].size();p1++){

			//MSLOUT.stream() << "Parse: "<<positionIds[chain1][p1]<<endl;
			vector<string> tokens1 = MslTools::tokenize(positionIds[chain1][p1],":");
			
			string pdb1 = tokens1[0];
			string ch1 = "";
			int res1   = 0;
			string icode1 = "";
			string identity1 = "";
			unsigned int conformation1 = 0;
			if (!MslTools::parseRotamerId(tokens1[1],ch1,res1,icode1,identity1,conformation1)){
				continue;
			}

			map<string,int> p1PositionMatches;
			map<string,bool> p1PairMatches;
			for (uint p2 = p1+1; p2 < positionIds[chain1].size();p2++){

				//MSLOUT.stream() << "Parse2: "<<positionIds[chain1][p2]<<endl;


				vector<string> tokens2 = MslTools::tokenize(positionIds[chain1][p2],":");
			
				string pdb2 = tokens2[0];
				string ch2 = "";
				int res2   = 0;
				string icode2 = "";
				string identity2 = "";
				unsigned int conformation2 = 0;
		
				if (!MslTools::parseRotamerId(tokens2[1],ch2,res2,icode2,identity2,conformation2)){
					continue;
				}

				// Don't look for VectorPairs between pdbs, or different chains, or same position COULD ALSO HAVE RESIDUE SEPARATION SKIP HERE abs(res2-res1) > CUTOFF
				// Or could include ONLY inter-chain but same PDB geometries.
				if (pdb1 != pdb2 || ch1 != ch2  || (res1 == res2 && icode1 == icode2)) continue;

				// Get the VectorPair data
				stringstream key;
				key << pdb1 << ":"<<MslTools::getRotamerId(ch1,res1,icode1,identity1,conformation1)<<"-"<<pdb2 << ":"<<MslTools::getRotamerId(ch2,res2,icode2,identity2,conformation2);

				//MSLOUT.stream() << "KEY: "<<key.str()<<endl;
				map<string,VectorPair>::iterator itPairs;
				itPairs = pairPositionHash.find(key.str());
				if (itPairs == pairPositionHash.end()){
					        MSLOUT.stream() << "KEY "<<key.str()<<" was NOT found in pairPositionHash"<<endl;
						continue;
				}
					
				VectorPair &vp = itPairs->second;

				// Search this vector pair against the input VectorHashing structure
				vector<VectorPair *> *matches = NULL;
				map<string,vector<VectorPair *> >::iterator itGeo;
				itGeo = _vh.geometricHash.find(getHashKey(vp));
				if (itGeo == _vh.geometricHash.end()){
					//MSLOUT.stream() << "Geometry: "<<getHashKey(vp)<<" was NOT found in geometricHash"<<endl;
					continue;
				} else {
					matches = &itGeo->second;
				}

				
				// Keep track of per-position matches
				for (uint m = 0; m < matches->size();m++){

					string matchKey1 = (*matches)[m]->getVectorPairId();
					//MSLOUT.stream() << "\tMATCHKEY1: "<<matchKey1<<endl;
					map<string,bool>::iterator it1;
					it1 = p1PairMatches.find(matchKey1);

					string matchKey2 = (*matches)[m]->getVectorPairIdReverse();
					//MSLOUT.stream() << "\tMATCHKEY2: "<<matchKey2<<endl;

					map<string,bool>::iterator it2;
					it2 = p1PairMatches.find(matchKey2);

					// We have not seen this pair before!
					if (it1 == p1PairMatches.end() || it2 == p1PairMatches.end()){
						p1PositionMatches[matchKey1]++;
						p1PositionMatches[matchKey2]++;

						//MSLOUT.stream() << "UNIQUE FOUND "<<vp.getVectorPairId()<<" and "<<(*matches)[m]->getVectorPairId()<<endl;
					} 
				} // MATCHES END


			} // P2 END


			// Now Analyze P1 for being a complete system
			map<string,int>::iterator it = p1PositionMatches.begin();
			for (;it != p1PositionMatches.end();it++){
				if (it->second >= _numAcceptableEdges){
					fprintf(stdout, "Position %s matches %s\n",positionIds[chain1][p1].c_str(),it->first.c_str());
				}
			}

			
		} // P1 END
	} // CHAIN1 END
}


bool VectorHashing::filterVectorPair(VectorPair &_vp, std::string residueName1, std::string residueName2){

	// TRUE  = REMOVE VECTOR PAIR
	// FALSE = KEEP VECTOR PAIR


	/*
	  Two filtering modes:
	  1) Generic filtering
	     -> Distance > 2 && < 15. Tor -60 to 60 or -160 to -180 or 160 to 180. Angle1,2 ?
	  2) Database filtering
	     -> read in probability function for each residue-residue pair + dist + ang1 + ang2 + tor
	     -> Distance bins 2-15 Angstroms
	     -> Torsion bins 60 degrees?
             -> Angle1,Angle2 ?

	     We may or may not care what types of residues they are.
	 */

	
	return false;
}


string VectorHashing::getHashKey(VectorPair &_vp){
	
	double distanceBin = MslTools::smartRound(_vp.getDistance(),distanceGridSize);

	
	double angle1 = _vp.getAngle1();
	if (angle1 > 180){
		angle1 = 360 - _vp.getAngle1();
	} else if (angle1 < 0){
		angle1 = 180 + _vp.getAngle1();
	}
	double angle2 = _vp.getAngle2();
	if (angle2 > 180){
		angle2 = 360 - _vp.getAngle2();
	} else if (angle2 < 0){
		angle2 = 180 + _vp.getAngle2();
	}

	double angleBin   = MslTools::smartRound(angle1+angle2, angleGridSize);

	double torsion = _vp.getTorsion();
	if (torsion > 180){
		torsion = 360 - _vp.getTorsion();
	} else if (torsion < 0){
		torsion = 180 + _vp.getTorsion();
	}
	double torsionBin  = MslTools::smartRound(torsion,dihedralGridSize);

	
	stringstream key;
	key << distanceBin <<":"<<angleBin<<":"<<torsionBin; 

	return (key.str());
	

}

void VectorHashing::setFilterVectorPairs(bool _flag){
	filterVectorPairs = _flag;
}

bool VectorHashing::getFilterVectorPairs(){
	return filterVectorPairs;
}




/* 
   OLD STUFF 

*/







