/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A \"Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)\" (2012) J. Comput. Chem, 33, 1645-61 
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

#include "VectorHashing.h"


#include "MslOut.h"
static MslOut MSLOUT("VectorHashing");

VectorHashing::VectorHashing(){
	archiveType = "binary";

	distanceGridSize = 0.5;
	angleGridSize    = 75;
	dihedralGridSize = 10;
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
	MSLOUT.stream() << "Working on chain "<<chainSpecificId.str()<<endl;

	// For each position
	for (uint i = 0; i < _ch.positionSize();i++){

		Position &posI = _ch.getPosition(i);
		
		MSLOUT.stream()<<"\tPosition: "<<posI.toString()<<endl;

		// For each rotamer
		for (uint rotI = 0; rotI < posI.getTotalNumberOfRotamers();rotI++){

			MSLOUT.stream()<<"\t\tRotamer: "<<rotI<<endl;
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
			MSLOUT.stream() << "Adding "<<posIid.str()<< " to chainPositionIds"<<endl;
			for (uint j = i+1; j < _ch.positionSize();j++){


				Position &posJ = _ch.getPosition(j);

				MSLOUT.stream()<<"\t\t\tPosition: "<<posJ.toString()<<endl;

				
				// For each rotamer
				for (uint rotJ = 0; rotJ < posJ.getTotalNumberOfRotamers();rotJ++){

					MSLOUT.stream()<<"\t\t\t\tRotamer: "<<rotJ<<endl;

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
					//cout << getHashKey(vp)<<endl;
	

				} // POSJ
			} // J

		} // POSI

	} // I

	positionIds.push_back(chainPositionIds);
	return true;
}

bool VectorHashing::getPositionId(VectorPair &_vp, string &_posId1, string &_posId2){
	string id = _vp.getVectorPairId();
	vector<string> tokens1 = MslTools::tokenize(id,"-");
	
	if (tokens1.size() < 2){
		cerr << "ERROR VectorHashing::getPositionId(VectorPair) doesn't have '-' character"<<endl;
		return false;
	}
	
	vector<string> tokens2 = MslTools::tokenize(tokens1[0],":");
	if (tokens2.size() < 2){
		cerr << "ERROR VectorHashing::getPositionId(VectorPair) doesn't have ':' character"<<endl;
		return false;
	}	
	string rotId1 = tokens2[1];
	
	if (tokens2[1].substr(0,1) == ",") {
		rotId1 = tokens2[1].substr(1,tokens2[1].size());
	}

	string ch1 = "";
	int res1   = 0;
	string icode1 = "";
	string identity1 = "";
	unsigned int conformation1 = 0;
	if (!MslTools::parseRotamerId(rotId1,ch1,res1,icode1,identity1,conformation1)){
		return false;
	}

	_posId1 = MslTools::getPositionId(ch1,res1,icode1);


	tokens2 = MslTools::tokenize(tokens1[1],":");
	if (tokens2.size() < 2){
		cerr << "ERROR VectorHashing::getPositionId(VectorPair) doesn't have ':' character2"<<endl;
		return false;
	}	
	string rotId2 = tokens2[1];
	
	if (tokens2[1].substr(0,1) == ",") {
		rotId2 = tokens2[1].substr(1,tokens2[1].size());
	}

	string ch2 = "";
	int res2   = 0;
	string icode2 = "";
	string identity2 = "";
	unsigned int conformation2 = 0;
	if (!MslTools::parseRotamerId(rotId2,ch2,res2,icode2,identity2,conformation2)){
		return false;
	}

	_posId2 = MslTools::getPositionId(ch2,res2,icode2);


	return true;
}
/*
  Discover "complete systems" as defined by:
				  
  Number of total edges = _vh.positionIds.size() - 1;
  Number of "acceptable" edges = _numAcceptableEdges;

  So _vh has position ids A,34   A,56    A,72 , then "complete systems" with 2 edges per position
				  
  Each VectorPair is complete if:
  A,34 - A,56 matches
  A,34 - A,72 matches
*/
vector<map<string,vector<string> > > VectorHashing::searchForVectorMatchAll(VectorHashing &_vh, int _numAcceptableEdges){

        vector<map<string,vector<string> > > cycles;
	map<string,map<string,map<string,vector<string> > > > edgeTrees;

	for (uint chain1 = 0; chain1 < positionIds.size();chain1++){
		
		for (uint p1 = 0; p1 < positionIds[chain1].size();p1++){

			//MSLOUT.stream() << "Parse: "<<positionIds[chain1][p1]<<" chain "<<chain1<<" pos "<<p1<<endl;
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

			// matched.position1,matched.position2 -> [ edge1, edge2, ...]  where edges are complete inverse rotamer descriptions
			map<string,map<string,map<string,vector<pair<string,string> > > > > p1EdgeTree;
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

				//MSLOUT.stream() << "Geometry: "<<getHashKey(vp)<<" was found in geometricHash"<<endl;				
				// Keep track of per-position matches
				for (uint m = 0; m < matches->size();m++){

					string matchedPos1 = "";
					string matchedPos2 = "";
					getPositionId(*(*matches)[m],matchedPos1,matchedPos2);

					string matchedEdge1 = matchedPos1 + "-" + matchedPos2;
					string matchedEdge2 = matchedPos2 + "-" + matchedPos1;

					p1EdgeTree[positionIds[chain1][p2]][matchedPos1][matchedPos2].push_back(pair<string,string>((*matches)[m]->getVectorAid(),(*matches)[m]->getVectorBid()));
					//p1EdgeTree[matchedPos2][matchedPos1].push_back((*matches)[m]->getVectorPairIdReverse());

					//MSLOUT.stream() << "UNIQUE FOUND "<<vp.getVectorPairId()<<" and "<<matchedEdge1<<" "<<p1PositionMatches[matchedEdge1]<<" "<<_numAcceptableEdges<<endl;
				} // MATCHES END


			} // P2 END


			// Now Analyze P1 for being a complete system

			// Does P1 have enough matched geometries to be considered a cycle?
			map<string, vector<string> > cycle;
			if (p1EdgeTree.size() >= _numAcceptableEdges){

				// For each P2
				map<string,map<string,map<string,vector<pair<string,string> > > > >::iterator it = p1EdgeTree.begin();
				for (;it != p1EdgeTree.end();it++){


					// For each matched position1
					map<string,map<string,vector<pair<string,string> > > >::iterator it2 = it->second.begin();
					for (;it2 != it->second.end();it2++){



						// For each matched position 2
						map<string,vector<pair<string, string> > >::iterator it3 = it2->second.begin();
						for (;it3 != it2->second.end();it3++){					

						  //fprintf(stdout, "%s = %s and %s = %s  OR  %s = %s and %s = %s with %6d rotamers\n",
						  //	positionIds[chain1][p1].c_str(),it2->first.c_str(),it->first.c_str(),it3->first.c_str(),
						  //	positionIds[chain1][p1].c_str(),it3->first.c_str(),it->first.c_str(),it2->first.c_str(), (int)it3->second.size());

							
							stringstream ss;
							ss << positionIds[chain1][p1] <<"--"<<it2->first;

							stringstream ss2;
							ss2 << it->first << "--"<<it3->first;
							cycle[ss.str()].push_back(ss2.str());


							ss.str("");
							ss << positionIds[chain1][p1] <<"--"<<it3->first;
							ss2.str("");
							ss2 << it->first << "--"<<it2->first;
							cycle[ss.str()].push_back(ss2.str());							


						}
					}
				}
			}

			map<string,vector<string> > completeCycles;
			map<string, vector<string> >::iterator cycleIt = cycle.begin();
			for (;cycleIt != cycle.end();cycleIt++){

			  if (cycleIt->second.size() == 2){
			    completeCycles[cycleIt->first] = cycleIt->second;

			    stringstream ss;
			    ss << "CYCLE: "<<cycleIt->first;
			    for (uint c = 0; c < cycleIt->second.size();c++){
			      ss << " AND "<< cycleIt->second[c];
			    }
			    //cout << ss.str()<<endl;
			  }
			}
				//edgeTrees[positionIds[chain1][p1]] = p1EdgeTree;

			cycles.push_back(completeCycles);
		} // P1 END
	} // CHAIN1 END

	return cycles;
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
	
	double distanceBin = MslTools::smartRound(_vp.getDistance1(),distanceGridSize);

	
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

	double torsion = _vp.getTorsion1();
	if (torsion > 180){
		torsion = 360 - _vp.getTorsion1();
	} else if (torsion < 0){
		torsion = 180 + _vp.getTorsion1();
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







