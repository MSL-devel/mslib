#include "CharmmSystemBuilder.h"
#include "SasaCalculator.h"
#include "MslTools.h"

#include <fstream>

using namespace MSL;
using namespace std;


void getMissingAtoms(System& sys, char* file, vector<Atom*>& missingAtoms,string pdbId) {
	// read missing atoms file
	map<Atom*,bool> missing;
	ifstream miss(file);
	if(!miss.is_open()) {
		cerr << "ERROR1 Unable to open " << file << endl;
		exit(0);
	}

	char line[10000];
	// Format 4DDP 1,389,ARG,CG 1,389,ARG,CD 1,389,ARG,NE 1,389,ARG,CZ 1,389,ARG,NH1 1,389,ARG,NH2
	while(miss.good()) {
		miss.getline(line,9999);
		string tmp = string(line);
		tmp = MslTools::uncomment(tmp);
		if(tmp.size() > 4 && tmp.substr(0,4) == pdbId)  {
			//cout << tmp << endl;
			vector<string> atoms = MslTools::tokenizeAndTrim(tmp);
			for(int i = 1; i < atoms.size(); i++) {
				string posId;
				vector<string> toks = MslTools::tokenizeAndTrim(atoms[i],",");
				if(sys.positionExists(toks[0] + "," + toks[1])) {
					Position& pos = sys.getLastFoundPosition();
					if(pos.atomExists("CA")) {
						missing[&pos.getLastFoundAtom()] = 1;
					} else {
						cerr << "ERROR2 " << pdbId << " " << pos.getPositionId() << "," << pos.getResidueName() << " CA not found." << endl;
					}
				} else {
					cerr << "ERROR3 " << pdbId << " " << atoms[i] << " not found." << endl;
				}
			}
			break;
		}
	}
	miss.close();

	for(map<Atom*,bool>::iterator it = missing.begin(); it != missing.end(); it++) {
		missingAtoms.push_back(it->first);
		cout << "Missing Atom: " << *(it->first) << endl;
	}
}

void getChainsToUse(char* file, map<string,bool>& chainsToUse,string pdbId) {
	// read the dedup_chainsfile and identify what chain to use
	ifstream dedup(file);
	if(!dedup.is_open()) {
		cerr << "ERROR4 Unable to open " << file << endl;
		exit(0);
	}
	char line[10000];

	while(dedup.good()) {
		dedup.getline(line,9999);
		string tmp = string(line);
		tmp = MslTools::uncomment(tmp);
		// Format : 1A3A A MANLFKLGAENIFLGRKAATKEEAIRFAGEQLVKGGYVEPEYVQAMLDREKLTPTYLGESIAVPHGTVEAKDRVLKTGVVFCQYPEGVRFGEEEDDIARLVIGIAARNNEHIQVITSLTNALDDESVIERLAHTTSVDEVLELLAGRK
		if(tmp.size() > 4 && tmp.substr(0,4) == pdbId) {
			chainsToUse[tmp.substr(5,1)] = true;
			cout << "Use oldChainId: " << tmp.substr(5,1) << endl;
		}

	}
	dedup.close();
}

void mapPositionIds(char* file, map<string,string>& oldPositionIds,map<string,string>& newPositionIds,string pdbId) {
	// read the reachinedResidues.txt file
	ifstream rechain(file);
	if(!rechain.is_open()) {
		cerr << "ERROR5 Unable to open " << file << endl;
		exit(0);
	}
	char line[10000];

	while(rechain.good()) {
		rechain.getline(line,9999);
		string tmp = string(line);
		tmp = MslTools::uncomment(tmp);
		// Format : 3QT5 B,327 1,327
		if(tmp.size() > 4 && tmp.substr(0,4) == pdbId) {
			vector<string> tokens = MslTools::tokenize(tmp);	
			oldPositionIds[tokens[2]] = tokens[1];
			newPositionIds[tokens[1]] = tokens[2];
		}

	}
	rechain.close();

}

void updateChainsToUse(map<string,bool>& chainsToUse,map<string,string>& oldPositionIds) {
	for(map<string,string>::iterator it = oldPositionIds.begin(); it != oldPositionIds.end(); it++) {
		string newChain = (MslTools::tokenize(it->first,","))[0];
		string oldChain = (MslTools::tokenize(it->second,","))[0];
		if(chainsToUse.find(oldChain) != chainsToUse.end()) {
			if(chainsToUse.find(newChain) == chainsToUse.end()) {
				chainsToUse[newChain] = 1;
				cout << "Use newChainId: " << newChain << endl;
			}
		}
	}
}

void getResidueBFactor(char* file,map<string,double>& residueBFactor, map<string,string>& newPositionIds) {
	// read the pdbId.bfac file
	ifstream bfac(file);
	if(!bfac.is_open()) {
		cerr << "ERROR6 Unable to open " << file << endl;
		exit(0);
	}

	char line[10000];
	while(bfac.good()) {
		bfac.getline(line,9999);
		string tmp = string(line);
		tmp = MslTools::uncomment(tmp);
		// Format : A,10,HIS 22.81
		if(tmp.size() > 4 ) {
			vector<string> tokens = MslTools::tokenize(tmp);	
			vector<string> tmp = MslTools::tokenize(tokens[0],",");
			string oldPosId = tmp[0] + "," + tmp[1]; 
			if(newPositionIds.find(oldPosId) != newPositionIds.end()) {
				residueBFactor[newPositionIds[oldPosId] + "," + tmp[2]] = MslTools::toDouble(tokens[1]);
			} else {
				residueBFactor[tokens[0]] = MslTools::toDouble(tokens[1]);
			}
		}

	}
	bfac.close();


}
void getDisulfides(char* file, map<string,bool>& disulfides,string pdbId) {
	// read the disulfides_reachined.txt file
	ifstream disu(file);
	if(!disu.is_open()) {
		cerr << "ERROR7 Unable to open " << file << endl;
		exit(0);
	}
	char line[10000];

	while(disu.good()) {
		disu.getline(line,9999);
		string tmp = string(line);
		tmp = MslTools::uncomment(tmp);
		// Format : 2Y32 D 39 D 69
		if(tmp.size() > 4 && tmp.substr(0,4) == pdbId) {
			vector<string> tokens = MslTools::tokenize(tmp);	
			if(tokens.size() < 5) {
				disulfides[tokens[1] + "," + tokens[2]] = 1;
				disulfides[tokens[3] + "," + tokens[4]] = 1;
			} else {
				cerr << "ERROR11 Wrong format <" << tmp << "> in disulfides file" << endl;
			}
		}
	}
	disu.close();

}


int main(int argc, char* argv[]) {
	if(argc != 9) {
		cerr << "Usage: selectCavities <pdbid> <missingAtomsFile> <dedup_chainsfile> <missing_atoms_dist_threshold> <bfactorfile> <rechainedResidues.txt> <disulfideIds> <outfile>" << endl;
		cerr << "	Expects files of the form pdbid.pdb in the current directory" << endl;
		cerr << "	Prints all residues that are atleast <missing_atoms_dist_threshold> far (CA-CA distance) from any residue with a missing atom in the original pdb." << endl;
		cerr << "	Prints only residues from unique chains in the pdb." << endl;
		cerr << "	Skips terminal residues in each chain" << endl;
		cerr << "	Excludes residues with bfactor >= 40" << endl;
		cerr << "	Excludes CYS with disulfide bonds" << endl;
		exit(0);
	}


	double distThresh = MslTools::toDouble(string(argv[4]));
	string pdbId = string(argv[1]);
	System sys;
	if(!sys.readPdb(pdbId + ".pdb")) {
		cerr << "ERROR8 Unable to build " << argv[1] << endl;
		exit(0);
	}

	vector<Atom*> missingAtoms; // really the CAs of residues with missing atoms
	getMissingAtoms(sys,argv[2],missingAtoms,pdbId);

	map<string,string> oldPositionIds; // map from new positionId to old positionId
	map<string,string> newPositionIds; // map from old positionId to new positionId
	mapPositionIds(argv[6],oldPositionIds,newPositionIds,pdbId);

	map<string,bool> chainsToUse;
	getChainsToUse(argv[3],chainsToUse,pdbId);

	// the original chain may be split up into multiple numbered chains so identify and consider residues from these chains too
	// get a list of residueids that are ok to use
	// convert old chainIds to the 0-based new chainIds
	// if chain A is good to use and was split into A,0,1 then A,0,1 should be usable
	updateChainsToUse(chainsToUse,oldPositionIds);

	map<string,double> residueBFactor;
	getResidueBFactor(argv[5],residueBFactor,newPositionIds);

	map<string,bool> disulfidePosId;
	getDisulfides(argv[7],disulfidePosId,pdbId);

	ofstream fout;
	fout.open(argv[8]);

	if(!fout.is_open()) {
		cerr << "ERROR9 Unable to open " << argv[8] << endl;
		exit(0);
	}

	

	for(int i = 0; i < sys.positionSize(); i++) {
		Residue& res = sys.getIdentity(i);
		string charmmResName = res.getResidueName();
		string pdbResName = res.getResidueName();
		if(charmmResName == "HSD" || charmmResName == "HSE" || charmmResName == "HSP" ) {
			pdbResName = "HIS";
		}

		string pdbIdentityId = res.getPositionId() + "," + pdbResName;
		string charmmIdentityId = res.getPositionId() + "," + charmmResName;

		// is bfactor <= 40?
		if (residueBFactor.find(pdbIdentityId) == residueBFactor.end()) {
			cerr << "ERROR10 " << pdbIdentityId << " b-factor not found" << endl;
			continue;
		}

		if(residueBFactor[pdbIdentityId]  >= 40) {
			cout << "REASON1 " << charmmIdentityId << " b-factor >= 40" << endl;
			continue;
		}

			
		Position* pPos = res.getParentPosition();
		Chain* pChain = pPos->getParentChain();

		// is it a duplicate chain?
		if(chainsToUse.find(res.getChainId()) == chainsToUse.end() ) {
			cout << "REASON2 " << charmmIdentityId << " on a duplicate chain " << endl;
			continue;
		}
		// is this the first or last position?
		int posIndex = pPos->getIndexInChain();
		int positionSize = pChain->positionSize(); 
		if(posIndex == 0 || posIndex == positionSize - 1) {
			cerr << "REASON3 " << charmmIdentityId << " has posIndex " << posIndex << " in chain of size " << positionSize << endl;
			continue;
		}

		// is it a disulfide? 
		if(disulfidePosId.find(pPos->getPositionId()) != disulfidePosId.end()) {
			cerr << "REASON4 " << charmmIdentityId << " is a disulfide " << endl;
			continue;
		}

		// is the CA of this residue sufficiently far from CA of residues with missing atoms?
		Atom* CA = NULL;
		if(!res.atomExists("CA")) {
			cerr <<  "REASON5 " << charmmIdentityId << " CA does not exist" << endl;
			continue;
		}
		CA = &res.getLastFoundAtom();
		int far = true;
		for(int j = 0; j < missingAtoms.size(); j++) {
			if(CA->distance(*missingAtoms[j]) < distThresh) {
				far = false;
				cerr <<  "REASON6 " << charmmIdentityId << " distance from " << missingAtoms[j]->getAtomId() << " less than " << distThresh << endl;
				break;
			}
		}
		if(far) {
			fout << pdbId << " " << charmmIdentityId << endl;
		}

	}
	cout << "Done " << pdbId << endl;
	fout.close();

}





