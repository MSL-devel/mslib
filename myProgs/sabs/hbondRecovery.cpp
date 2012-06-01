#include "HydrogenBondBuilder.h"
#include "CharmmSystemBuilder.h"
#include "SasaCalculator.h"
#include "System.h"

using namespace std;
using namespace MSL;

map<string,double> refSasa;
string pdbName ;
string topfile;
string parfile;
string hbondfile;

/* to check for symmetric bonds

ARG : HH21 = HH22 = HH11 = HH12
ASN : HD21 = HD22
ASP : OD1 = OD2
GLN : HE21 = HE22
GLU : OE1 = OE2
LYS : HZ1 = HZ2 = HZ3

*/

string getCorrectedAtomId(Atom* _a ) {
	string resName = _a->getResidueName();
	string atomName = _a->getName();
	if(resName == "ARG" && atomName.substr(0,2) == "HH") {
		return atomName.substr(0,2);
	}
	if(resName == "ASN" && atomName.substr(0,3) == "HD2") {
		return atomName.substr(0,3);
	}
	if(resName == "ASP" && atomName.substr(0,2) == "OD") {
		return atomName.substr(0,2);
	}
	if(resName == "GLN" && atomName.substr(0,3) == "HE2") {
		return atomName.substr(0,3);
	}
	if(resName == "GLU" && atomName.substr(0,2) == "OE") {
		return atomName.substr(0,2);
	}
	if(resName == "LYS" && atomName.substr(0,2) == "HZ") {
		return atomName.substr(0,2);
	}
	return atomName;
}
      

bool isBackboneAtom(string name) {
	if(name == "N" || name == "HN" || name == "C" || name == "CA" || name == "O") {
		return true;
	}
	return false;
}

void collectHB(char* pdb, map<string,double>& bonds,map<string,double>& buried ) {
	System sys;
	CharmmSystemBuilder csb(sys,topfile,parfile);
	csb.setBuildNonBondedInteractions(false);
	if(!csb.buildSystemFromPDB(string(pdb))) {
		cerr << "Unable to build " << pdb << endl;
		exit(0);
	}
	sys.buildAllAtoms();
	HydrogenBondBuilder hb(sys,hbondfile);
	hb.buildInteractions(10.0);

	SasaCalculator sc(sys.getAtomPointers());
	sc.calcSasa();

	sys.calcEnergy() ;
	//cout << sys.calcEnergy() << endl;
	//sys.printEnergySummary();

	EnergySet *eSet = sys.getEnergySet();
	//cout << eSet->isTermActive("SCWRL4_HBOND") << endl;

	vector<Interaction*> interactions = (*(eSet->getEnergyTerms()))["SCWRL4_HBOND"];
	//cout << pdb << " " << interactions.size() << endl;
	vector<Atom*> atoms;

	for(int i = 0; i < interactions.size(); i++) {
		double e = interactions[i]->getEnergy();
		//cout << e << endl;
		if(e < 0) {
			//cout << e << endl;
			atoms = interactions[i]->getAtomPointers();
			//cout << atoms[0]->getIdentityId() << "," << atoms[0]->getName() << " " << atoms[2]->getIdentityId() << "," << atoms[2]->getName() << " " << e << endl;
			if(!isBackboneAtom(atoms[0]->getName()) || !isBackboneAtom(atoms[2]->getName())) {
				bonds[atoms[0]->getIdentityId() + "," + getCorrectedAtomId(atoms[0]) + "," + atoms[2]->getIdentityId() + "," + getCorrectedAtomId(atoms[2])] = e;
				Residue* res1 = atoms[0]->getParentResidue();
				Residue* res2 = atoms[2]->getParentResidue();
				if(res1->getSasa() < 0.25 * refSasa[res1->getResidueName()] || res2->getSasa() < 0.25 * refSasa[res2->getResidueName()]) {
					buried[atoms[0]->getIdentityId() + "," + getCorrectedAtomId(atoms[0]) + "," + atoms[2]->getIdentityId() + "," + getCorrectedAtomId(atoms[2])] = e;
				}
				//bonds[atoms[0]->getAtomId() + " " + atoms[2]->getAtomId()] = e;
				//cout << atoms[0]->getIdentityId() << "," << atoms[0]->getName() << " " << atoms[2]->getIdentityId() << "," << atoms[2]->getName() << " " << e << endl;
			}
		}
	}

}

void printStats(string prefix,map<string,double>& origHB, map<string,double>& recoveredHB) {
	map<string,int> orig;
	map<string,int> repacked;
	map<string,int> recovered;

	int totRecovered = 0;

	for(map<string,double>::iterator it = recoveredHB.begin(); it != recoveredHB.end(); it++) {
		vector<string> toks = MslTools::tokenize(it->first,",");
		// toks[2] - acc toks[6] - donor
		// toks[3] - acc Atom toks[7] - donAtom
		if(!isBackboneAtom(toks[3])) {
			if(repacked.find(toks[2]) != repacked.end()) {
				repacked[toks[2]]++;
			} else {
				repacked[toks[2]] = 1;
			}
		}
		if(!isBackboneAtom(toks[7])) {
			if(repacked.find(toks[6]) != repacked.end()) {
				repacked[toks[6]]++;
			} else {
				repacked[toks[6]] = 1;
			}
		}
	}

	for(map<string,double>::iterator it = origHB.begin(); it != origHB.end(); it++) {
		vector<string> toks = MslTools::tokenize(it->first,",");
		if(!isBackboneAtom(toks[3])) {
			if(orig.find(toks[2]) != orig.end()) {
				orig[toks[2]]++;
			} else {
				orig[toks[2]] = 1;
			}
		}
		if(!isBackboneAtom(toks[7])) {
			if(orig.find(toks[6]) != orig.end()) {
				orig[toks[6]]++;
			} else {
				orig[toks[6]] = 1;
			}
		}
		if(recoveredHB.find(it->first) != recoveredHB.end() ) {
			totRecovered++;
			if(!isBackboneAtom(toks[3])) {
				if(recovered.find(toks[2]) != recovered.end()) {
					recovered[toks[2]]++;
				} else {
					recovered[toks[2]] = 1;
				}
			}
			if(!isBackboneAtom(toks[7])) {
				if(recovered.find(toks[6]) != recovered.end()) {
					recovered[toks[6]]++;
				} else {
					recovered[toks[6]] = 1;
				}
			}
		}
	}

	for(map<string,int>::iterator it = orig.begin(); it != orig.end(); it++) {
		if(it->second == 0) {
			continue;
		}
		cout << prefix << "Recovery: " << pdbName << " " << it->first << " " << it->second << " "  ;
		if(recovered.find(it->first) != recovered.end()) {
			if(it->second != 0) {
				cout << recovered[it->first] << endl;
			} else {
				cout << "0" << endl;
			}
		} else {
			cout << "0" << endl;
		}
	}
	
	cout << prefix << "Total: " << pdbName << " " << origHB.size() << " " << totRecovered << " " << recoveredHB.size() << endl;
}


int main (int argc, char* argv[]) {

	if(argc != 7) {
		cerr << "Usage: hbondRecovery <origPDB> <repackedPDB> <pdbName> <topfile> <parfile> <hbondfile>" << endl;
		cerr << "Prints out the number of hbonds recovered" << endl;
		exit(0);
	}

	topfile = string(argv[4]);
	parfile = string(argv[5]);
	hbondfile = string(argv[6]);

	refSasa["ALA"] = 116.40;
	refSasa["ARG"] = 249.26;
	refSasa["ASN"] = 168.87;
	refSasa["ASP"] = 155.37;
	refSasa["CYS"] = 141.48;
	refSasa["GLU"] = 187.16;
	refSasa["GLN"] = 189.17;
	refSasa["GLY"] = 83.91;
	refSasa["HSD"] = 198.51;
	refSasa["HSP"] = 198.51;
	refSasa["HSE"] = 198.51;
	refSasa["ILE"] = 189.95;
	refSasa["LEU"] = 197.99;
	refSasa["LYS"] = 207.49;
	refSasa["MET"] = 210.55;
	refSasa["PHE"] = 223.29;
	refSasa["PRO"] = 144.80;
	refSasa["SER"] = 125.68;
	refSasa["THR"] = 148.06;
	refSasa["TRP"] = 265.42;
	refSasa["TYR"] = 238.30;
	refSasa["VAL"] = 162.24;


	cout << topfile << " " << parfile << " " << hbondfile << endl;
	pdbName = string(argv[3]);

	map<string,double> origHB;
	map<string,double> recoveredHB;

	map<string,double> origBuriedHB;
	map<string,double> recoveredBuriedHB;


	collectHB(argv[1],origHB,origBuriedHB);
	//cout << "**************************************************************************************************************" << endl;
	collectHB(argv[2],recoveredHB,recoveredBuriedHB);

	printStats("",origHB,recoveredHB);
	printStats("Buried",origBuriedHB,recoveredHB);

	/* old */
/*
	int tot = 0;
	int recovered = 0;
	
	for(map<string,double>::iterator it = origHB.begin(); it != origHB.end(); it++) {
		tot++;
		if(recoveredHB.find(it->first) != recoveredHB.end() ) {
			cout << it->first << endl;
			recovered++;
		}
	}

	cout << argv[3] << " " << tot << " " << recovered << endl;
	exit(0);


	*/

	

}
