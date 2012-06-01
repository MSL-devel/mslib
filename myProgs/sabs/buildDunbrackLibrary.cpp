#include "MslTools.h"
#include "System.h"
#include "AtomPointerVector.h"
#include "PolymerSequence.h"
#include "Position.h"
#include "CharmmSystemBuilder.h"
#include "AtomSelection.h"
#include "RotamerLibrary.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"

#include <vector>

using namespace std;

using namespace MSL;

RotamerLibrary rotlib;

map<string,int> chi1pos;
map<string,int> chi2pos;
map<string,int> chi3pos;
map<string,int> chi4pos;

void addConformation(string library, string resName, vector<double>& values) {
	if(resName == "LYS") {
	//	values[14] = 30;
	//	rotlib.addConformation(library,resName,values);
	//	values[14] = 60;
	//	rotlib.addConformation(library,resName,values);
	//	values[14] = 90;
		rotlib.addConformation(library,resName,values);
	} else if(resName == "SER" || resName == "THR") {
		values[5] = 30;
		rotlib.addConformation(library,resName,values);
		values[5] = 60;
		rotlib.addConformation(library,resName,values);
		values[5] = 90;
		rotlib.addConformation(library,resName,values);
		values[5] = 150;
		rotlib.addConformation(library,resName,values);
		values[5] = 180;
		rotlib.addConformation(library,resName,values);
		values[5] = -150;
		rotlib.addConformation(library,resName,values);
		values[5] = -90;
		rotlib.addConformation(library,resName,values);
		values[5] = -60;
		rotlib.addConformation(library,resName,values);
		values[5] = -30;
		rotlib.addConformation(library,resName,values);
	} else if(resName == "TYR") {
		values[15] = 0;
		rotlib.addConformation(library,resName,values);
		values[15] = 45;
		rotlib.addConformation(library,resName,values);
		values[15] = 90;
		rotlib.addConformation(library,resName,values);
		values[15] = 135;
		rotlib.addConformation(library,resName,values);
		values[15] = 180;
		rotlib.addConformation(library,resName,values);
		values[15] = -135;
		rotlib.addConformation(library,resName,values);
		values[15] = -90;
		rotlib.addConformation(library,resName,values);
		values[15] = -45;
		rotlib.addConformation(library,resName,values);
		
	} else {
		rotlib.addConformation(library,resName,values);
	}

}

	
void buildLibrary(string resName, vector<vector<double> >& _chis) {	
	
	int resNum = 2; // position of residue in the sequence built (GLY-X-GLY)
	if(resName == "PRO") {
		return;
	}
	PolymerSequence seq("A: GLY " + resName + " GLY ");
	
	System sys;
	CharmmSystemBuilder csb(sys,"/data00/sabs/pdb/top_all22_prot.inp", "/data00/sabs/pdb/par_all22_prot.inp");

	if(!csb.buildSystem(seq)) {
		cerr << "Unable to build " << endl;
		exit(0);
	}
	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	
	sys.buildAtoms();


	if(!sys.positionExists("A",resNum,"")) {
		cerr <<  "Residue not found " << endl;
		exit(0);
	}	

	Position p = sys.getLastFoundPosition();

	int num_chis = _chis[0].size();
	//cout << resName << " NumChis: " << num_chis << " NumRots: " << num_rots << endl;
	
	int no_rots = _chis.size();
	
	vector<string> initAtomNames(rotlib.getMobileAtoms("",resName));
	vector<Atom*> initAtoms;
//	cout << "init atoms: " << endl;
	for(int j = 0; j < initAtomNames.size(); j++) {
		if(!sys.atomExists("A", resNum,"",initAtomNames[j])) {
			cerr << initAtomNames[j] << " not found " << endl;
			exit(0);
		}
		Atom* atom = &(sys.getLastFoundAtom());
	//	cout << atom->getName() << " ";
		initAtoms.push_back(atom);
	}
//	cout << endl << "Defi " << endl;

	vector<string> defi(rotlib.getInternalCoorDefinitionLines("",resName));
	vector <double> values;

	for(int i = 0; i < defi.size(); i++) {
		vector<string> atom_names = MslTools::tokenizeAndTrim(defi[i]," ",false,"\t\n\r *");
		vector<Atom*> atoms;
		double val = -10000 ;
		for(int j = 1; j < atom_names.size(); j++) {
			if(!sys.atomExists("A" , resNum,"" , atom_names[j])) {
				cerr << atom_names[j] << " not found " << endl;
				exit(0);
			}
			Atom* atom = &(sys.getLastFoundAtom());
		//	cout << atom->getName() << " ";
		//	cout << *atom << endl;
			atoms.push_back(atom);
		}

		if(atoms.size() == 2)  {
			//bond
			val = atoms[0]->distance(*(atoms[1]));

		} else if (atoms.size() == 3)  {
			//angle
			val = atoms[0]->angle(*(atoms[1]),*(atoms[2]));
		} else if (atoms.size() == 4) {
			//dihedral
			val = atoms[0]->dihedral(*(atoms[1]),*(atoms[2]),*(atoms[3]));
		} else {
			cerr << "Error" << endl;
			exit(0);
		}


	//	cout << " " << val << endl;
		values.push_back(val);
		
	}


	
	for(int i = 0; i < no_rots; i++) {
		values[chi1pos[resName]]= _chis[i][0];
		if(num_chis > 1) {
			if(resName != "TRP") {
				values[chi2pos[resName]] = _chis[i][1];
			} else {
				if(_chis[i][1] > 0) {
					values[chi2pos[resName]] = _chis[i][1] - 180;
				} else {
					values[chi2pos[resName]] = _chis[i][1] + 180;
				}
			}

		}
		if(num_chis > 2) {
			values[chi3pos[resName]] = _chis[i][2];
		}
		if(num_chis > 3) {
			values[chi4pos[resName]] = _chis[i][3];
		}
		addConformation("",resName,values);
	}

}

void initializeChiPos() {
	// in our rotamer library these are the corresponding coordinates
	chi1pos["ARG"] = 2;
	chi1pos["ASN"] = 2;
	chi1pos["ASP"] = 2;
	chi1pos["CYS"] = 2;
	chi1pos["GLN"] = 2;
	chi1pos["GLU"] = 2;
	chi1pos["HSD"] = 2;
	chi1pos["HSE"] = 2;
	chi1pos["HSP"] = 2;
	chi1pos["ILE"] = 2;
	chi1pos["LEU"] = 2;
	chi1pos["LYS"] = 2;
	chi1pos["MET"] = 2;
	chi1pos["PHE"] = 2;
	chi1pos["SER"] = 2;
	chi1pos["THR"] = 2;
	chi1pos["TRP"] = 2;
	chi1pos["TYR"] = 2;
	chi1pos["VAL"] = 2;

	chi2pos["ARG"] = 5;
	chi2pos["ASN"] = 5;
	chi2pos["ASP"] = 5;
	chi2pos["GLN"] = 5;
	chi2pos["GLU"] = 5;
	chi2pos["HSD"] = 5;
	chi2pos["HSE"] = 5;
	chi2pos["HSP"] = 5;
	chi2pos["ILE"] = 8; // ILE is different in our library
	chi2pos["LEU"] = 5;
	chi2pos["LYS"] = 5;
	chi2pos["MET"] = 5;
	chi2pos["PHE"] = 5;
	chi2pos["TRP"] = 5;
	chi2pos["TYR"] = 5;

	chi3pos["ARG"] = 8;
	chi3pos["GLN"] = 8;
	chi3pos["GLU"] = 8;
	chi3pos["LYS"] = 8;
	chi3pos["MET"] = 8;

	chi4pos["ARG"] = 11;
	chi4pos["LYS"] = 11;
}

int main(int argc, char* argv[]) {

	if(argc != 3) {
		cerr << "Usage: buildDunbrackLibrary chisFile outFile\n";
		cerr << "Format of chisFile: " << endl;
		cerr << "	ARG 109.6 105.8 60.0 61.9" << endl;
		cerr << "	LEU 60.0 61.9" << endl;
		exit(0);
	}

	initializeChiPos();

	map<string,vector<vector<double> > > chis;
	ifstream chiFile;
	chiFile.open(argv[1]);
	if(!chiFile.is_open()) {
		cerr << "Unable to open " << argv[1] << endl;
		exit(0);
	}

	char temp[1000];
	while(chiFile.getline(temp,999)) {
		vector<string> toks = MslTools::tokenizeAndTrim(string(temp));
		if(toks.size() < 2) {
			continue;
		}
		chis[toks[0]].push_back(vector<double>());
		for(int i = 1; i < toks.size(); i++) {
			chis[toks[0]].back().push_back(MslTools::toDouble(toks[i]));
		}
	}
/*
	for(map<string,vector<vector<double> > >::iterator it = chis.begin(); it != chis.end(); it++) {
		for(int i = 0; i < it->second.size(); i++) {
			cout << it->first;
			for(int j = 0; j < (it->second)[i].size(); j++) {
				cout << " " << (it->second)[i][j] ;
			}
			cout << endl;
		}
	}
*/
	rotlib.readFile("/library/rotlib/balanced/rotlib-balanced-200.txt");
	rotlib.removeAllConformations();

	for(map<string,vector<vector<double> > >::iterator it = chis.begin(); it != chis.end(); it++) {
		if(it->first == "HIS") {
			buildLibrary("HSD",it->second);
			buildLibrary("HSE",it->second);
			buildLibrary("HSP",it->second);
		} else {
			buildLibrary(it->first,it->second);
		}
	}

	if(!rotlib.writeFile(string(argv[2]))) {
		cerr << "Unable to write " << argv[2] << endl;
	}
	return 0;
}

