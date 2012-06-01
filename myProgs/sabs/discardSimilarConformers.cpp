#include "System.h"
#include "AtomPointerVector.h"
#include "PolymerSequence.h"
#include "Position.h"
#include "Residue.h"
#include "PDBReader.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "RotamerLibrary.h"
#include "RotamerLibraryBuilder.h"
#include "RotamerLibraryReader.h"
#include "RotamerLibraryWriter.h"
#include "AtomSelection.h"
#include <math.h>
#include <stdlib.h>

using namespace std;
using namespace MSL;


void discardSimilarConformers (double threshold,string residue, string inpFile, string outFile ) {
	PolymerSequence ps("A: " + residue + "\nB: " + residue + "\n");
	System sys;
	CharmmSystemBuilder csb(sys,"/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp");
	if(!csb.buildSystem(ps)) {
		cerr << "Unable to build CharmmSystem" << endl;
		exit(0);
	}
	SystemRotamerLoader rotLoader;
	rotLoader.setSystem(sys);
	
	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	sys.buildAllAtoms();

	if(!sys.positionExists("A,1")) {
		cerr << "unable to find residue " << residue << " in chain A" << endl;
		exit(0);
	}
	Position * pPosA1 = &(sys.getLastFoundPosition()); 


	if(!sys.positionExists("B,1")) {
		cerr << "unable to find residue " << residue << " in chain B" << endl;
		exit(0);
	}


	Position * pPosB1 = &(sys.getLastFoundPosition()); 

	RotamerLibrary pRotlib_base; // Will eventually have all unique conformers
	pRotlib_base.readFile(inpFile); // some sample rotamerlibrary doesn't matter which one except it should use the same library name
	pRotlib_base.removeAllConformations();

	
	RotamerLibrary pRotlib_All(pRotlib_base); 
	if(!pRotlib_All.readFile(inpFile)) { // Use Filtered or just Minimised
		cerr << "Unable to read " << inpFile << endl;
		exit(0);
	}

	unsigned int conf = pRotlib_All.size("",residue);

	cout << "size : " << conf << endl;
	unsigned int random = rand() % conf;
	cout << "Random : " << random << endl;
	vector<vector < double > > conformers = pRotlib_All.getInternalCoor("",residue);

	if(!pRotlib_base.addConformation("",residue,conformers[random])) {
		cerr << "Unable to add " << residue <<" Conformation " << endl;
		exit(0);
	}
	cout << "Base size " << pRotlib_base.size("",residue) << endl;
	cout << "Threshold: " << threshold << endl;

	pRotlib_All.removeRotamer("",residue,random);

	rotLoader.setRotamerLibrary(&pRotlib_All);
	if(!rotLoader.loadRotamers(pPosA1,residue,pRotlib_All.size("", pPosA1->getResidueName()))) {
		cerr << "Unable to load rotamers " << endl;
		exit(0);
	}
	rotLoader.setRotamerLibrary(&pRotlib_base);
	if(!rotLoader.loadRotamers(pPosB1,residue,pRotlib_base.size("", residue))) {
		cerr << "Unable to load rotamers " << endl;
		exit(0);
	}
	conformers.erase(conformers.begin() + random);
	int k = 0;

	while( pRotlib_All.size("",residue) > 0 && pRotlib_base.size("",residue) < 5001) {
	
		unsigned int i = rand() % pRotlib_All.size("", residue);

		pPosA1->getCurrentIdentity().setActiveConformation(i);
		AtomPointerVector a1 = pPosA1->getAtomPointers();
		//for(vector<Atom*>::iterator it = a1.begin(); it != a1.end(); it++) {
		//	cout << **it << endl;
		//}			
		//cout << endl;
		
		// for all conformers in chain B (pRotlib_base)
		//rotLoader.loadRotamers(pPosB1,library,pPosB1->getResidueName(),0,pRotlib_base.size(library, pPosB1->getResidueName())-1,true);
		bool add = true;
		//for (unsigned int j  = 0; j < pRotlib_base.size("",residue);j++) {
		for (unsigned int j  = 0; j < pPosB1->getTotalNumberOfRotamers() ;j++) {
			pPosB1->getCurrentIdentity().setActiveConformation(j);
			//cout << j << endl;
			sys.buildAtoms();
			//cout << "Built" << endl;
			AtomPointerVector b1 = pPosB1->getAtomPointers();

		//	for(vector<Atom*>::iterator it = a1.begin(); it != a1.end(); it++) {
		//		cout << **it << endl;
		//	}			
		//	cout << endl;

		//	for(vector<Atom*>::iterator it = b1.begin(); it != b1.end(); it++) {
		//		cout << **it << endl;
		//	}			
		//	cout << endl;
			double trmsd = a1.rmsd(b1);
		//	cout << j << "RMSD: " << trmsd << endl;
	
			if(trmsd < threshold) {
				add = false;
				break;
		      }
		
		} 

		if(add == true) {
			//cout << "Added " << pPosB1->getTotalNumberOfRotamers() << endl;

			if(!pRotlib_base.addConformation("",residue,conformers[i])) {
				cerr << "Unable to add " << residue <<" Conformation " << endl;
				exit(0);
			}
			//cout << "Added to base" << endl;
			rotLoader.addRotamers(pPosB1,residue,pPosB1->getTotalNumberOfRotamers(),pPosB1->getTotalNumberOfRotamers());
			//cout << "Added to chain B" << endl;
		} else {
			//cout << "Dont add" << endl;
		}

		pRotlib_All.removeRotamer("",residue,i);
		conformers.erase(conformers.begin()+i);
		k++;
		//cout << "In " << k++ << "th conformer" << " Size is " << pRotlib_All.size("",residue) << " " << conformers.size() << endl;
		//cout << i << ":Base size " << pRotlib_base.size("",residue) << endl;
	}	

	pRotlib_base.writeFile(outFile);
	
	cout << "Total No of conformers " << pRotlib_All.size("",residue)+1 << endl;
	cout << "Residue:  " << residue << endl;
	cout << "Threshold:  " << threshold << endl;
	cout << "No of representative conformers at " << threshold << " "  << pRotlib_base.size("",residue) << endl;
	cout << "Considered  " << k << " "  << residue << " conformers" << endl;
}

int main(int argc, char* argv[]) {
	if(argc !=5) {
		cerr << "Usage: discardSimilarConformers <inpRotFile ><ResidueName> <threshold> <outRotFile>" << endl;
		cerr << "Reads all rotamers of residue from the database and creates a list of rotamers that differ from each other by atleast the threshold RMSD" << endl;
		exit(0);
	}

	srand(unsigned(time(0)));

	string residue = string(argv[2]);
	discardSimilarConformers(MslTools::toDouble(string(argv[3])),residue,string(argv[1]),string(argv[4]));
	cout << "Done" << endl;
}
