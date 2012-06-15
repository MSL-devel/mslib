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

#include <iostream>
#include <cstdlib>
#include <map>
#include <getopt.h>

#include "AtomContainer.h"
#include "Transforms.h"
#include "System.h"
#include "OptimalRMSDCalculator.h"
#include <time.h>

using namespace std;
using namespace MSL;

double optimalRMSD(AtomContainer& C1, AtomContainer& C2, Transforms& align) {

	if (!align.rmsdAlignment(C1.getAtomPointers(), C2.getAtomPointers())) {
		cerr << "could not superimpose pair of AtomContaners" << endl; exit(-1);
	}
	double rmsd = C1.getAtomPointers().rmsd(C2.getAtomPointers());
	align.revertRmsdAlignment(C1.getAtomPointers(), C2.getAtomPointers(), C1.getAtomPointers());

	return rmsd;
}

// read structure into a system without funny re-ordering stuff
void mySystemReadPDB(System& M, string pdbf) {
	AtomContainer C; C.readPdb(pdbf);
	int ri = 1; vector<int> newResIndices(C.atomSize(), 0);
	for (int i = 0; i < C.atomSize(); i++) {
		if ((i > 0) && ((C[i-1].getResidueNumber() != C[i].getResidueNumber()) || (C[i-1].getResidueName().compare(C[i].getResidueName()) != 0) || (C[i-1].getChainId().compare(C[i].getChainId()) != 0))) ri++;
		newResIndices[i] = ri;
	}
	for (int i = 0; i < C.atomSize(); i++) C[i].setResidueNumber(newResIndices[i]);
	M.addAtoms(C.getAtomPointers());
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cerr << "\n./testOptimalRMSDCalculator [query.pdb] [target.pdb]\n\n"; exit(-1);
	}
	vector<string> bbn; bbn.push_back("N"); bbn.push_back("CA"); bbn.push_back("C"); bbn.push_back("O");

	// read all input structures
	System M1;
	AtomContainer C1;
	mySystemReadPDB(M1, argv[1]);
	for (int pi = 0; pi < M1.positionSize(); pi++) {
		Position& p = M1.getPosition(pi);
		if (p.atomExists("CA")) {
			C1.addAtom(p.getAtom("CA"));
		}
	}
	System M2;
	AtomContainer C2;
	mySystemReadPDB(M2, argv[2]);
	for (int pi = 0; pi < M2.positionSize(); pi++) {
		Position& p = M2.getPosition(pi);
		if (p.atomExists("CA")) {
			C2.addAtom(p.getAtom("CA"));
		}
	}
	// testing time
	double rmsd;
	OptimalRMSDCalculator rr;
	AtomPointerVector& CP1 = C1.getAtomPointers();
	AtomPointerVector& CP2 = C2.getAtomPointers(); 

	clock_t start = clock();
	Transforms align;
	for(int i=0;i<1000000;i++){
		rmsd = optimalRMSD(C1,C2,align);
	}
	clock_t ends = clock();
	cout << "RMSD " << rmsd << " Running Time for MSL_RMSD: " << (double)(ends-start)/CLOCKS_PER_SEC << endl;

	start = clock();
	for(int i=0;i<1000000;i++){
		rmsd = rr.bestRMSD(C1.getAtomPointers(),C2.getAtomPointers());
	}
	ends = clock();
	cout << "RMSD " << rmsd << " Running Time for rmsdKabsch: " << (double)(ends-start)/CLOCKS_PER_SEC << endl;

	AtomPointerVector CP3(CP1);
	start = clock();
	for(int i=0;i<1000000;i++){
		rr.align(CP1,CP2,CP3);
		rmsd = rr.lastRMSD();
	}
	ends = clock();
	cout << "RMSD " << rmsd << " Running Time for rmsdAlignmentKabsch: " << (double)(ends-start)/CLOCKS_PER_SEC << endl;

	AtomContainer C4(C1);
	rr.align(C1.getAtomPointers(),C2.getAtomPointers(),C4.getAtomPointers());
	vector<double> trans = rr.lastTranslation();
	vector<vector<double> > rot = rr.lastRotation();
	cout << trans[0] << " " << trans[0] << " " << trans[0] << endl;
	cout << rot[0][0] << " " << rot[0][1] << " " << rot[0][2] << endl;
	cout << rot[1][0] << " " << rot[1][1] << " " << rot[1][2] << endl;
	cout << rot[2][0] << " " << rot[2][1] << " " << rot[2][2] << endl;
	C4=C1;
	C4.writePdb("align.pdb");
}


