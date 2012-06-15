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

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "AtomSelection.h"
#include "SelfPairManager.h"
#include "Enumerator.h"
#include "Timer.h"
#include "PairwiseEnergyCalculator.h"
#include "Transforms.h"
using namespace std;

using namespace MSL;


#include "SysEnv.h"
static SysEnv SYSENV;

int main() {

	Timer t;
	double start = t.getWallTime();	
	System sys;

	/*************************************************************
	 *
	 *   Create a sequence from scratch with a couple positions
	 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
	 *
	 *************************************************************/
	PolymerSequence seq("A: ALA ARG ASN ILE CYS GLU GLN           B: GLY HSE ILE LEU LYS MET          C: MET PHE PRO SER THR TRP TYR VAL");


	/*************************************************************
	 *
	 *   Create a system from the sequence, seed the chains, and
	 *   build the chains, translate them so that they do not overlap, 
	 *
	 *************************************************************/
	cout << "Sequence:" << endl;
	cout << seq.toString();
	cout << endl;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	cout << "Use toppar " << topFile << ", " << parFile << endl;

	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.setBuildNonBondedInteractions(true); // Don't build non-bonded terms.
	CSB.buildSystem(seq);
	cout << endl;

	// seed the chains (give cartesian coordinates to 3 atoms in each)
	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}

	if (!sys.seed("C 1 C", "C 1 CA", "C 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}

	// build the active atoms
	sys.buildAtoms(); 
	AtomSelection as(sys.getAllAtomPointers());


	Transforms tr;

	// translate chains B and C so that they do not clash
	AtomPointerVector chainB = as.select("chain B");
	cout << "Select and translate chain B by (13, 4, 9)" << endl;
	cout << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
//	chainB.translate(CartesianPoint(13, 4, 9));
	tr.translate(chainB, CartesianPoint(13, 4, 9));
	
	AtomPointerVector chainC = as.select("chain C");
	cout << "Select and translate chain C by (-5, -10, -8)" << endl;
	cout << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
//	chainC.translate(CartesianPoint(-5, -10, -8));
	tr.translate(chainC, CartesianPoint(-5, -10, -8));

	// write a PDB of the system
	string filename = "/tmp/initialBuild.pdb";
	cout << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();



	/************************************************
	 *   S T A R T : ADD ROTAMERS
	 ************************************************/
	cout << "Add rotamers to the system" << endl;
	string rotlib = SYSENV.getEnv("MSL_ROTLIB");
	cout << "Read rotamer library " << rotlib << endl;
	cout << endl;
	SystemRotamerLoader sysRot(sys, rotlib);

	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosC2 = &(sys.getPosition("C,2"));

	sysRot.loadRotamers(pPosA4, "BALANCED-200", "ILE", 0, 0); // ILE rotamers at A 4 (rotamer and identity variable position)
	//sysRot.loadRotamers(pPosA4, "BALANCED-200", "ASP", 0, 2); // ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosB4, "BALANCED-200", "LEU", 0, 0); // LEU rotamers at B 4 (identity only variable position)
	//sysRot.loadRotamers(pPosB4, "BALANCED-200", "ALA", 0, 0); // ALA rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosC2, "BALANCED-200", "PHE", 0, 0); // PHE rotamers at C 2 (rotamer only variable position)

	/*************************************************************
	 *
	 *   List the chain, positions and identities of the system
	 *
	 ************************************************************* /
	cout << "The systems has " << sys.size() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
	for (unsigned int i=0; i<sys.size(); i++) {
		Chain * pChain = &(sys.getChain(i));
		cout << "  Chain " << pChain->getChainId() << " has " << pChain->size() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
		for (unsigned int j=0; j<pChain->size(); j++) {
			Position * pPos = &(pChain->getPositionByIndex(j));
			cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->size() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
			for (unsigned int k=0; k<pPos->size(); k++) {
				Residue * pRes = &(pPos->getIdentity(k));
				cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
			}
		}
	}
	*/




	fprintf(stdout,"DONE BUILDING: %8.3f\n",t.getWallTime() - start);

	sys.calcEnergy();
	sys.printEnergySummary();

	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	// this should rebuild the full table and look like the one before
	CSB.updateNonBonded();

	sys.calcEnergy();
	sys.printEnergySummary();
	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	// this should rebuild a shorter table with a 8A cutoff
	cout << "Calculate the energies with cutoffs: on 8.0 - off 9.0 - list exclusion 12.0" << endl;
	CSB.updateNonBonded( 8.0, 9.0, 12.0);

	sys.calcEnergy();
	sys.printEnergySummary();
	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	// this should rebuild a shorter table with a 8A cutoff
	cout << "Calculate the energies with cutoffs: on 8.0 - off 9.0 - list exclusion 11.0" << endl;
	CSB.updateNonBonded(8.0, 9.0, 11.0);

	sys.calcEnergy();
	sys.printEnergySummary();
	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	// this should rebuild a shorter table with a 8A cutoff
	cout << "Calculate the energies with cutoffs: on 8.0 - off 9.0 - list exclusion 10.0" << endl;
	CSB.updateNonBonded( 8.0, 9.0, 10.0);

	sys.calcEnergy();
	sys.printEnergySummary();
	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	// this should rebuild a shorter table with a 8A cutoff
	cout << "Calculate the energies with cutoffs: on 8.0 - off 9.0 - list exclusion 9.0" << endl;
	CSB.updateNonBonded(8.0, 9.0, 9.0);

	sys.calcEnergy();
	sys.printEnergySummary();
	fprintf(stdout,"DONE CALCULATING ENERGY: %8.3f\n",t.getWallTime() - start);

	return 0;
}


