/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

void test1();
void test2();

int main() {

	/***************************************************
	 *  Test 1 & 2:
	 *
	 *  A: [ALA GLY] ARG ASN [ILE(2) ASP(2)] CYS GLU GLN
	 *  B: GLY HSE ILE [LEU ALA] LYS [MET SER]
	 *  C: MET PHE(2) PRO SER THR TRP [TYR THR] VAL
	 *
	 *  5 positions with double (or triple) identity (and
	 *  a position (C2) with multiple rotamers).
	 *
	 *  In test 1 all identities are created when the System
	 *  is built
	 *  In test 2 the additional identities are added with a
	 *  call to CharmmSystemBuilder::addIdentity
	 *
	 *  Energies are computed using the SelfPairManager and
	 *  with direct calculation.
	 *
	 *  Diff to verify output files 
	 *    /tmp/testAddCharmmIdentities-01.txt and
	 *    /tmp/testAddCharmmIdentities-02.txt
	 *  
	 ***************************************************/
	test1();
	test2();

	return 0;
}


void test1(){

	string baseNum = "01";

	string outfile = "/tmp/testAddCharmmIdentities-" + baseNum + ".txt";
	ofstream out_fs;
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}

	out_fs << " ***********************     TEST 1 ****************************"<<endl;

	System sys;

	/*************************************************************
	 *
	 *   Create a sequence from scratch with a couple positions
	 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
	 *
	 *************************************************************/
	PolymerSequence seq("\
A: [ALA GLY] ARG ASN [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE [LEU ALA LYS] LYS [MET SER]\n\
C: MET PHE PRO SER THR TRP [TYR THR] VAL");


	/*************************************************************
	 *
	 *   Create a system from the sequence, seed the chains, and
	 *   build the chains, translate them so that they do not overlap, 
	 *
	 *************************************************************/
	out_fs << "Sequence:" << endl;
	out_fs << seq.toString();
	out_fs << endl;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	out_fs << "Use toppar " << topFile << ", " << parFile << endl;

	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.buildSystem(seq);
	out_fs << endl;

	if (!sys.seed()) {
		cerr << "cannot seed the system" << endl;
	}
	/*
	if (!sys.seed("A,1,C", "A,1,CA", "A,1,N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	if (!sys.seed("B,1,C", "B,1,CA", "B,1,N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}
	if (!sys.seed("C,1,C", "C,1,CA", "C,1,N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 C" << endl;
	}
	*/

	//sys.buildAtoms(); // build the active atoms
	sys.buildAllAtoms(); // build the active atoms
	AtomSelection as(sys.getAllAtomPointers());
	Transforms tr;

	AtomPointerVector chainB = as.select("chain B", true); // the true is needed to select also the atoms of the inactive identities
	out_fs << "Select and translate chain B by (13, 4, 9)" << endl;
	out_fs << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
	tr.translate(chainB, CartesianPoint(13, 4, 9));
	
	AtomPointerVector chainC = as.select("chain C", true);
	out_fs << "Select and translate chain C by (-5, -10, -8)" << endl;
	out_fs << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
	tr.translate(chainC, CartesianPoint(-5, -10, -8));

/*
	// copy the backbone atoms from the active identties to all alternative
	// and build the inactive identities
	vector<string> bbAtoms;
	bbAtoms.push_back("N");
	bbAtoms.push_back("HN");
	bbAtoms.push_back("CA");
	bbAtoms.push_back("C");
	bbAtoms.push_back("O");
	sys.copyCoordinatesOfAtomsInPosition(bbAtoms);
	sys.buildAllAtoms();
*/
	string filename = "/tmp/initialBuild-" + baseNum + ".pdb";
	out_fs << "Write pdb " << filename << endl;
	if (!sys.writePdb(filename)) {
		cerr << "ERROR writing " << filename << endl;
		return;
	}

	/************************************************
	 *   S T A R T : ADD ROTAMERS
	 ************************************************/
	out_fs << "Add rotamers to the system" << endl;
	string rotlib = SYSENV.getEnv("MSL_ROTLIB");
	out_fs << "Read rotamer library " << rotlib << endl;
	out_fs << endl;
	SystemRotamerLoader sysRot(sys, rotlib);

	Position * pPosA1 = &(sys.getPosition("A,1"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));
	Position * pPosC2 = &(sys.getPosition("C,2"));
	Position * pPosC7 = &(sys.getPosition("C,7"));

	sysRot.loadRotamers(pPosA1, "ALA", 1); // 1 ALA rotamers at A 1 (identity only variable position)
	sysRot.loadRotamers(pPosA4, "ILE", 3); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosB4, "LEU", 1); // 1 LEU rotamers at B 4 (identity only variable position)
	sysRot.loadRotamers(pPosB4, "ALA", 1); // 1 ALA rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB4, "LYS", 1); // 1 LYS rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB6, "MET", 1); // 1 MET rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosB6, "SER", 1); // 1 SER rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosC2, "PHE", 3); // 3 PHE rotamers at C 2 (rotamer only variable position)
	sysRot.loadRotamers(pPosC7, "TYR", 1); // 1 TYR rotamers at C 7 (identity only variable position)
	sysRot.loadRotamers(pPosC7, "THR", 1); // 1 THR rotamers at C 7 (identity only variable position)

	/*************************************************************
	 *
	 *   List the chain, positions and identities of the system
	 *
	 *************************************************************/
	out_fs << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
	for (unsigned int i=0; i<sys.chainSize(); i++) {
		Chain * pChain = &(sys.getChain(i));
		out_fs << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
		for (unsigned int j=0; j<pChain->positionSize(); j++) {
			Position * pPos = &(pChain->getPosition(j));
			out_fs << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
			for (unsigned int k=0; k<pPos->identitySize(); k++) {
				Residue * pRes = &(pPos->getIdentity(k));
				out_fs << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
			}
		}
	}

	// check if all atoms have been built
	AtomSelection sel(sys.getAllAtomPointers());
	sel.select("noCoor, HASCOOR 0");
	if (sel.selectionSize("noCoor") > 0) {
		AtomPointerVector noCoorAtoms = sel.getSelection("noCoor");
		cerr << "The following atoms do not have coordinates" << endl;
		cerr << noCoorAtoms;
		exit(1);
	}

	SelfPairManager SPM(&sys);
	SPM.saveEnergiesByTerm(true);
	SPM.calculateEnergies();

	vector<unsigned int> rots = SPM.getNumberOfRotamers();
	out_fs << "The system has " << rots.size() << " variable positions" << endl;
	for (unsigned int i=0; i<rots.size(); i++) {
		out_fs << "   Variable position " << i << " has " << rots[i] << " rotamers" << endl;
	}

	out_fs << endl;
	out_fs << "START COMPARING SELF & PAIR ENERGIES WITH CALCULATED ENERGIES" << endl;
	out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
	out_fs << endl;

	Enumerator stateEnum(rots);

	for (unsigned int i=0; i<stateEnum.size(); i++) {
		// for (unsigned int i=0; i<2; i++) {
		vector<unsigned int> state = stateEnum[i];
		out_fs << "--------------- NEW STATE TO COMPUTE ENERGY OF -----------------"<<endl;
		out_fs << "Self & Pair state energy of";
		for (unsigned int j=0; j<state.size(); j++) {
			out_fs << " " << state[j];
		}
		out_fs << " = " << SPM.getStateEnergy(state) << endl;
		vector<string> stateDescriptors = SPM.getStateDescriptors(state);
		vector<vector<unsigned int> > stateIndeces = SPM.getStatePositionIdentityRotamerIndeces(state);
		for (unsigned int j=0; j<stateDescriptors.size(); j++) {
			out_fs << " * " << stateDescriptors[j] << " -- " << stateIndeces[j][0] << "/" << stateIndeces[j][1] << "/" << stateIndeces[j][2] << endl;
		}

		out_fs << SPM.getSummary(state);

		sys.setActiveRotamers(state);
		out_fs << endl;
		out_fs << "Direct calculation of energy (serial " << i << ") = " << sys.calcEnergy() << endl;
		out_fs << sys.getEnergySummary();

		char c[1000];
		sprintf(c, "/tmp/serial-%2s-%04u.pdb", baseNum.c_str(), i);
		sys.writePdb(c);
		out_fs << "Written pdb " << c << endl;
		out_fs << " - - - - - - - " << endl;
		out_fs << endl;

	}

	out_fs.close();
	cout << "Test " << baseNum << ": created output file " << outfile << endl;
}

void test2(){

	string baseNum = "02";

	string outfile = "/tmp/testAddCharmmIdentities-" + baseNum + ".txt";
	ofstream out_fs;
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}

	out_fs << " ***********************     TEST 2 ****************************"<<endl;

	System sys;

	/*************************************************************
	 *
	 *   Create a sequence from scratch with a couple positions
	 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
	 *
	 *************************************************************/
	PolymerSequence seq("\
A: ALA ARG ASN ILE CYS GLU GLN\n\
B: GLY HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");


	/*************************************************************
	 *
	 *   Create a system from the sequence, seed the chains, and
	 *   build the chains, translate them so that they do not overlap, 
	 *
	 *************************************************************/
	out_fs << "Sequence:" << endl;
	out_fs << seq.toString();
	out_fs << endl;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	out_fs << "Use toppar " << topFile << ", " << parFile << endl;

	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.buildSystem(seq);
	out_fs << endl;

	if (!sys.seed()) {
		cerr << "cannot seed the system" << endl;
	}
	/*
	if (!sys.seed("A,1,C", "A,1,CA", "A,1,N")) {
		cerr << "cannot seed atoms c, ca, n on residue 1 a" << endl;
	}
	if (!sys.seed("B,1,C", "B,1,CA", "B,1,N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}
	if (!sys.seed("C,1,C", "C,1,CA", "C,1,N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 C" << endl;
	}
	*/

	sys.buildAtoms(); // build the active atoms
	AtomSelection as(sys.getAllAtomPointers());
	Transforms tr;

	AtomPointerVector chainB = as.select("chain B");
	out_fs << "Select and translate chain B by (13, 4, 9)" << endl;
	out_fs << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
	tr.translate(chainB, CartesianPoint(13, 4, 9));
	
	AtomPointerVector chainC = as.select("chain C");
	out_fs << "Select and translate chain C by (-5, -10, -8)" << endl;
	out_fs << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
	tr.translate(chainC, CartesianPoint(-5, -10, -8));

	// copy the backbone atoms from the active identties to all alternative
	// and build the inactive identities
	/*
	vector<string> bbAtoms;
	bbAtoms.push_back("N");
	bbAtoms.push_back("HN");
	bbAtoms.push_back("CA");
	bbAtoms.push_back("C");
	bbAtoms.push_back("O");
	sys.copyCoordinatesOfAtomsInPosition(bbAtoms);
	*/
	sys.buildAllAtoms();

	string filename = "/tmp/initialBuild-" + baseNum + ".pdb";
	out_fs << "Write pdb " << filename << endl;
	if (!sys.writePdb(filename)) {
		cerr << "ERROR writing " << filename << endl;
		return;
	}

	/************************************************
	 *   S T A R T : ADD ROTAMERS
	 ************************************************/
	out_fs << "Add rotamers to the system" << endl;
	string rotlib = SYSENV.getEnv("MSL_ROTLIB");
	out_fs << "Read rotamer library " << rotlib << endl;
	out_fs << endl;
	SystemRotamerLoader sysRot(sys, rotlib);

	Position * pPosA1 = &(sys.getPosition("A,1"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));
	Position * pPosC2 = &(sys.getPosition("C,2"));
	Position * pPosC7 = &(sys.getPosition("C,7"));

	CSB.addIdentity("A,1", "GLY");
	CSB.addIdentity("A,4", "ASP");
	vector<string> posB4Ids;
	posB4Ids.push_back("ALA");
	posB4Ids.push_back("LYS");
	CSB.addIdentity("B,4", posB4Ids);
	CSB.addIdentity("B,6", "SER");
	CSB.addIdentity("C,7", "THR");
	//exit(0);
	//sys.buildAllAtoms();
	CSB.updateNonBonded();

	sysRot.loadRotamers(pPosA1, "ALA", 1); // 1 ALA rotamers at A 1 (identity only variable position)
	sysRot.loadRotamers(pPosA4, "ILE", 3); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosB4, "LEU", 1); // 1 LEU rotamers at B 4 (identity only variable position)
	sysRot.loadRotamers(pPosB4, "ALA", 1); // 1 ALA rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB4, "LYS", 1); // 1 LYS rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB6, "MET", 1); // 1 MET rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosB6, "SER", 1); // 1 SER rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosC2, "PHE", 3); // 3 PHE rotamers at C 2 (rotamer only variable position)
	sysRot.loadRotamers(pPosC7, "TYR", 1); // 1 TYR rotamers at C 7 (identity only variable position)
	sysRot.loadRotamers(pPosC7, "THR", 1); // 1 THR rotamers at C 7 (identity only variable position)

	/*************************************************************
	 *
	 *   List the chain, positions and identities of the system
	 *
	 *************************************************************/
	out_fs << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
	for (unsigned int i=0; i<sys.chainSize(); i++) {
		Chain * pChain = &(sys.getChain(i));
		out_fs << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
		for (unsigned int j=0; j<pChain->positionSize(); j++) {
			Position * pPos = &(pChain->getPosition(j));
			out_fs << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
			for (unsigned int k=0; k<pPos->identitySize(); k++) {
				Residue * pRes = &(pPos->getIdentity(k));
				out_fs << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
			}
		}
	}

	SelfPairManager SPM(&sys);
	SPM.saveEnergiesByTerm(true);
	SPM.calculateEnergies();

	vector<unsigned int> rots = SPM.getNumberOfRotamers();
	out_fs << "The system has " << rots.size() << " variable positions" << endl;
	for (unsigned int i=0; i<rots.size(); i++) {
		out_fs << "   Variable position " << i << " has " << rots[i] << " rotamers" << endl;
	}

	out_fs << endl;
	out_fs << "START COMPARING SELF & PAIR ENERGIES WITH CALCULATED ENERGIES" << endl;
	out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
	out_fs << endl;

	Enumerator stateEnum(rots);

	for (unsigned int i=0; i<stateEnum.size(); i++) {
		// for (unsigned int i=0; i<2; i++) {
		vector<unsigned int> state = stateEnum[i];
		out_fs << "--------------- NEW STATE TO COMPUTE ENERGY OF -----------------"<<endl;
		out_fs << "Self & Pair state energy of";
		for (unsigned int j=0; j<state.size(); j++) {
			out_fs << " " << state[j];
		}
		out_fs << " = " << SPM.getStateEnergy(state) << endl;
		vector<string> stateDescriptors = SPM.getStateDescriptors(state);
		vector<vector<unsigned int> > stateIndeces = SPM.getStatePositionIdentityRotamerIndeces(state);
		for (unsigned int j=0; j<stateDescriptors.size(); j++) {
			out_fs << " * " << stateDescriptors[j] << " -- " << stateIndeces[j][0] << "/" << stateIndeces[j][1] << "/" << stateIndeces[j][2] << endl;
		}

		out_fs << SPM.getSummary(state);

		sys.setActiveRotamers(state);
		out_fs << endl;
		out_fs << "Direct calculation of energy (serial " << i << ") = " << sys.calcEnergy() << endl;
		out_fs << sys.getEnergySummary();

		char c[1000];
		sprintf(c, "/tmp/serial-%2s-%04u.pdb", baseNum.c_str(), i);
		sys.writePdb(c);
		out_fs << "Written pdb " << c << endl;
		out_fs << " - - - - - - - " << endl;
		out_fs << endl;

	}

	out_fs.close();
	cout << "Test " << baseNum << ": created output file " << outfile << endl;
}



/*
void test2(){


	cout << " ***********************     TEST 2 ****************************"<<endl;

		System sys;

		/ *************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 ************************************************************* /
		PolymerSequence seq("\
A: ALA ARG ASN ILE CYS GLU GLN\n\
B: GLY HSE ILE [LEU ALA] LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");
	

		/ *************************************************************
		 *
		 *   Create a system from the sequence, seed the chains, and
		 *   build the chains, translate them so that they do not overlap, 
		 *
		 ************************************************************* /
		cout << "Sequence:" << endl;
		cout << seq.toString();
		cout << endl;
		string topFile = "/library/charmmTopPar/top_all22_prot.inp";
		string parFile = "/library/charmmTopPar/par_all22_prot.inp";
		cout << "Use toppar " << topFile << ", " << parFile << endl;

		CharmmSystemBuilder CSB(sys, topFile, parFile);
		//CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
		CSB.buildSystem(seq);
		cout << endl;

		if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
			cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
		}
		if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
			cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
		}

		if (!sys.seed("C 1 C", "C 1 CA", "C 1 N")) {
			cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
		}

		sys.buildAtoms(); // build the active atoms
		AtomSelection as(sys.getAllAtomPointers());
		Transforms tr;

		AtomPointerVector chainB = as.select("chain B");
		cout << "Select and translate chain B by (13, 4, 9)" << endl;
		cout << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
	//	chainB.translate(CartesianPoint(13, 4, 9));
		tr.translate(chainB, CartesianPoint(13, 4, 9));
		
		AtomPointerVector chainC = as.select("chain C");
		cout << "Select and translate chain C by (-5, -10, -8)" << endl;
		cout << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
		//chainC.translate(CartesianPoint(-5, -10, -8));
		tr.translate(chainC, CartesianPoint(-5, -10, -8));

		// copy the backbone atoms from the active identties to all alternative
		// and build the inactive identities
		vector<string> bbAtoms;
		bbAtoms.push_back("N");
		bbAtoms.push_back("HN");
		bbAtoms.push_back("CA");
		bbAtoms.push_back("C");
		bbAtoms.push_back("O");
		sys.copyCoordinatesOfAtomsInPosition(bbAtoms);
		cout << endl;
		sys.buildAllAtoms();
		cout << endl;

		string filename = "/tmp/initialBuild-2.pdb";
		cout << "Write pdb " << filename << endl;
		if (!sys.writePdb(filename)) {
			cerr << "ERROR writing " << filename << endl;
			return;
		}


		/ ************************************************
		 *   S T A R T : ADD ROTAMERS
		 ************************************************ /
		cout << "Add rotamers to the system" << endl;
		string rotlib = "/library/rotlib/balanced/rotlib-balanced-200.txt";
		cout << "Read rotamer library " << rotlib << endl;
		cout << endl;
		SystemRotamerLoader sysRot(sys, rotlib);

		Position * pPosA4 = &(sys.getPosition("A,4"));
		Position * pPosB4 = &(sys.getPosition("B,4"));
		Position * pPosC2 = &(sys.getPosition("C,2"));

		CSB.addIdentity(*pPosA4, "ASP", bbAtoms);
		sys.buildAllAtoms();
		CSB.updateNonBonded();

		sysRot.loadRotamers(pPosA4, "BALANCED-200", "ILE", 0, 3); // ILE rotamers at A 4 (rotamer and identity variable position)
		sysRot.loadRotamers(pPosA4, "BALANCED-200", "ASP", 0, 2); // ASP rotamers at A 4    "      "      "        "        "
		sysRot.loadRotamers(pPosB4, "BALANCED-200", "LEU", 0, 0); // LEU rotamers at B 4 (identity only variable position)
		sysRot.loadRotamers(pPosB4, "BALANCED-200", "ALA", 0, 0); // ALA rotamers at B 4    "        "      "       "
		sysRot.loadRotamers(pPosC2, "BALANCED-200", "PHE", 0, 3); // PHE rotamers at C 2 (rotamer only variable position)

		/ *************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 ************************************************************* /
		cout << "The systems has " << sys.size() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
		for (unsigned int i=0; i<sys.size(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->size() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->size(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->size() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				for (unsigned int k=0; k<pPos->size(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
				}
			}
		}

		Timer t;

	
		SelfPairManager SPM(&sys);

		double start = t.getWallTime();	
		SPM.calculateEnergies();

		fprintf(stdout,"SPM TIME: %8.3f\n",t.getWallTime() - start);


		PairwiseEnergyCalculator pec("/library/charmmTopPar/par_all27_prot_lipid.inp");

		start = t.getWallTime();	

		// Calculate the complete energy table.
		pec.calculateEnergyTable(sys);

		fprintf(stdout,"PEC TIME: %8.3f\n",t.getWallTime() - start);


		vector<unsigned int> rots = SPM.getNumberOfRotamers();
		cout << "The system has " << rots.size() << " variable positions" << endl;
		for (unsigned int i=0; i<rots.size(); i++) {
			cout << "   Variable position " << i << " has " << rots[i] << " rotamers" << endl;
		}

		cout << endl;
		cout << "START COMPARING SELF & PAIR ENERGIES WITH CALCULATED ENERGIES" << endl;
		cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		cout << endl;
	
		Enumerator stateEnum(rots);

		for (unsigned int i=0; i<stateEnum.size(); i++) {
			// for (unsigned int i=0; i<2; i++) {
			vector<unsigned int> state = stateEnum[i];
			cout << "--------------- NEW STATE TO COMPUTE ENERGY OF -----------------"<<endl;
			cout << "Self & Pair state energy of";
			for (unsigned int j=0; j<state.size(); j++) {
				cout << " " << state[j];
			}
			cout << " = " << SPM.getStateEnergy(state) << endl;
			vector<string> stateDescriptors = SPM.getStateDescriptors(state);
			vector<vector<unsigned int> > stateIndeces = SPM.getStatePositionIdentityRotamerIndeces(state);
			for (unsigned int j=0; j<stateDescriptors.size(); j++) {
				cout << " * " << stateDescriptors[j] << " -- " << stateIndeces[j][0] << "/" << stateIndeces[j][1] << "/" << stateIndeces[j][2] << endl;
			}

			cout << SPM.getSummary(state);

			sys.setActiveRotamers(state);
			cout << endl;
			cout << "Direct calculation of energy (serial " << i << ") = " << sys.calcEnergy() << endl;
			(sys.getEnergySet())->printSummary();

			char c[1000];
			sprintf(c, "/tmp/serial-2-%04u.pdb", i);
			sys.writePdb(c);
			cout << "Written pdb " << c << endl;
			cout << " - - - - - - - " << endl;
			cout << endl;

			// Now do it Dan's way...
			cout << "DANS ON-THE-FLY WAY"<<endl;

			// Define a local 
			PairwiseEnergyCalculator pecState("/library/charmmTopPar/par_all27_prot_lipid.inp");
			
		
			pecState.calculateStateEnergy(sys,state);
			pecState.printSummary();

			double e = pec.getStateEnergy(sys, state);
			fprintf(stdout, "\nDANS PRE-COMPUTED WAY ( using energy table ): %8.3f\n\n\n",e);



		}
}
*/

