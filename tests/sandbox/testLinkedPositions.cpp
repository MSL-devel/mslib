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



#include <string>
#include <map>
#include <fstream>
#include <cmath>
using namespace std;



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


#include "PolymerSequence.h"
#include "PairwiseEnergyCalculator.h"

using namespace MSL;

#include  "SysEnv.h"
static SysEnv SYSENV;


void test1(vector<double> & _energies, vector<unsigned int> & _counts);
void test2(vector<double> & _energies, vector<unsigned int> & _counts);

int main() {

	/***************************************************
	 *  Test 1 & 2:
	 *
	 *  A: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
	 *  B: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
	 *  C: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN");
	 *
	 *  In test 1 positions 1, 3, 4 and 6 are linked between
	 *  chains.
	 * 
	 *  In test 2 the positions are not linked but the
	 *  states in which the rotamers do not correspond are
	 *  not calculated.
	 * 
	 *  Energies are printed to files:
	 *     /tmp/testLinkedPositions-01.txt
	 *     /tmp/testLinkedPositions-02.txt
	 *  
	 *  PDB files are saved as: 
	 *    /tmp/serial-01-nnnn.pdb
	 *    /tmp/serial-02-nnnnnnnn.pdb
	 *  where nnnn and nnnnnnnn are the states.
	 *  
	 ***************************************************/

	vector<double> e1;
	vector<double> e2;
	vector<unsigned int> c1;
	vector<unsigned int> c2;
	test1(e1, c1);
	test2(e2, c2);

	cout << "Compare the energies of the states (tolerance 1E-10):" << endl;
	for (unsigned int i=0; i<e1.size(); i++) {
		if (abs(e1[i] - e2[i]) < 1E-10) {
			cout << i << " OK" << endl;
		} else {
			cout << i << " NOT OK (diff = " << e1[i]-e2[i] << ")" << endl;
		}
	}

	cout << "Compare the number of interactions in the states:" << endl;
	for (unsigned int i=0; i<c1.size(); i++) {
		if (c1[i] == c2[i]) {
			cout << i << " OK" << endl;
		} else {
			cout << i << " NOT OK (diff = " << c1[i] << " " << c2[i] << ")" << endl;
		}
	}


	return 0;
}


void test1(vector<double> & _energies, vector<unsigned int> & _counts) {


	string baseNum = "01";

	string outfile = "/tmp/testLinkedPositions-" + baseNum + ".txt";
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
A: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
B: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
C: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN");


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

	sys.buildAllAtoms(); // build the active atoms
	AtomSelection as(sys.getAllAtomPointers());
	Transforms tr;

	out_fs << "Translate everything by 16 on X" << endl;
	tr.translate(sys.getAllAtomPointers(), CartesianPoint(16,0,0));

	AtomPointerVector chainB = as.select("chain B", true); // the true is needed to select also the atoms of the inactive identities
	out_fs << "Rotate chain B by 120 around Z" << endl;
	tr.Zrotate(chainB, 120);
	
	AtomPointerVector chainC = as.select("chain C", true);
	out_fs << "Rotate chain C by 240 around Z" << endl;
	tr.Zrotate(chainC, 240);

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

	Position * pPosA3 = &(sys.getPosition("A,3"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosA6 = &(sys.getPosition("A,6"));

	Position * pPosB3 = &(sys.getPosition("B,3"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));

	Position * pPosC3 = &(sys.getPosition("C,3"));
	Position * pPosC4 = &(sys.getPosition("C,4"));
	Position * pPosC6 = &(sys.getPosition("C,6"));

	sysRot.loadRotamers(pPosA3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosA4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosA6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)
	
	sysRot.loadRotamers(pPosB3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosB3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosB4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosB4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosB6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)

	sysRot.loadRotamers(pPosC3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosC3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosC4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosC4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosC6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)

	// link the positions by symmetry
	vector<vector<string> > linkedPositions;
	linkedPositions.push_back(vector<string>());
	linkedPositions.back().push_back("A,1");
	linkedPositions.back().push_back("B,1");
	linkedPositions.back().push_back("C,1");
	linkedPositions.push_back(vector<string>());
	linkedPositions.back().push_back("A,3");
	linkedPositions.back().push_back("B,3");
	linkedPositions.back().push_back("C,3");
	linkedPositions.push_back(vector<string>());
	linkedPositions.back().push_back("A,4");
	linkedPositions.back().push_back("B,4");
	linkedPositions.back().push_back("C,4");
	linkedPositions.push_back(vector<string>());
	linkedPositions.back().push_back("A,6");
	linkedPositions.back().push_back("B,6");
	linkedPositions.back().push_back("C,6");
	sys.setLinkedPositions(linkedPositions);

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
			out_fs << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers, and is of type " << pPos->getLinkedPositionType() << endl;
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
	SPM.saveInteractionCounts(true);
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

	_energies.clear();
	_counts.clear();
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

		_energies.push_back(SPM.getStateEnergy(state));
		_counts.push_back(SPM.getStateInteractionCount(state));
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


void test2(vector<double> & _energies, vector<unsigned int> & _counts) {


	string baseNum = "02";

	string outfile = "/tmp/testLinkedPositions-" + baseNum + ".txt";
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
A: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
B: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN\n\
C: [ALA GLY] ARG [ASN LYS] [ILE ASP] CYS GLU GLN");


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

	out_fs << "Translate everything by 16 on X" << endl;
	tr.translate(sys.getAllAtomPointers(), CartesianPoint(16,0,0));

	AtomPointerVector chainB = as.select("chain B", true); // the true is needed to select also the atoms of the inactive identities
	out_fs << "Rotate chain B by 120 around Z" << endl;
	tr.Zrotate(chainB, 120);
	
	AtomPointerVector chainC = as.select("chain C", true);
	out_fs << "Rotate chain C by 240 around Z" << endl;
	tr.Zrotate(chainC, 240);

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

	Position * pPosA3 = &(sys.getPosition("A,3"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosA6 = &(sys.getPosition("A,6"));

	Position * pPosB3 = &(sys.getPosition("B,3"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));

	Position * pPosC3 = &(sys.getPosition("C,3"));
	Position * pPosC4 = &(sys.getPosition("C,4"));
	Position * pPosC6 = &(sys.getPosition("C,6"));

	sysRot.loadRotamers(pPosA3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosA4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosA4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosA6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)
	
	sysRot.loadRotamers(pPosB3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosB3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosB4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosB4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosB6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)

	sysRot.loadRotamers(pPosC3, "ASN", 2); // 2 ASN rotamers at A 3 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosC3, "LYS", 2); // 3 LYS rotamers at A 3    "      "      "        "        "
	sysRot.loadRotamers(pPosC4, "ILE", 1); // 3 ILE rotamers at A 4 (rotamer and identity variable position)
	sysRot.loadRotamers(pPosC4, "ASP", 3); // 3 ASP rotamers at A 4    "      "      "        "        "
	sysRot.loadRotamers(pPosC6, "GLU", 4); // 4 GLU rotamers at B 4 (rotamer only variable position)

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
			out_fs << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers, and is of type " << pPos->getLinkedPositionType() << endl;
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
	SPM.saveInteractionCounts(true);
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

	_energies.clear();
	_counts.clear();
	for (unsigned int i=0; i<stateEnum.size(); i++) {
		// for (unsigned int i=0; i<2; i++) {
		vector<unsigned int> state = stateEnum[i];
		if (state[0] != state[4] || state[0] != state[8] || state[1] != state[5] || state[1] != state[9] || state[2] != state[6] || state[2] != state[10] || state[3] != state[7] || state[3] != state[11]) {
			/***************************************************************
			 *  This IF statement ensures to compute only the states in 
			 *  which the positions that are linked in test1 have the 
			 *  same rotamers
			 ***************************************************************/
			continue;
		}
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

		_energies.push_back(SPM.getStateEnergy(state));
		_counts.push_back(SPM.getStateInteractionCount(state));
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
