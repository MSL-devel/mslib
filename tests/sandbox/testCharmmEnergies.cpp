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




void test1();
void test2();
void test3();
void test4();


int main() {

	/***************************************************

	  TEST 1
	  Simple System #1
             3 Positions variable, 2 in idenitity and 1 in conformation only.

	*****************************************************/
	test1();


	/***************************************************

	  TEST 2
	  Simple System #2 
               1 Variable Position at terminal of chain 2 conseuctive variable positions

	*****************************************************/
	//test2();



	/***************************************************

	  TEST 3 
	  Complex System #1 
               85 variable positions, 2 rotamers each (I had 100 rotamers each, but too large)
               use SelfPairManager to compute energy table

	*****************************************************/
	//test3();


	/***************************************************

	  TEST 3 
	  Complex System #1 
               85 variable positions, 2 rotamers each (I had 100 rotamers each, but too large)
               use PairwiseEnergyCalculator to compute energy table

	*****************************************************/
	//test4();


	
	return 0;
}


void test1(){


	cout << " ***********************     TEST 1 ****************************"<<endl;

		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 *************************************************************/
		PolymerSequence seq("\
A: ALA ARG ASN [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE [LEU ALA] LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");
	

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

		sysRot.loadRotamers(pPosA4, "ILE", 0, 3); // ILE rotamers at A 4 (rotamer and identity variable position)
		sysRot.loadRotamers(pPosA4, "ASP", 0, 2); // ASP rotamers at A 4    "      "      "        "        "
		sysRot.loadRotamers(pPosB4, "LEU", 0, 0); // LEU rotamers at B 4 (identity only variable position)
		sysRot.loadRotamers(pPosB4, "ALA", 0, 0); // ALA rotamers at B 4    "        "      "       "
		sysRot.loadRotamers(pPosC2, "PHE", 0, 3); // PHE rotamers at C 2 (rotamer only variable position)

		/*************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 *************************************************************/
		cout << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
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


		PairwiseEnergyCalculator pec(SYSENV.getEnv("MSL_CHARMM_PAR"));

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
			sprintf(c, "/tmp/serial-%04u.pdb", i);
			sys.writePdb(c);
			cout << "Written pdb " << c << endl;
			cout << " - - - - - - - " << endl;
			cout << endl;

			// Now do it Dan's way...
			cout << "DANS ON-THE-FLY WAY"<<endl;

			// Define a local 
			PairwiseEnergyCalculator pecState(SYSENV.getEnv("MSL_CHARMM_PAR"));
			
		
			pecState.calculateStateEnergy(sys,state);
			pecState.printSummary();

			double e = pec.getStateEnergy(sys, state);
			fprintf(stdout, "\nDANS PRE-COMPUTED WAY ( using energy table ): %8.3f\n\n\n",e);



		}
}


void test2(){



	cout << " ***********************     TEST 2 ****************************"<<endl;
		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 *************************************************************/
		PolymerSequence seq("\
A: ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");
	

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
		//chainB.translate(CartesianPoint(13, 4, 9));
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

		Position * pPosA3 = &(sys.getPosition("A,3"));
		Position * pPosA4 = &(sys.getPosition("A,4"));
		Position * pPosC1 = &(sys.getPosition("C,1"));

		sysRot.loadRotamers(pPosA3, "ASN", 0, 3); // ASN rotamers at A 3 (rotamer and identity variable position)
		sysRot.loadRotamers(pPosA3, "TYR", 0, 2); // TYR rotamers at A 3    "      "      "        "        "
		sysRot.loadRotamers(pPosA4, "ILE", 0, 3); // ILE rotamers at A 4 (rotamer and identity variable position)
		sysRot.loadRotamers(pPosA4, "ASP", 0, 2); // ASP rotamers at A 4    "      "      "        "        "
		sysRot.loadRotamers(pPosC1, "MET", 0, 3); // MET rotamers at C 1 (rotamer only variable position)

		/*************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 *************************************************************/
		cout << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
				}
			}
		}

	
		SelfPairManager SPM(&sys);
		SPM.calculateEnergies();

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
			sprintf(c, "/tmp/serial-%04u.pdb", i);
			sys.writePdb(c);
			cout << "Written pdb " << c << endl;
			cout << " - - - - - - - " << endl;
			cout << endl;

			// Now do it Dan's way...
			cout << "DANS WAY"<<endl;

			PairwiseEnergyCalculator pec(SYSENV.getEnv("MSL_CHARMM_PAR"));
			pec.calculateStateEnergy(sys,state);

			pec.printSummary();


		}
}


void test3(){



	cout << " ***********************     TEST 3 ****************************"<<endl;
		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 *************************************************************/
		PolymerSequence seq("\
A: ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL");
	

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
		CSB.setCreatePairwiseTable(false);
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
		//chainB.translate(CartesianPoint(13, 4, 9));
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
		for (uint i = 0; i < sys.positionSize();i++){
			Position * posVar = &(sys.getPosition(i));
			if (posVar->getNumberOfIdentities() > 1){
				for (uint j = 0; j < posVar->getNumberOfIdentities();j++){
					string identityName = posVar->getIdentity(j).getResidueName();
					sysRot.loadRotamers(posVar, identityName, 0, 100); // Load 100 rotamers
				}
			}
		}

		/*************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 *************************************************************/
		cout << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
				}
			}
		}

		Timer t;
		double start = t.getWallTime();	
		SelfPairManager SPM(&sys);
		SPM.calculateEnergies();	
		fprintf(stdout,"SPM TIME: %8.3f\n",t.getWallTime() - start);
		int foo =0;
		cin >> foo;

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

		start = t.getWallTime();
		//for (unsigned int i=0; i<stateEnum.size(); i++) {
		  for (unsigned int i=0; i<0; i++) {
			vector<unsigned int> state = stateEnum[i];
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

		}
		fprintf(stdout,"STATE TIME: %8.3f\n",t.getWallTime() - start);
		cin >> foo;

}


void test4(){



	cout << " ***********************     TEST 4 ****************************"<<endl;
		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 *****************************n********************************/
		PolymerSequence seq("\
A: ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL");
	

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
		CSB.setCreatePairwiseTable(false);
		CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
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
		//chainB.translate(CartesianPoint(13, 4, 9));
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
		for (uint i = 0; i < sys.positionSize();i++){
			Position * posVar = &(sys.getPosition(i));
			if (posVar->getNumberOfIdentities() > 1){
				for (uint j = 0; j < posVar->getNumberOfIdentities();j++){
					string identityName = posVar->getIdentity(j).getResidueName();
					sysRot.loadRotamers(posVar, identityName, 0, 100, "BALANCED-200"); // Load 100 rotamers
				}
			}
		}

		/*************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 *************************************************************/
		cout << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
				}
			}
		}

		Timer t;
		double start = t.getWallTime();	

		PairwiseEnergyCalculator pec(SYSENV.getEnv("MSL_CHARMM_PAR"));
		pec.calculateTotalEnergy(sys);
		pec.printSummary();

		fprintf(stdout,"PEC TIME: %8.3f\n",t.getWallTime() - start);

		
		int foo =0;
		cin >> foo;

}
			/*

			map<Interaction*,int> &pecMap = pec.getInteractions();

			map<Interaction*,int>::iterator pecInt;

			map<string, vector<Interaction *> > *intMap = sys.getEnergySet()->getEnergyTerms();
			for (uint i = 0; i < (*intMap)["CHARMM_VDW"].size();i++){
				if (!(*intMap)["CHARMM_VDW"][i]->isActive()) continue;
		    
				//if ( ((*intMap)["CHARMM_VDW"][i]->getAtomPointers()[0]->getResidueNumber() == 1 && (*intMap)["CHARMM_VDW"][i]->getAtomPointers()[0]->getChainId() == "A") ||
				//     ((*intMap)["CHARMM_VDW"][i]->getAtomPointers()[1]->getResidueNumber() == 1 && (*intMap)["CHARMM_VDW"][i]->getAtomPointers()[1]->getChainId() == "A")){
				//cout << (*intMap)["CHARMM_VDW"][i]<<" "<<(*intMap)["CHARMM_VDW"][i]->toString()<<endl;
				//}

		
				pecInt = pecMap.find((*intMap)["CHARMM_VDW"][i]);
				if (pecInt == pecMap.end()){
					cout << " - "<<(*intMap)["CHARMM_VDW"][i]<<" "<<(*intMap)["CHARMM_VDW"][i]->toString()<<endl;
					pecMap[(*intMap)["CHARMM_VDW"][i]] = -1;
				} else {
					pecMap[(*intMap)["CHARMM_VDW"][i]] = -2;
			
				}
			}

			for (pecInt = pecMap.begin(); pecInt != pecMap.end();pecInt++){
				if (pecInt->second == 1){
					cout << " + "<<pecInt->first<<" "<<pecInt->first->toString()<<endl;
				}
			}
			*/
