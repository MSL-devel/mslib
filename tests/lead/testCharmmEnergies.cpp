/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
Kulp DW et al. "Structural informatics, modeling and design with a open 
source Molecular Software Library (MSL)" (2012) J. Comp. Chem, in press
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
#include "OnTheFlyManager.h"
#include "Transforms.h"
using namespace std;

using namespace MSL;

#include "SysEnv.h"
static SysEnv SYSENV;

#include "MslOut.h"
static MslOut MSLOUT("testCharmmEnergies");


//const double EPS = 0.000000000000001; // Test will fail with this Epsilon 05/24/12 dwkulp
const double EPS = 0.000001;

void test0();
void test1();
void test2();
void test3();
void test4_onthefly();
void test4_preCompute();
void setupSystem(PolymerSequence &seq, System &sys, int _numRotamers, bool _rotamersOnMultipleIdentitiesOnly=true);
void runTest(System &sys, string _testName, bool _enumerateStates=true);

int main() {

	/***************************************************

	  TEST 0
	  Simple System #0
             2 Positions variable, 1 in idenitity and 1 in conformation only.

	*****************************************************/
        test0();


	/***************************************************

	  TEST 1
	  Simple System #1
             3 Positions variable, 2 in idenitity and 1 in conformation only.

	*****************************************************/
	test1();


	/***************************************************

	  TEST 2
	  Simple System #2 
               1 Variable Position at terminal of chain, 2 conseuctive variable positions

	*****************************************************/
	test2();



	/***************************************************

	  TEST 3 
	  Complex System #1 
               85 variable positions, 20 rotamers each 

         The system has 3 chains, 133 positions, 143 identities, 563 rotamers, 2225 active atoms, 2390 total atoms

	*****************************************************/
	test3();


	/***************************************************

	  TEST 4 
             Time test , computing energy for a System with the following:


	*****************************************************/
	test4_onthefly();
	test4_preCompute();


	cout << "LEAD"<<endl;
	return 0;
}


void test0(){
	cout << " ***********************     TEST 1 ****************************"<<endl;



		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple rotamers 
		 *
		 *************************************************************/
		System sys;
		PolymerSequence seq("\
A: ILE\n\
B: LEU\n");


		setupSystem(seq,sys,2,false);


		runTest(sys, "test0");

}
void test1(){


	cout << " ***********************     TEST 1 ****************************"<<endl;

		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA),
		 *   each has 2 rotamers.
		 *************************************************************/
		PolymerSequence seq("\
A: ALA ARG ASN [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE [LEU ALA] LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");
	

		setupSystem(seq,sys,2);


		runTest(sys, "test1");
}


void test2(){



	cout << " ***********************     TEST 2 ****************************"<<endl;
		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ASN TYR) and (A 5 ILE ASP)
		 *
		 *************************************************************/
		PolymerSequence seq("\
A: ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: [SER THR] HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL");
	
		setupSystem(seq,sys,3);


		runTest(sys, "test2");
}


void test3(){



	cout << " ***********************     TEST 3 ****************************"<<endl;
		System sys;

		PolymerSequence seq("\
A: ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET\n\
C: MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL");
	

		
		setupSystem(seq,sys,2);

		runTest(sys, "test3",false); // Do not enumerate the states, exceeds max limit of Enumerator::calcEnumeration
}


void test4_onthefly(){

		Timer t;
		double start = t.getWallTime();	

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
B: GLY HSE [ILE PHE ARG GLU] LEU LYS [ARG LYS GLU ASP] GLY HSE [ILE LEU PHE VAL] LEU LYS MET GLY [GLN ASN SER THR] ILE LEU [LYS ARG GLU ASP] MET GLY HSE [ILE LEU VAL] LEU LYS MET [GLY PHE TRP TYR] HSE ILE LEU [LYS ARG GLU ASP] [MET TRP PHE TYR] GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET\n\
C: [MET ASN] [PHE TRP] PRO [SER THR TYR] THR [TRP PHE TYR] TYR [VAL ARG LEU] MET [PHE TRP] PRO [SER THR] THR TRP [TYR PHE] VAL [MET ILE LEU] PHE PRO SER [ASN GLN THR] TRP TYR [ARG LYS ASN] MET PHE PRO [SER THR] THR TRP [TYR PHE] VAL MET PHE PRO [SER THR] THR TRP TYR [VAL ILE LEU] MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL");
	

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




		OnTheFlyManager pec(&sys,SYSENV.getEnv("MSL_CHARMM_PAR"));
		start = t.getWallTime();
		pec.calculateTotalEnergy(sys);
		fprintf(stdout,"Test4 ON-THE-FLY TIME: %8.3f\n",t.getWallTime() - start);
	
}

void test4_preCompute(){

		Timer t;
		double start = t.getWallTime();	

		cout << " ***********************     TEST 4 ****************************"<<endl;
		System sys;

		/*************************************************************
		 *
		 *   Create a sequence from scratch with a couple positions
		 *   having multiple identities (A 4 ILE-ASP and B 4 LEU-ALA
		 *
		 *****************************n********************************/
		// [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ]
		PolymerSequence seq("\
A: ALA [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ASN TYR] [ILE ASP] CYS GLU [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN ALA ARG [ASN TYR] [ILE ASP] CYS GLU GLN\n\
B: GLY HSE [ILE PHE ARG GLU] LEU LYS [ARG LYS GLU ASP] GLY HSE [ILE LEU PHE VAL] LEU LYS MET GLY [GLN ASN SER THR] ILE LEU [LYS ARG GLU ASP] MET GLY HSE [ILE LEU VAL] LEU LYS MET [GLY PHE TRP TYR] HSE [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [LYS ARG GLU ASP] [MET TRP PHE TYR] GLY HSE ILE LEU LYS MET GLY HSE ILE LEU LYS MET\n\
C: [MET ARG ASN] [PHE TRP] PRO [SER THR TYR] THR [TRP PHE TYR] TYR [VAL ARG LEU] MET [PHE TRP] PRO [SER THR] THR TRP [TYR PHE] VAL [MET ILE LEU] PHE PRO SER [ASN GLN THR] TRP TYR [ARG LYS ASN] MET PHE PRO [SER THR] THR TRP [TYR PHE] VAL MET PHE PRO [SER THR] THR TRP TYR [VAL ILE LEU] MET PHE PRO SER THR TRP TYR VAL MET PHE PRO SER THR TRP TYR VAL [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ] [ CYS ASP GLU PHE HSD HSE ILE LYS LEU MET ASN GLN ARG SER THR VAL TRP TYR ]");

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
		CSB.setBuildNonBondedInteractions(true); // BUILD non-bonded terms.
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
		int numTotalRotamers = 0;
		int numTotalIdentities = 0;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				numTotalIdentities += pPos->identitySize();
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
					numTotalRotamers += pRes->getNumberOfAltConformations();
				}
			}
		}
		cout << "The systems has " << sys.chainSize() << " chains, " <<  sys.positionSize() << " positions, " << numTotalIdentities<<" identities, "<<numTotalRotamers<<" rotamers, "<< sys.atomSize() << " active atoms, " << sys.allAtomSize() << " total atoms" << endl;

		
		start = t.getWallTime();
		sys.calcEnergy();
		fprintf(stdout,"TEST4 ENE-SET TIME: %8.3f\n",t.getWallTime() - start);
		sys.printEnergySummary();
}


void setupSystem(PolymerSequence &seq, System &sys, int _numRotamers, bool _rotamersOnMultipleIdentitiesOnly){


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
		for (uint c = 0; c < sys.chainSize();c++){
		  if (!sys.seed(MslTools::stringf("%s 1 C",sys.getChain(c).getChainId().c_str()).c_str(),MslTools::stringf("%s 1 CA",sys.getChain(c).getChainId().c_str()).c_str(),MslTools::stringf("%s 1 N",sys.getChain(c).getChainId().c_str()).c_str())){
		      cerr << "Cannot seed atoms C, CA, N on residue 1 " <<sys.getChain(c).getChainId()<< endl;
		    }
		}

		sys.buildAtoms(); // build the active atoms
		AtomSelection as(sys.getAllAtomPointers());
		Transforms tr;
		if (sys.chainExists("B")){
		    AtomPointerVector chainB = as.select("chain B");
		    cout << "Select and translate chain B by (5, 4, 9)" << endl;
		    cout << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
		    tr.translate(chainB, CartesianPoint(5, 4, 9));
		}

		if (sys.chainExists("C")){
		  AtomPointerVector chainC = as.select("chain C");
		  cout << "Select and translate chain C by (-5, -10, -8)" << endl;
		  cout << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
		  //chainC.translate(CartesianPoint(-5, -10, -8));
		  tr.translate(chainC, CartesianPoint(-5, -10, -8));
		}
		
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
			if (!_rotamersOnMultipleIdentitiesOnly || posVar->getNumberOfIdentities() > 1){
				for (uint j = 0; j < posVar->getNumberOfIdentities();j++){
					string identityName = posVar->getIdentity(j).getResidueName();
					if (identityName != "GLY" || identityName != "PRO")
					  sysRot.loadRotamers(posVar, identityName, 0, _numRotamers+1, "BALANCED-200"); 
				}
			}
		}




		/*************************************************************
		 *
		 *   List the chain, positions and identities of the system
		 *
		 *************************************************************/
		int numTotalRotamers = 0;
		int numTotalIdentities = 0;
		int numTotalVariablePositions = 0;
		for (unsigned int i=0; i<sys.chainSize(); i++) {
			Chain * pChain = &(sys.getChain(i));
			cout << "  Chain " << pChain->getChainId() << " has " << pChain->positionSize() << " positions, " << pChain->atomSize() << " active atoms, " << pChain->allAtomSize() << " total atoms" << endl;
			for (unsigned int j=0; j<pChain->positionSize(); j++) {
				Position * pPos = &(pChain->getPosition(j));
				cout << "     Position " << pPos->getResidueNumber() << " has " << pPos->identitySize() << " identities, " << pPos->atomSize() << " active atoms, " << pPos->allAtomSize() << " total atoms, and " << pPos->getTotalNumberOfRotamers() << " total number of rotamers" << endl;
				numTotalIdentities += pPos->identitySize();
				if (pPos->getTotalNumberOfRotamers() > 1){
				  numTotalVariablePositions++;
				}
				for (unsigned int k=0; k<pPos->identitySize(); k++) {
					Residue * pRes = &(pPos->getIdentity(k));
					cout << "        Identity " << pRes->getResidueName() << " has " << pRes->size() << " atoms, and " << pRes->getNumberOfAltConformations() << " rotamers" << endl;
					numTotalRotamers += pRes->getNumberOfAltConformations();

				}
			}
		}
		fprintf(stdout,"\nThe system has:\n\t%6d chains\n\t%6d positions\n\t%6d variable positions\n\t%6d identities\n\t%6d rotamers\n\t%6d active atoms\n\t%6d total atoms\n\n",
			sys.chainSize(), sys.positionSize(),numTotalVariablePositions,numTotalIdentities,numTotalRotamers,sys.atomSize() ,sys.allAtomSize());


}

void runTest(System &sys, string _testName, bool _enumerateStates){


		SelfPairManager SPM(&sys);
	        OnTheFlyManager pec(&sys, SYSENV.getEnv("MSL_CHARMM_PAR"));


		if (_enumerateStates){

		  Timer t;
		  double start = t.getWallTime();	

		  SPM.saveEnergiesByTerm(true);
		  //SPM.setOnTheFly(true);
		  SPM.calculateEnergies();

		  fprintf(stdout,"SPM TIME: %8.3f\n",t.getWallTime() - start);



		  
		  start = t.getWallTime();	
		  pec.calculateEnergyTable(sys);

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
		    sys.printEnergySummary();

		    char c[1000];
		    sprintf(c, "/tmp/serial-%04u.pdb", i);
		    sys.writePdb(c);
		    cout << "Written pdb " << c << endl;
		    cout << " - - - - - - - " << endl;
		    cout << endl;

		    // Now do it Dan's way...
		    cout << "ON-THE-FLY WAY"<<endl;

		    // Define a local 
		    OnTheFlyManager pecState(&sys,SYSENV.getEnv("MSL_CHARMM_PAR"));
		    pecState.getCharmmEnergyCalculator()->setCollectNumberOfInteractions(true);
		    pecState.getCharmmEnergyCalculator()->setEnergyByType(true);
		    pecState.calculateStateEnergy(sys,state);
		    pecState.printSummary();

		    double e = pec.getStateEnergy(sys, state);
		    //fprintf(stdout, "\nPRE-COMPUTED WAY ( using energy table ): %8.3f\n\n\n",e);


		    /*
		      Pass/fail
		    */
		    double onTheFly_calc    = pec.getStateEnergy(sys, state);
		    double onTheFly_lookup  = pecState.calculateStateEnergy(sys,state);
		    double SPM_lookup       = SPM.getStateEnergy(state);
		    double SYS_calc         = sys.calcEnergy();


		    if (fabs(SYS_calc-SPM_lookup) > EPS) {
		      fprintf(stderr,"FAILURE %s 0001 %8.12f %8.12f %8.12f\n",_testName.c_str(),SYS_calc,SPM_lookup,abs(SYS_calc-SPM_lookup));
		      exit(1001);
		    }

		    if (fabs(onTheFly_calc-onTheFly_lookup) > EPS) {
		      fprintf(stderr,"FAILURE %s 0002 %8.12f %8.12f %8.12f\n",_testName.c_str(),onTheFly_calc,onTheFly_lookup,abs(onTheFly_calc-onTheFly_lookup));
		      exit(1002);
		    }

		    if (fabs(onTheFly_calc-SPM_lookup) > EPS) {
		      fprintf(stderr,"FAILURE %s 0003 %8.12f %8.12f %8.12f\n",_testName.c_str(),onTheFly_calc,SPM_lookup,abs(onTheFly_calc-SPM_lookup));
		      exit(1003);
		    }

		    if (fabs(onTheFly_calc-SYS_calc) > EPS) {
		      fprintf(stderr,"FAILURE %s 0004 %8.12f %8.12f %8.12f\n",_testName.c_str(),onTheFly_calc,SYS_calc,abs(onTheFly_calc-SYS_calc));
		      exit(1004);
		    }
			
		  }
		} else {
		  cout << "Calculate total energy ON-THE-FLY"<<endl;
		  double onTheFly_calc    = pec.calculateTotalEnergy(sys);
		  cout << "Calculate total energy SYS-CALC"<<endl;
		  double SYS_calc         = sys.calcEnergy();

		  if (fabs(onTheFly_calc-SYS_calc) > EPS) {
		      fprintf(stderr,"FAILURE %s 0005 %8.12f %8.12f %8.12f\n",_testName.c_str(),onTheFly_calc,SYS_calc,abs(onTheFly_calc-SYS_calc));
		      exit(1005);
		  }
		}




}
