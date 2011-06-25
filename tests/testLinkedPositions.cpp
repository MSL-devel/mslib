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



#include <string>
#include <map>
#include <fstream>
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

/*
  Results from last run (5/13/09)

0,2,4 = resi 6  (chains A,B,C)
1,3,5 = resi 9  (chains A,B,C)

M1,M2     = Master positions
S1-X,S2-X = Slave positions

Compare UNLINKED energy table to LINKED energy table

Precision of addition outside program is 0.000 this gives small errors in .000 or .00 

SuperRotamer 0 - SuperRotamer 0
1.     0      0      1      0   -3.516 M1 to M2
2.     0      0      3      0   -0.062 M1 to S2-0
3.     0      0      5      0   -0.003 M1 to S2-1
4.     1      0      2      0   -0.003 M2 to S1-0
5.     1      0      4      0   -0.062 M2 to S1-1
6.     2      0      3      0   -3.515 S1-0 to S2-0
7.     2      0      5      0   -0.062 S1-0 to S2-1
8.     3      0      4      0   -0.003 S2-0 to S1-1
9.     4      0      5      0   -3.517 S1-1 to S2-1

 Total:                 -10.743
 Linked:                -10.745

SuperRotamer 0 - SuperRotamer 1
1.     0      0      1      1    4.695 M1 to M2
2.     0      0      3      1   -0.051 M1 to S2-0
3.     0      0      5      1   -0.005 M1 to S2-1
4.     1      1      2      0   -0.005 M2 to S1-0
5.     1      1      4      0   -0.051 M2 to S2-1
6.     2      0      3      1    4.620 S1-0 to S2-0
7.     2      0      5      1   -0.051 S1-0 to S2-1
8.     3      1      4      0   -0.005 S2-0 to S1-1
9.     4      0      5      1    4.674 S1-1 to S2-1
 Total:                  13.821
 Linked:                 13.819




Better Precision :
SuperRotamer 0 - SuperRotamer 0
1.     0      0      1      0 -3.516221 M1 to M2	   
2.     0      0      3      0 -0.061966	M1 to S2-0  
3.     0      0      5      0 -0.003346	M1 to S2-1  
4.     1      0      2      0 -0.003274	M2 to S1-0  
5.     1      0      4      0 -0.061934	M2 to S2-1  
x6.     2      0      3      0 -3.515245	S1-0 to S2-0
7.     2      0      5      0 -0.061833	S1-0 to S2-1
8.     3      0      4      0 -0.003467	S2-0 to S1-1
9.     4      0      5      0 -3.517256	S1-1 to S2-1


Total :            -10.744542
Linked:            -10.744543


 */

void test1(vector<double> & _energies, vector<unsigned int> & _counts);
void test2(vector<double> & _energies, vector<unsigned int> & _counts);

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

	vector<double> e1;
	vector<double> e2;
	vector<unsigned int> c1;
	vector<unsigned int> c2;
	test1(e1, c1);
	test2(e2, c2);

	for (unsigned int i=0; i<e1.size(); i++) {
		if (e1[i] == e2[i]) {
			cout << i << " E OK" << endl;
		} else {
			cout << i << " E NOT OK " << e1[i]-e2[i] << endl;
		}
	}

	for (unsigned int i=0; i<c1.size(); i++) {
		if (c1[i] == c2[i]) {
			cout << i << " Count OK" << endl;
		} else {
			cout << i << " Count NOT OK " << c1[i] << " " << c2[i] << endl;
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

//	Position * pPosA1 = &(sys.getPosition("A,1"));
	Position * pPosA3 = &(sys.getPosition("A,3"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosA6 = &(sys.getPosition("A,6"));

//	Position * pPosB1 = &(sys.getPosition("B,1"));
	Position * pPosB3 = &(sys.getPosition("B,3"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));

//	Position * pPosC1 = &(sys.getPosition("C,1"));
	Position * pPosC3 = &(sys.getPosition("C,3"));
	Position * pPosC4 = &(sys.getPosition("C,4"));
	Position * pPosC6 = &(sys.getPosition("C,6"));

	//sysRot.loadRotamers(pPosA1, "ALA", 1); // 1 ALA rotamers at A 1 (identity only variable position)
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

	/*
	sysRot.loadRotamers(pPosB4, "ALA", 1); // 1 ALA rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB4, "LYS", 1); // 1 LYS rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB6, "MET", 1); // 1 MET rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosB6, "SER", 1); // 1 SER rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosC2, "PHE", 3); // 3 PHE rotamers at C 2 (rotamer only variable position)
	sysRot.loadRotamers(pPosC7, "TYR", 1); // 1 TYR rotamers at C 7 (identity only variable position)
	sysRot.loadRotamers(pPosC7, "THR", 1); // 1 THR rotamers at C 7 (identity only variable position)
	*/

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


	/*
	writePdbFile();


	System initialSystem;
	initialSystem.readPdb("/tmp/symmetricTrimer.pdb");

	stringstream seq;
	for (uint c = 0; c< initialSystem.chainSize();c++){
		Chain ch = initialSystem.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.positionSize();p++){
			Position & pos = ch.getPosition(p);

			cout << pos <<endl;
			string chainId = pos.getChainId();
			int resNum     = pos.getResidueNumber();


			if (resNum == 9 || resNum == 6){
				seq << " [ VAL THR ]";
			} else {

				seq << " "<<pos.getCurrentIdentity().getResidueName();
			}

		}
		seq << "\n";
		
	}

	cout <<seq.str() << endl;


	PolymerSequence pseq(seq.str());

	cout << pseq << endl;
	

	System sys;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	cout << "Use toppar " << topFile << ", " << parFile << endl;
	CharmmSystemBuilder CSB(sys,topFile,parFile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);


	sys.assignCoordinates(initialSystem.getAtomPointers(),false);

	sys.buildAllAtoms(); 

	string filename = "/tmp/initialBuild.pdb";
	cout << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();


	SystemRotamerLoader sysRot(sys, SYSENV.getEnv("MSL_ROTLIB"));
		
	for (uint i = 0; i < sys.positionSize();i++){
		Position * posVar = &(sys.getPosition(i));

		if (posVar->getTotalNumberOfRotamers() > 1){
			sysRot.loadRotamers(posVar, "BALANCED-200", "VAL", 0, 2); 
			sysRot.loadRotamers(posVar, "BALANCED-200", "THR", 0, 2); 

			cout << "Position: "<<posVar->getChainId()<<" "<<posVar->getResidueNumber()<<" has "<<posVar->getTotalNumberOfRotamers()<<endl;
		}
	}


	PairwiseEnergyCalculator pec(parFile);
	pec.calculateEnergyTable(sys);

	vector<vector<double> > selfEnergy = pec.getSelfTable();
	vector<vector<double> > templateEnergy = pec.getTemplateTable();
	vector<vector<vector<vector<double> > > > pairEnergy = pec.getPairTable();


	
	ofstream eout;
	eout.open("/tmp/energyTable.unlinked.txt");
	double fixedE = 0.0;
	int variableIndexI = 0;
	for (uint i = 0; i< selfEnergy.size();i++){
		for (uint j = 0; j < selfEnergy[i].size();j++){
			if (selfEnergy[i].size() > 1){
				char t[40];
				sprintf(t, "%6d %6d %8.6f\n", variableIndexI,j,(selfEnergy[i][j]+templateEnergy[i][j]));
				eout << t;
			} else {
				fixedE += (selfEnergy[i][j]+templateEnergy[i][j]);
			}
		}

		if (selfEnergy[i].size() > 1){
			variableIndexI++;
		}
	}

	variableIndexI = 0;
	for (uint i = 0; i< pairEnergy.size();i++){
		if (selfEnergy[i].size() <= 1) continue;

		for (uint j = 0; j < pairEnergy[i].size();j++){
			int variableIndexJ = 0;
			for (uint ii = 0; ii < pairEnergy[i][j].size();ii++){
				if (selfEnergy[ii].size() <= 1) continue;
				for (uint jj = 0; jj < pairEnergy[i][j][ii].size();jj++){
					char t[100];
					sprintf(t, "%6d %6d %6d %6d %8.6f\n", variableIndexI,j,variableIndexJ,jj,pairEnergy[i][j][ii][jj]);
					eout << t;
				}

				variableIndexJ++;
			}
		}
		variableIndexI++;
	}
	char t[40];
	sprintf(t, "Fixed: %8.3f\n",fixedE);
	eout << t;
	eout.close();




	vector<vector<Position *> > linkedPos;
	linkedPos.resize(2);
	for (uint i = 0; i < sys.positionSize();i++){
		Position * posVar = &(sys.getPosition(i));

		if (posVar->getTotalNumberOfRotamers() > 1){

			if (posVar->getResidueNumber() == 6){
				linkedPos[0].push_back(posVar);
			}

			if (posVar->getResidueNumber() == 9){
				linkedPos[1].push_back(posVar);
			}
			
		}

	}

	/ *
	linkedPos[0][0]->setLinkedPositionType(Position::MASTER);
	for (uint i = 0; i < linkedPos[0].size();i++){

		if (i > 0){
			linkedPos[0][i]->setLinkedPositionType(Position::SLAVE);
		}

		for (uint j = 0;j < linkedPos[0].size();j++){
			if (i == j) continue;

			linkedPos[0][i]->addLinkedPosition(*linkedPos[0][j]);
		}
	}
	linkedPos[1][0]->setLinkedPositionType(Position::MASTER);
	for (uint i = 0; i < linkedPos[1].size();i++){

		if (i > 0){
			linkedPos[1][i]->setLinkedPositionType(Position::SLAVE);
		}

		for (uint j = 0;j < linkedPos[1].size();j++){
			if (i == j) continue;

			linkedPos[1][i]->addLinkedPosition(*linkedPos[1][j]);
		}
	}
	* /
	for (uint i = 1; i < linkedPos[0].size();i++){
		linkedPos[0][0]->addLinkedPosition(*linkedPos[0][i]);
	}
	for (uint i = 1; i < linkedPos[1].size();i++){
		linkedPos[1][0]->addLinkedPosition(*linkedPos[1][i]);
	}



	for (uint i = 0; i < sys.positionSize();i++){
		cout << sys.getPosition(i).getChainId()<<" " <<sys.getPosition(i).getResidueNumber()<<" "<<sys.getPosition(i).getLinkedPositionType()<<endl;
	}

	PairwiseEnergyCalculator pecLinked(parFile);
	pecLinked.calculateEnergyTable(sys);

	selfEnergy     = pecLinked.getSelfTable();
	templateEnergy = pecLinked.getTemplateTable();
	pairEnergy     = pecLinked.getPairTable();

	eout.open("/tmp/energyTable.linked.txt");
	fixedE = 0.0;
	variableIndexI = 0;
	for (uint i = 0; i< selfEnergy.size();i++){
		for (uint j = 0; j < selfEnergy[i].size();j++){
			if (selfEnergy[i].size() > 1){
				char t[40];
				sprintf(t, "%6d %6d %8.6f\n", variableIndexI,j,(selfEnergy[i][j]+templateEnergy[i][j]));
				eout << t;
			} else {
				fixedE += (selfEnergy[i][j]+templateEnergy[i][j]);
			}
		}

		if (selfEnergy[i].size() > 1){
			variableIndexI++;
		}
	}

	variableIndexI = 0;
	for (uint i = 0; i< pairEnergy.size();i++){
		if (selfEnergy[i].size() <= 1) continue;

		for (uint j = 0; j < pairEnergy[i].size();j++){
			int variableIndexJ = 0;
			for (uint ii = 0; ii < pairEnergy[i][j].size();ii++){
				if (selfEnergy[ii].size() <= 1) continue;
				for (uint jj = 0; jj < pairEnergy[i][j][ii].size();jj++){
					char t[100];
					sprintf(t, "%6d %6d %6d %6d %8.6f\n", variableIndexI,j,variableIndexJ,jj,pairEnergy[i][j][ii][jj]);
					eout << t;
				}

				variableIndexJ++;
			}
		}
		variableIndexI++;
	}
	char g[40];
	sprintf(g, "Fixed: %8.3f\n",fixedE);
	eout << g;
	eout.close();

	*/	

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

//	Position * pPosA1 = &(sys.getPosition("A,1"));
	Position * pPosA3 = &(sys.getPosition("A,3"));
	Position * pPosA4 = &(sys.getPosition("A,4"));
	Position * pPosA6 = &(sys.getPosition("A,6"));

//	Position * pPosB1 = &(sys.getPosition("B,1"));
	Position * pPosB3 = &(sys.getPosition("B,3"));
	Position * pPosB4 = &(sys.getPosition("B,4"));
	Position * pPosB6 = &(sys.getPosition("B,6"));

//	Position * pPosC1 = &(sys.getPosition("C,1"));
	Position * pPosC3 = &(sys.getPosition("C,3"));
	Position * pPosC4 = &(sys.getPosition("C,4"));
	Position * pPosC6 = &(sys.getPosition("C,6"));

	//sysRot.loadRotamers(pPosA1, "ALA", 1); // 1 ALA rotamers at A 1 (identity only variable position)
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

	/*
	sysRot.loadRotamers(pPosB4, "ALA", 1); // 1 ALA rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB4, "LYS", 1); // 1 LYS rotamers at B 4    "        "      "       "
	sysRot.loadRotamers(pPosB6, "MET", 1); // 1 MET rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosB6, "SER", 1); // 1 SER rotamers at B 6    "        "      "       "
	sysRot.loadRotamers(pPosC2, "PHE", 3); // 3 PHE rotamers at C 2 (rotamer only variable position)
	sysRot.loadRotamers(pPosC7, "TYR", 1); // 1 TYR rotamers at C 7 (identity only variable position)
	sysRot.loadRotamers(pPosC7, "THR", 1); // 1 THR rotamers at C 7 (identity only variable position)
	*/

	/*
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
	*/

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


	/*
	writePdbFile();


	System initialSystem;
	initialSystem.readPdb("/tmp/symmetricTrimer.pdb");

	stringstream seq;
	for (uint c = 0; c< initialSystem.chainSize();c++){
		Chain ch = initialSystem.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.positionSize();p++){
			Position & pos = ch.getPosition(p);

			cout << pos <<endl;
			string chainId = pos.getChainId();
			int resNum     = pos.getResidueNumber();


			if (resNum == 9 || resNum == 6){
				seq << " [ VAL THR ]";
			} else {

				seq << " "<<pos.getCurrentIdentity().getResidueName();
			}

		}
		seq << "\n";
		
	}

	cout <<seq.str() << endl;


	PolymerSequence pseq(seq.str());

	cout << pseq << endl;
	

	System sys;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	cout << "Use toppar " << topFile << ", " << parFile << endl;
	CharmmSystemBuilder CSB(sys,topFile,parFile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);


	sys.assignCoordinates(initialSystem.getAtomPointers(),false);

	sys.buildAllAtoms(); 

	string filename = "/tmp/initialBuild.pdb";
	cout << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();


	SystemRotamerLoader sysRot(sys, SYSENV.getEnv("MSL_ROTLIB"));
		
	for (uint i = 0; i < sys.positionSize();i++){
		Position * posVar = &(sys.getPosition(i));

		if (posVar->getTotalNumberOfRotamers() > 1){
			sysRot.loadRotamers(posVar, "BALANCED-200", "VAL", 0, 2); 
			sysRot.loadRotamers(posVar, "BALANCED-200", "THR", 0, 2); 

			cout << "Position: "<<posVar->getChainId()<<" "<<posVar->getResidueNumber()<<" has "<<posVar->getTotalNumberOfRotamers()<<endl;
		}
	}


	PairwiseEnergyCalculator pec(parFile);
	pec.calculateEnergyTable(sys);

	vector<vector<double> > selfEnergy = pec.getSelfTable();
	vector<vector<double> > templateEnergy = pec.getTemplateTable();
	vector<vector<vector<vector<double> > > > pairEnergy = pec.getPairTable();


	
	ofstream eout;
	eout.open("/tmp/energyTable.unlinked.txt");
	double fixedE = 0.0;
	int variableIndexI = 0;
	for (uint i = 0; i< selfEnergy.size();i++){
		for (uint j = 0; j < selfEnergy[i].size();j++){
			if (selfEnergy[i].size() > 1){
				char t[40];
				sprintf(t, "%6d %6d %8.6f\n", variableIndexI,j,(selfEnergy[i][j]+templateEnergy[i][j]));
				eout << t;
			} else {
				fixedE += (selfEnergy[i][j]+templateEnergy[i][j]);
			}
		}

		if (selfEnergy[i].size() > 1){
			variableIndexI++;
		}
	}

	variableIndexI = 0;
	for (uint i = 0; i< pairEnergy.size();i++){
		if (selfEnergy[i].size() <= 1) continue;

		for (uint j = 0; j < pairEnergy[i].size();j++){
			int variableIndexJ = 0;
			for (uint ii = 0; ii < pairEnergy[i][j].size();ii++){
				if (selfEnergy[ii].size() <= 1) continue;
				for (uint jj = 0; jj < pairEnergy[i][j][ii].size();jj++){
					char t[100];
					sprintf(t, "%6d %6d %6d %6d %8.6f\n", variableIndexI,j,variableIndexJ,jj,pairEnergy[i][j][ii][jj]);
					eout << t;
				}

				variableIndexJ++;
			}
		}
		variableIndexI++;
	}
	char t[40];
	sprintf(t, "Fixed: %8.3f\n",fixedE);
	eout << t;
	eout.close();




	vector<vector<Position *> > linkedPos;
	linkedPos.resize(2);
	for (uint i = 0; i < sys.positionSize();i++){
		Position * posVar = &(sys.getPosition(i));

		if (posVar->getTotalNumberOfRotamers() > 1){

			if (posVar->getResidueNumber() == 6){
				linkedPos[0].push_back(posVar);
			}

			if (posVar->getResidueNumber() == 9){
				linkedPos[1].push_back(posVar);
			}
			
		}

	}

	/ *
	linkedPos[0][0]->setLinkedPositionType(Position::MASTER);
	for (uint i = 0; i < linkedPos[0].size();i++){

		if (i > 0){
			linkedPos[0][i]->setLinkedPositionType(Position::SLAVE);
		}

		for (uint j = 0;j < linkedPos[0].size();j++){
			if (i == j) continue;

			linkedPos[0][i]->addLinkedPosition(*linkedPos[0][j]);
		}
	}
	linkedPos[1][0]->setLinkedPositionType(Position::MASTER);
	for (uint i = 0; i < linkedPos[1].size();i++){

		if (i > 0){
			linkedPos[1][i]->setLinkedPositionType(Position::SLAVE);
		}

		for (uint j = 0;j < linkedPos[1].size();j++){
			if (i == j) continue;

			linkedPos[1][i]->addLinkedPosition(*linkedPos[1][j]);
		}
	}
	* /
	for (uint i = 1; i < linkedPos[0].size();i++){
		linkedPos[0][0]->addLinkedPosition(*linkedPos[0][i]);
	}
	for (uint i = 1; i < linkedPos[1].size();i++){
		linkedPos[1][0]->addLinkedPosition(*linkedPos[1][i]);
	}



	for (uint i = 0; i < sys.positionSize();i++){
		cout << sys.getPosition(i).getChainId()<<" " <<sys.getPosition(i).getResidueNumber()<<" "<<sys.getPosition(i).getLinkedPositionType()<<endl;
	}

	PairwiseEnergyCalculator pecLinked(parFile);
	pecLinked.calculateEnergyTable(sys);

	selfEnergy     = pecLinked.getSelfTable();
	templateEnergy = pecLinked.getTemplateTable();
	pairEnergy     = pecLinked.getPairTable();

	eout.open("/tmp/energyTable.linked.txt");
	fixedE = 0.0;
	variableIndexI = 0;
	for (uint i = 0; i< selfEnergy.size();i++){
		for (uint j = 0; j < selfEnergy[i].size();j++){
			if (selfEnergy[i].size() > 1){
				char t[40];
				sprintf(t, "%6d %6d %8.6f\n", variableIndexI,j,(selfEnergy[i][j]+templateEnergy[i][j]));
				eout << t;
			} else {
				fixedE += (selfEnergy[i][j]+templateEnergy[i][j]);
			}
		}

		if (selfEnergy[i].size() > 1){
			variableIndexI++;
		}
	}

	variableIndexI = 0;
	for (uint i = 0; i< pairEnergy.size();i++){
		if (selfEnergy[i].size() <= 1) continue;

		for (uint j = 0; j < pairEnergy[i].size();j++){
			int variableIndexJ = 0;
			for (uint ii = 0; ii < pairEnergy[i][j].size();ii++){
				if (selfEnergy[ii].size() <= 1) continue;
				for (uint jj = 0; jj < pairEnergy[i][j][ii].size();jj++){
					char t[100];
					sprintf(t, "%6d %6d %6d %6d %8.6f\n", variableIndexI,j,variableIndexJ,jj,pairEnergy[i][j][ii][jj]);
					eout << t;
				}

				variableIndexJ++;
			}
		}
		variableIndexI++;
	}
	char g[40];
	sprintf(g, "Fixed: %8.3f\n",fixedE);
	eout << g;
	eout.close();

	*/	

}
