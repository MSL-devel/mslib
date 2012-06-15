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

// the backbone atoms to be read (PRO is given entirely)
std::string initialPdb = "\n\
ATOM      1  N   ALA A   1       2.143   1.328   0.000  1.00  0.00           N\n\
ATOM      2  CA  ALA A   1       1.539   0.000   0.000  1.00  0.00           C\n\
ATOM     11  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n\
ATOM     12  O   ALA A   1      -0.661   1.037   0.000  1.00  0.00           O\n\
ATOM     13  N   ARG A   2      -0.612  -1.210  -0.000  1.00  0.00           N\n\
ATOM     15  CA  ARG A   2      -2.052  -1.462  -0.000  1.00  0.00           C\n\
ATOM     35  C   ARG A   2      -2.314  -2.962  -0.000  1.00  0.00           C\n\
ATOM     36  O   ARG A   2      -1.393  -3.772  -0.000  1.00  0.00           O\n\
ATOM     37  N   ASN A   3      -3.605  -3.361  -0.000  1.00  0.00           N\n\
ATOM     39  CA  ASN A   3      -4.042  -4.750  -0.000  1.00  0.00           C\n\
ATOM     49  C   ASN A   3      -5.566  -4.690  -0.000  1.00  0.00           C\n\
ATOM     50  O   ASN A   3      -6.144  -3.607   0.000  1.00  0.00           O\n\
ATOM     51  N   ILE A   4      -6.232  -5.861  -0.000  1.00  0.00           N\n\
ATOM     53  CA  ILE A   4      -7.678  -5.994   0.000  1.00  0.00           C\n\
ATOM     68  C   ILE A   4      -7.971  -7.484  -0.000  1.00  0.00           C\n\
ATOM     69  O   ILE A   4      -7.052  -8.303  -0.000  1.00  0.00           O\n\
ATOM     70  N   CYS A   5      -9.259  -7.875  -0.000  1.00  0.00           N\n\
ATOM     72  CA  CYS A   5      -9.690  -9.256  -0.000  1.00  0.00           C\n\
ATOM     79  C   CYS A   5     -11.210  -9.218  -0.000  1.00  0.00           C\n\
ATOM     80  O   CYS A   5     -11.809  -8.143   0.000  1.00  0.00           O\n\
ATOM     81  N   GLU A   6     -11.880 -10.390  -0.000  1.00  0.00           N\n\
ATOM     83  CA  GLU A   6     -13.330 -10.510   0.000  1.00  0.00           C\n\
ATOM     94  C   GLU A   6     -13.660 -11.995  -0.000  1.00  0.00           C\n\
ATOM     95  O   GLU A   6     -12.768 -12.844  -0.000  1.00  0.00           O\n\
ATOM     96  N   GLN A   7     -14.965 -12.339  -0.000  1.00  0.00           N\n\
ATOM     98  CA  GLN A   7     -15.452 -13.708  -0.000  1.00  0.00           C\n\
ATOM    111  C   GLN A   7     -16.968 -13.628   0.000  1.00  0.00           C\n\
ATOM    112  OT1 GLN A   7     -17.500 -12.486   0.000  1.00  0.00           O\n\
ATOM    113  OT2 GLN A   7     -17.617 -14.708  -0.000  1.00  0.00           O\n\
TER     114      GLN A   7                                                    \n\
ATOM    115  N   GLY B   1      14.969   5.377   9.000  1.00  0.00           N\n\
ATOM    116  CA  GLY B   1      14.497   4.000   9.000  1.00  0.00           C\n\
ATOM    122  C   GLY B   1      13.000   4.000   9.000  1.00  0.00           C\n\
ATOM    123  O   GLY B   1      12.370   5.055   9.000  1.00  0.00           O\n\
ATOM    124  N   HSE B   2      12.376   2.805   9.000  1.00  0.00           N\n\
ATOM    126  CA  HSE B   2      10.929   2.641   9.000  1.00  0.00           C\n\
ATOM    139  C   HSE B   2      10.667   1.147   9.000  1.00  0.00           C\n\
ATOM    140  O   HSE B   2      11.599   0.347   9.000  1.00  0.00           O\n\
ATOM    141  N   ILE B   3       9.381   0.754   9.000  1.00  0.00           N\n\
ATOM    143  CA  ILE B   3       8.935  -0.626   9.000  1.00  0.00           C\n\
ATOM    158  C   ILE B   3       7.416  -0.584   9.000  1.00  0.00           C\n\
ATOM    159  O   ILE B   3       6.820   0.491   9.000  1.00  0.00           O\n\
ATOM    160  N   LEU B   4       6.752  -1.756   9.000  1.00  0.00           N\n\
ATOM    162  CA  LEU B   4       5.310  -1.873   9.000  1.00  0.00           C\n\
ATOM    177  C   LEU B   4       5.010  -3.362   9.000  1.00  0.00           C\n\
ATOM    178  O   LEU B   4       5.925  -4.184   9.000  1.00  0.00           O\n\
ATOM    179  N   LYS B   5       3.720  -3.745   9.000  1.00  0.00           N\n\
ATOM    181  CA  LYS B   5       3.279  -5.123   9.000  1.00  0.00           C\n\
ATOM    199  C   LYS B   5       1.761  -5.112   9.000  1.00  0.00           C\n\
ATOM    200  O   LYS B   5       1.140  -4.052   9.000  1.00  0.00           O\n\
ATOM    201  N   MET B   6       1.134  -6.305   9.000  1.00  0.00           N\n\
ATOM    203  CA  MET B   6      -0.303  -6.487   9.000  1.00  0.00           C\n\
ATOM    216  C   MET B   6      -0.543  -7.987   9.000  1.00  0.00           C\n\
ATOM    217  OT1 MET B   6       0.462  -8.747   9.000  1.00  0.00           O\n\
ATOM    218  OT2 MET B   6      -1.735  -8.395   9.000  1.00  0.00           O\n\
TER     219      MET B   6                                                    \n\
ATOM    220  N   MET C   1      -3.073  -8.607  -8.000  1.00  0.00           N\n\
ATOM    221  CA  MET C   1      -3.481 -10.000  -8.000  1.00  0.00           C\n\
ATOM    237  C   MET C   1      -5.000 -10.000  -8.000  1.00  0.00           C\n\
ATOM    238  O   MET C   1      -5.626  -8.943  -8.000  1.00  0.00           O\n\
ATOM    239  N   PHE C   2      -5.627 -11.192  -8.000  1.00  0.00           N\n\
ATOM    241  CA  PHE C   2      -7.064 -11.363  -8.000  1.00  0.00           C\n\
ATOM    257  C   PHE C   2      -7.318 -12.865  -8.000  1.00  0.00           C\n\
ATOM    258  O   PHE C   2      -6.378 -13.656  -8.000  1.00  0.00           O\n\
ATOM    259  N   PRO C   3      -8.600 -13.282  -8.000  1.00  0.00           N\n\
ATOM    260  CD  PRO C   3      -9.739 -12.366  -8.035  1.00  0.00           C\n\
ATOM    261  HD1 PRO C   3      -9.788 -11.799  -7.078  1.00  0.00           H\n\
ATOM    262  HD2 PRO C   3      -9.679 -11.666  -8.900  1.00  0.00           H\n\
ATOM    263  CA  PRO C   3      -9.001 -14.678  -8.000  1.00  0.00           C\n\
ATOM    264  HA  PRO C   3      -8.423 -15.247  -7.281  1.00  0.00           H\n\
ATOM    265  CB  PRO C   3     -10.486 -14.600  -7.598  1.00  0.00           C\n\
ATOM    266  HB1 PRO C   3     -10.559 -14.566  -6.487  1.00  0.00           H\n\
ATOM    267  HB2 PRO C   3     -11.077 -15.459  -7.972  1.00  0.00           H\n\
ATOM    268  CG  PRO C   3     -10.966 -13.264  -8.174  1.00  0.00           C\n\
ATOM    269  HG1 PRO C   3     -11.856 -12.863  -7.650  1.00  0.00           H\n\
ATOM    270  HG2 PRO C   3     -11.200 -13.394  -9.256  1.00  0.00           H\n\
ATOM    271  C   PRO C   3      -8.821 -15.301  -9.397  1.00  0.00           C\n\
ATOM    272  O   PRO C   3      -8.437 -14.612 -10.343  1.00  0.00           O\n\
ATOM    273  N   SER C   4      -9.174 -16.608  -9.486  1.00  0.00           N\n\
ATOM    275  CA  SER C   4      -9.094 -17.409 -10.694  1.00  0.00           C\n\
ATOM    296  C   THR C   5      -9.917 -21.795 -12.402  1.00  0.00           C\n\
ATOM    297  O   THR C   5      -9.489 -21.178 -13.376  1.00  0.00           O\n\
ATOM    298  N   TRP C   6     -10.270 -23.090 -12.478  1.00  0.00           N\n\
ATOM    300  CA  TRP C   6     -10.193 -23.899 -13.682  1.00  0.00           C\n\
ATOM    320  C   TRP C   6     -10.684 -25.295 -13.334  1.00  0.00           C\n\
ATOM    321  O   TRP C   6     -11.057 -25.571 -12.194  1.00  0.00           O\n\
ATOM    322  N   TYR C   7     -10.696 -26.217 -14.320  1.00  0.00           N\n\
ATOM    324  CA  TYR C   7     -11.132 -27.594 -14.165  1.00  0.00           C\n\
ATOM    341  C   TYR C   7     -10.972 -28.256 -15.527  1.00  0.00           C\n\
ATOM    342  O   TYR C   7     -10.542 -27.619 -16.486  1.00  0.00           O\n\
ATOM    343  N   VAL C   8     -11.320 -29.555 -15.629  1.00  0.00           N\n\
ATOM    345  CA  VAL C   8     -11.229 -30.332 -16.851  1.00  0.00           C\n\
ATOM    357  C   VAL C   8     -11.722 -31.721 -16.487  1.00  0.00           C\n\
ATOM    358  OT1 VAL C   8     -12.087 -31.926 -15.299  1.00  0.00           O\n\
ATOM    359  OT2 VAL C   8     -11.742 -32.598 -17.392  1.00  0.00           O\n\
TER     360      VAL C   8                                                    \n\
END\n";


int main() {

	/***************************************************
	 *  Test 1, 2 and 3:
	 *
	 *  A: [ALA GLY] ARG ASN [ILE(3) ASP(3)] CYS GLU GLN
	 *  B: GLY HSE ILE [LEU ALA LYS] LYS [MET SER]
	 *  C: MET PHE(3) PRO SER THR TRP [TYR THR] VAL
	 *
	 *  5 positions with double (or triple) identity (and
	 *  a position (C2) with multiple rotamers).
	 *
	 *  In test 1 all identities are created when the System
	 *  is built
	 *  In test 2 the additional identities are added with a
	 *  call to CharmmSystemBuilder::addIdentity
	 *  In test 3 all identities are created when the System
	 *  is built, but one is different and is replaced 
	 *
	 *  Energies are computed using the SelfPairManager and
	 *  with direct calculation.
	 *
	 *  Diff to verify output files 
	 *    /tmp/testAddCharmmIdentities-01.txt and
	 *    /tmp/testAddCharmmIdentities-02.txt and
	 *    /tmp/testAddCharmmIdentities-03.txt
	 *  
	 ***************************************************/
	test1();
	test2();
	test3();

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

	PDBReader pdbin;
	pdbin.read(initialPdb);
	pdbin.close();
	sys.assignCoordinates(pdbin.getAtomPointers(), false);

	//if (!sys.seed()) {
	//	cerr << "cannot seed the system" << endl;
	//}
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
	/*
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
	*/
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

	PDBReader pdbin;
	pdbin.read(initialPdb);
	pdbin.close();
	sys.assignCoordinates(pdbin.getAtomPointers(), false);

//	if (!sys.seed()) {
//		cerr << "cannot seed the system" << endl;
//	}
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

	//sys.buildAtoms(); // build the active atoms
	sys.buildAllAtoms(); // build the active atoms
	/*
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
	*/
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

void test3(){

	string baseNum = "03";

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
C: MET LEU PRO SER THR TRP [TYR THR] VAL");


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

	PDBReader pdbin;
	pdbin.read(initialPdb);
	pdbin.close();
	sys.assignCoordinates(pdbin.getAtomPointers(), false);

	//if (!sys.seed()) {
	//	cerr << "cannot seed the system" << endl;
	//}
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
	/*
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
	*/
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

	CSB.addIdentity("C,2", "PHE");
	CSB.removeIdentity("C,2", "LEU");
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

