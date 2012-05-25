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


vector<double> runTest(unsigned int _type, bool _writePdb);
bool validate(vector<double> & _t1, vector<double> & _t2, vector<double> & _t3, vector<double> & _t4, double _epsilon);


int main() {
	ofstream out_fs;
	vector<double> allE;
	string outfile;

	vector<double> t1 = runTest(1, false);
	outfile = "/tmp/testCharmmEnergies-energies-1.txt";
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}
	for (unsigned int i=0; i<t1.size(); i++) {
		out_fs <<  setprecision(15) << t1[i] << endl;
	}
	out_fs.close();

	vector<double> t2 = runTest(2, false);
	outfile = "/tmp/testCharmmEnergies-energies-2.txt";
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}
	for (unsigned int i=0; i<t1.size(); i++) {
		out_fs <<  setprecision(15) << t2[i] << endl;
	}
	out_fs.close();

	vector<double> t3 = runTest(3, false);
	outfile = "/tmp/testCharmmEnergies-energies-3.txt";
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}
	for (unsigned int i=0; i<t1.size(); i++) {
		out_fs <<  setprecision(15) << t3[i] << endl;
	}
	out_fs.close();

	vector<double> t4 = runTest(4, false);
	outfile = "/tmp/testCharmmEnergies-energies-4.txt";
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}
	for (unsigned int i=0; i<t1.size(); i++) {
		out_fs <<  setprecision(15) << t4[i] << endl;
	}
	out_fs.close();

	cout << "============================================" << endl;
	cout << "Check the results" << endl;
	cout << endl;
	if (validate(t1, t2, t3, t4, 1e-8)) {
		cout << "GOLD" << endl;
	} else {
		cout << "LEAD" << endl;
	}

	return 0;
}


vector<double> runTest(unsigned int _type, bool _writePdb) {


	string outfile = "/tmp/testCharmmEnergies-output-" + MslTools::intToString(_type) + ".txt";
	ofstream out_fs;
	out_fs.open(outfile.c_str());
	if (out_fs.fail()) {
		cerr << "Error writing test output " << outfile << endl;
		exit(1);
	}

	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	cout << " ***********************     TEST " << _type << " ****************************"<<endl;
	out_fs << " ***********************     TEST " << _type << " ****************************"<<endl;

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
	out_fs << "Sequence:" << endl;
	out_fs << seq.toString();
	out_fs << endl;
	string topFile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parFile = SYSENV.getEnv("MSL_CHARMM_PAR");
	out_fs << "Use toppar " << topFile << ", " << parFile << endl;

	CharmmSystemBuilder CSB(sys, topFile, parFile);
	CSB.buildSystem(seq);
	out_fs << endl;

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
	out_fs << "Select and translate chain B by (13, 4, 9)" << endl;
	out_fs << "Selection chain B has " << as.size("chain B") << "atoms" << endl;
	tr.translate(chainB, CartesianPoint(13, 4, 9));
	
	AtomPointerVector chainC = as.select("chain C");
	out_fs << "Select and translate chain C by (-5, -10, -8)" << endl;
	out_fs << "Selection chain C has " << as.size("chain C") << "atoms" << endl;
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
	out_fs << endl;
	sys.buildAllAtoms();
	out_fs << endl;

	string filename = "/tmp/initialBuild.pdb";
	out_fs << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();



	/************************************************
	 *   S T A R T : ADD ROTAMERS
	 ************************************************/
	out_fs << "Add rotamers to the system" << endl;
	string rotlib = SYSENV.getEnv("MSL_ROTLIB");
	out_fs << "Read rotamer library " << rotlib << endl;
	out_fs << endl;
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

	vector<double> out;
	if (_type == 1) {

		cout << "ENERGIES CALCULATED WITH THE SYSTEM/ENERGY SET" << endl;
		cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;

		out_fs << "ENERGIES CALCULATED WITH THE SYSTEM/ENERGY SET" << endl;
		out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		out_fs << endl;

		vector<unsigned int> rots;
		rots.push_back(7);
		rots.push_back(2);
		rots.push_back(4);
		Enumerator stateEnum(rots);

		for (unsigned int i=0; i<stateEnum.size(); i++) {
			// for (unsigned int i=0; i<2; i++) {
			vector<unsigned int> state = stateEnum[i];
			out_fs << endl;
			out_fs << "STATE:";
			for (unsigned int j=0; j<state.size(); j++) {
				out_fs << " " << state[j];
			}
			out_fs << endl;
			sys.setActiveRotamers(state);
			out_fs << endl;
			double E = sys.calcEnergy();
			out.push_back(E);
			out_fs << "Direct calculation of energy (serial " << i << ") = " << E << endl;
			out_fs << sys.getEnergySummary();

			if (_writePdb) {
				char c[1000];
				sprintf(c, "/tmp/serial-%04u.pdb", i);
				sys.writePdb(c);
				out_fs << "Written pdb " << c << endl;
				out_fs << " - - - - - - - " << endl;
				out_fs << endl;
			}
		}
	} else if (_type == 2) {
		cout << "ENERGIES CALCULATED WITH THE SELF-PAIR MANAGER" << endl;
		cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;

		out_fs << "ENERGIES CALCULATED WITH THE SELF-PAIR MANAGER" << endl;
		out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		out_fs << endl;

		SelfPairManager SPM(&sys);
		SPM.saveEnergiesByTerm(true);
		SPM.saveInteractionCounts(true);
		SPM.calculateEnergies();

		vector<unsigned int> rots = SPM.getNumberOfRotamers();
		out_fs << "The system has " << rots.size() << " variable positions" << endl;
		for (unsigned int i=0; i<rots.size(); i++) {
			out_fs << "   Variable position " << i << " has " << rots[i] << " rotamers" << endl;
		}

		Enumerator stateEnum(rots);

		for (unsigned int i=0; i<stateEnum.size(); i++) {
			// for (unsigned int i=0; i<2; i++) {
			vector<unsigned int> state = stateEnum[i];
			out_fs << endl;
			out_fs << "STATE:";
			for (unsigned int j=0; j<state.size(); j++) {
				out_fs << " " << state[j];
			}
			out_fs << endl;
			double E = SPM.getStateEnergy(state);
			out.push_back(E);
			out_fs << "Energy = " << E << endl;
			vector<string> stateDescriptors = SPM.getStateDescriptors(state);
			vector<vector<unsigned int> > stateIndeces = SPM.getStatePositionIdentityRotamerIndeces(state);
			for (unsigned int j=0; j<stateDescriptors.size(); j++) {
				out_fs << " * " << stateDescriptors[j] << " -- " << stateIndeces[j][0] << "/" << stateIndeces[j][1] << "/" << stateIndeces[j][2] << endl;
			}

			out_fs << SPM.getSummary(state);

		}
	} else if (_type == 3) {

		cout << "ENERGIES WITH THE PAIRWISE ENERGY CALCULATOR: ON THE FLY" << endl;
		cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;

		out_fs << "ENERGIES WITH THE PAIRWISE ENERGY CALCULATOR: ON THE FLY" << endl;
		out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		out_fs << endl;

		OnTheFlyManager pec(&sys, SYSENV.getEnv("MSL_CHARMM_PAR"));

		vector<unsigned int> rots;
		rots.push_back(7);
		rots.push_back(2);
		rots.push_back(4);
		Enumerator stateEnum(rots);

		for (unsigned int i=0; i<stateEnum.size(); i++) {
			// for (unsigned int i=0; i<2; i++) {
			vector<unsigned int> state = stateEnum[i];
			out_fs << endl;
			out_fs << "STATE:";
			for (unsigned int j=0; j<state.size(); j++) {
				out_fs << " " << state[j];
			}
			out_fs << endl;
		
			double E = pec.calculateStateEnergy(sys,state);
			out.push_back(E);
			out_fs << "Energy = " << E << endl;
			//pec.printSummary();
			out_fs << pec.getSummary();

		}
	} else if (_type == 4) {

		cout << "ENERGIES WITH THE PAIRWISE ENERGY CALCULATOR: PRECOMPUTED TABLES" << endl;
		cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;

		out_fs << "ENERGIES WITH THE PAIRWISE ENERGY CALCULATOR: PRECOMPUTED TABLES" << endl;
		out_fs << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
		out_fs << endl;

		OnTheFlyManager pec(&sys, SYSENV.getEnv("MSL_CHARMM_PAR"));

		// Calculate the complete energy table.
		pec.calculateEnergyTable(sys);

		vector<unsigned int> rots;
		rots.push_back(7);
		rots.push_back(2);
		rots.push_back(4);
		Enumerator stateEnum(rots);

		for (unsigned int i=0; i<stateEnum.size(); i++) {
			// for (unsigned int i=0; i<2; i++) {
			vector<unsigned int> state = stateEnum[i];
			out_fs << endl;
			out_fs << "STATE:";
			for (unsigned int j=0; j<state.size(); j++) {
				out_fs << " " << state[j];
			}
			out_fs << endl;
		
			double E = pec.getStateEnergy(sys, state);
			out.push_back(E);
			out_fs << "Energy = " << E << endl;

		}
	}
	out_fs.close();
	cout << "DONE: output witten to " << outfile << endl;
	return out;
}


bool validate(vector<double> & _t1, vector<double> & _t2, vector<double> & _t3, vector<double> & _t4, double _epsilon) {

	vector<double> orig1;
	orig1.push_back(1909.60827668261);
	orig1.push_back(1193.53081845169);
	orig1.push_back(1096.08455442612);
	orig1.push_back(1265.30118768309);
	orig1.push_back(1103.15698559712);
	orig1.push_back(976.635403576887);
	orig1.push_back(979.460701678256);
	orig1.push_back(1760.33597509929);
	orig1.push_back(1044.25865347544);
	orig1.push_back(946.812228929392);
	orig1.push_back(1116.02888724301);
	orig1.push_back(953.887430794695);
	orig1.push_back(827.364393629392);
	orig1.push_back(830.191142895793);
	orig1.push_back(32625.6638418553);
	orig1.push_back(31909.5722416549);
	orig1.push_back(31812.280759274);
	orig1.push_back(31981.3629158052);
	orig1.push_back(31819.1917540136);
	orig1.push_back(31692.674949865);
	orig1.push_back(31695.4813793109);
	orig1.push_back(32476.3914457988);
	orig1.push_back(31760.2999822055);
	orig1.push_back(31663.0083393041);
	orig1.push_back(31832.0905208919);
	orig1.push_back(31669.922104738);
	orig1.push_back(31543.4038454443);
	orig1.push_back(31546.2117260553);
	orig1.push_back(1774.56058981724);
	orig1.push_back(1058.46593237396);
	orig1.push_back(961.20089636412);
	orig1.push_back(1130.25928773802);
	orig1.push_back(968.111503728769);
	orig1.push_back(841.626024016152);
	orig1.push_back(844.398395579353);
	orig1.push_back(1625.28819239328);
	orig1.push_back(909.19367155708);
	orig1.push_back(811.928475026759);
	orig1.push_back(980.986891457309);
	orig1.push_back(818.841853085706);
	orig1.push_back(692.354918228022);
	orig1.push_back(695.128740956255);
	orig1.push_back(2432.52561977074);
	orig1.push_back(1716.45284640629);
	orig1.push_back(1618.9666947077);
	orig1.push_back(1788.21893065533);
	orig1.push_back(1625.97115072215);
	orig1.push_back(1499.48302478068);
	orig1.push_back(1502.26039833597);
	orig1.push_back(2283.25330897927);
	orig1.push_back(1567.1806722219);
	orig1.push_back(1469.69436000282);
	orig1.push_back(1638.9466210071);
	orig1.push_back(1476.70158671158);
	orig1.push_back(1350.21200562504);
	orig1.push_back(1352.99083034536);

	vector<double> orig2;
	orig2.push_back(1909.60827668258);
	orig2.push_back(1193.53081845169);
	orig2.push_back(1096.08455442612);
	orig2.push_back(1265.30118768309);
	orig2.push_back(1103.15698559713);
	orig2.push_back(976.635403576885);
	orig2.push_back(979.460701678255);
	orig2.push_back(1760.33597509925);
	orig2.push_back(1044.25865347544);
	orig2.push_back(946.812228929388);
	orig2.push_back(1116.02888724302);
	orig2.push_back(953.887430794696);
	orig2.push_back(827.36439362939);
	orig2.push_back(830.191142895791);
	orig2.push_back(32625.6638418553);
	orig2.push_back(31909.5722416549);
	orig2.push_back(31812.280759274);
	orig2.push_back(31981.3629158052);
	orig2.push_back(31819.1917540136);
	orig2.push_back(31692.674949865);
	orig2.push_back(31695.481379311);
	orig2.push_back(32476.3914457989);
	orig2.push_back(31760.2999822055);
	orig2.push_back(31663.0083393041);
	orig2.push_back(31832.090520892);
	orig2.push_back(31669.9221047381);
	orig2.push_back(31543.4038454443);
	orig2.push_back(31546.2117260554);
	orig2.push_back(1774.56058981721);
	orig2.push_back(1058.46593237396);
	orig2.push_back(961.200896364116);
	orig2.push_back(1130.25928773802);
	orig2.push_back(968.111503728771);
	orig2.push_back(841.626024016148);
	orig2.push_back(844.39839557935);
	orig2.push_back(1625.28819239325);
	orig2.push_back(909.19367155708);
	orig2.push_back(811.92847502675);
	orig2.push_back(980.986891457309);
	orig2.push_back(818.841853085705);
	orig2.push_back(692.354918228016);
	orig2.push_back(695.128740956249);
	orig2.push_back(2432.52561977072);
	orig2.push_back(1716.45284640629);
	orig2.push_back(1618.96669470769);
	orig2.push_back(1788.21893065532);
	orig2.push_back(1625.97115072215);
	orig2.push_back(1499.48302478067);
	orig2.push_back(1502.26039833596);
	orig2.push_back(2283.25330897925);
	orig2.push_back(1567.18067222189);
	orig2.push_back(1469.69436000281);
	orig2.push_back(1638.9466210071);
	orig2.push_back(1476.70158671157);
	orig2.push_back(1350.21200562503);
	orig2.push_back(1352.99083034535);

	vector<double> orig3;
	orig3.push_back(1909.60827668258);
	orig3.push_back(1193.53081845169);
	orig3.push_back(1096.08455442612);
	orig3.push_back(1265.30118768309);
	orig3.push_back(1103.15698559712);
	orig3.push_back(976.635403576884);
	orig3.push_back(979.460701678253);
	orig3.push_back(1760.33597509926);
	orig3.push_back(1044.25865347544);
	orig3.push_back(946.812228929387);
	orig3.push_back(1116.02888724302);
	orig3.push_back(953.887430794695);
	orig3.push_back(827.364393629389);
	orig3.push_back(830.19114289579);
	orig3.push_back(32625.6638418553);
	orig3.push_back(31909.5722416549);
	orig3.push_back(31812.280759274);
	orig3.push_back(31981.3629158052);
	orig3.push_back(31819.1917540136);
	orig3.push_back(31692.674949865);
	orig3.push_back(31695.4813793109);
	orig3.push_back(32476.3914457988);
	orig3.push_back(31760.2999822055);
	orig3.push_back(31663.0083393041);
	orig3.push_back(31832.090520892);
	orig3.push_back(31669.922104738);
	orig3.push_back(31543.4038454443);
	orig3.push_back(31546.2117260553);
	orig3.push_back(1774.56058981721);
	orig3.push_back(1058.46593237396);
	orig3.push_back(961.200896364114);
	orig3.push_back(1130.25928773802);
	orig3.push_back(968.111503728769);
	orig3.push_back(841.626024016147);
	orig3.push_back(844.398395579348);
	orig3.push_back(1625.28819239325);
	orig3.push_back(909.19367155708);
	orig3.push_back(811.928475026748);
	orig3.push_back(980.98689145731);
	orig3.push_back(818.841853085704);
	orig3.push_back(692.354918228016);
	orig3.push_back(695.128740956248);
	orig3.push_back(2432.52561977071);
	orig3.push_back(1716.45284640629);
	orig3.push_back(1618.96669470768);
	orig3.push_back(1788.21893065532);
	orig3.push_back(1625.97115072215);
	orig3.push_back(1499.48302478067);
	orig3.push_back(1502.26039833596);
	orig3.push_back(2283.25330897925);
	orig3.push_back(1567.18067222189);
	orig3.push_back(1469.69436000281);
	orig3.push_back(1638.9466210071);
	orig3.push_back(1476.70158671157);
	orig3.push_back(1350.21200562503);
	orig3.push_back(1352.99083034535);

	vector<double> orig4;
	orig4.push_back(1909.60827668258);
	orig4.push_back(1193.53081845169);
	orig4.push_back(1096.08455442612);
	orig4.push_back(1265.30118768309);
	orig4.push_back(1103.15698559712);
	orig4.push_back(976.635403576884);
	orig4.push_back(979.460701678253);
	orig4.push_back(1760.33597509926);
	orig4.push_back(1044.25865347544);
	orig4.push_back(946.812228929387);
	orig4.push_back(1116.02888724302);
	orig4.push_back(953.887430794695);
	orig4.push_back(827.364393629389);
	orig4.push_back(830.19114289579);
	orig4.push_back(32625.6638418553);
	orig4.push_back(31909.5722416549);
	orig4.push_back(31812.280759274);
	orig4.push_back(31981.3629158052);
	orig4.push_back(31819.1917540136);
	orig4.push_back(31692.674949865);
	orig4.push_back(31695.481379311);
	orig4.push_back(32476.3914457989);
	orig4.push_back(31760.2999822055);
	orig4.push_back(31663.0083393041);
	orig4.push_back(31832.090520892);
	orig4.push_back(31669.922104738);
	orig4.push_back(31543.4038454443);
	orig4.push_back(31546.2117260553);
	orig4.push_back(1774.56058981721);
	orig4.push_back(1058.46593237396);
	orig4.push_back(961.200896364114);
	orig4.push_back(1130.25928773802);
	orig4.push_back(968.111503728769);
	orig4.push_back(841.626024016147);
	orig4.push_back(844.398395579348);
	orig4.push_back(1625.28819239325);
	orig4.push_back(909.193671557079);
	orig4.push_back(811.928475026748);
	orig4.push_back(980.986891457309);
	orig4.push_back(818.841853085704);
	orig4.push_back(692.354918228016);
	orig4.push_back(695.128740956248);
	orig4.push_back(2432.52561977071);
	orig4.push_back(1716.45284640629);
	orig4.push_back(1618.96669470768);
	orig4.push_back(1788.21893065532);
	orig4.push_back(1625.97115072215);
	orig4.push_back(1499.48302478067);
	orig4.push_back(1502.26039833596);
	orig4.push_back(2283.25330897925);
	orig4.push_back(1567.18067222189);
	orig4.push_back(1469.69436000281);
	orig4.push_back(1638.9466210071);
	orig4.push_back(1476.70158671157);
	orig4.push_back(1350.21200562503);
	orig4.push_back(1352.99083034535);

	if (_t1.size() != orig1.size()) {
		return false;
	}

	bool out = true;
	for (unsigned int i=0; i<_t1.size(); i++) {
		if (abs(_t1[i] - orig1[i]) < _epsilon) {
			cout << "1." << i << " OK" << endl;
		} else {
			cout << "1." << i << " OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t2.size(); i++) {
		if (abs(_t2[i] - orig2[i]) < _epsilon) {
			cout << "2." << i << " OK" << endl;
		} else {
			cout << "2." << i << " OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t3.size(); i++) {
		if (abs(_t3[i] - orig3[i]) < _epsilon) {
			cout << "3." << i << " OK" << endl;
		} else {
			cout << "3." << i << " OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t4.size(); i++) {
		if (abs(_t4[i] - orig4[i]) < _epsilon) {
			cout << "4." << i << " OK" << endl;
		} else {
			cout << "4." << i << " OK" << endl;
			out = false;
		}
	}
	return out;
}

