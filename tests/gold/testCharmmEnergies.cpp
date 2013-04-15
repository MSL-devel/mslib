/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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
#include "OnTheFlyManager.h"
#include "Transforms.h"
using namespace std;

using namespace MSL;

#include "SysEnv.h"
static SysEnv SYSENV;


vector<double> runTest(unsigned int _type, bool _writePdb);
bool validate(vector<double> & _t1, vector<double> & _t2, vector<double> & _t3, vector<double> & _t4, double _epsilon);


int main(int argc, char *argv[]) {

	bool writePDB_flag = false;
	if (argc > 1) {
		string argument = (string)argv[1];
		if (argument == "--writePdb") {
			writePDB_flag = true;
			cout << "Write PDB flag is on, will write PDB file to temp" << endl;
		}
	}

	ofstream out_fs;
	vector<double> allE;
	string outfile;

	vector<double> t1 = runTest(1, writePDB_flag);
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

	vector<double> t2 = runTest(2, writePDB_flag);
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

	vector<double> t3 = runTest(3, writePDB_flag);
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

	vector<double> t4 = runTest(4, writePDB_flag);
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
	orig1.push_back(5738.72964026776);
	orig1.push_back(4506.77714543232);
	orig1.push_back(4489.21437709076);
	orig1.push_back(4532.83527180681);
	orig1.push_back(4281.82749308834);
	orig1.push_back(4275.71469593009);
	orig1.push_back(4312.08304006534);
	orig1.push_back(5615.692504122);
	orig1.push_back(4383.74000890034);
	orig1.push_back(4366.17737344503);
	orig1.push_back(4409.79809997772);
	orig1.push_back(4158.79282165969);
	orig1.push_back(4152.67867529556);
	orig1.push_back(4189.04861547748);
	orig1.push_back(69057.2811548337);
	orig1.push_back(67825.3350131501);
	orig1.push_back(67807.7589940671);
	orig1.push_back(67851.495123945);
	orig1.push_back(67600.3391510038);
	orig1.push_back(67594.2048305189);
	orig1.push_back(67630.5843137652);
	orig1.push_back(68934.2439368011);
	orig1.push_back(67702.2977947313);
	orig1.push_back(67684.7219085346);
	orig1.push_back(67728.4578702292);
	orig1.push_back(67477.3043976884);
	orig1.push_back(67471.1687279976);
	orig1.push_back(67507.5498072905);
	orig1.push_back(5830.66394408087);
	orig1.push_back(4598.71005228592);
	orig1.push_back(4581.15360165196);
	orig1.push_back(4624.70452633226);
	orig1.push_back(4373.69609467866);
	orig1.push_back(4367.55705089062);
	orig1.push_back(4403.94020021866);
	orig1.push_back(5707.6268133217);
	orig1.push_back(4475.67292114053);
	orig1.push_back(4458.11660339281);
	orig1.push_back(4501.66735988976);
	orig1.push_back(4250.66142863659);
	orig1.push_back(4244.52103564269);
	orig1.push_back(4280.9057810174);
	orig1.push_back(72590.0677075907);
	orig1.push_back(71358.1214857456);
	orig1.push_back(71340.544615113);
	orig1.push_back(71384.2912090589);
	orig1.push_back(71133.1354437028);
	orig1.push_back(71126.9977570644);
	orig1.push_back(71163.381538374);
	orig1.push_back(72467.0304853565);
	orig1.push_back(71235.0842631252);
	orig1.push_back(71217.5075253789);
	orig1.push_back(71261.2539511414);
	orig1.push_back(71010.1006861857);
	orig1.push_back(71003.9616503414);
	orig1.push_back(71040.3470276977);

	vector<double> orig2;
	orig2.push_back(5738.72964026773);
	orig2.push_back(4506.77714543234);
	orig2.push_back(4489.21437709079);
	orig2.push_back(4532.83527180683);
	orig2.push_back(4281.82749308837);
	orig2.push_back(4275.71469593012);
	orig2.push_back(4312.08304006536);
	orig2.push_back(5615.69250412197);
	orig2.push_back(4383.74000890035);
	orig2.push_back(4366.17737344504);
	orig2.push_back(4409.79809997772);
	orig2.push_back(4158.7928216597);
	orig2.push_back(4152.67867529557);
	orig2.push_back(4189.04861547749);
	orig2.push_back(69057.2811548334);
	orig2.push_back(67825.3350131499);
	orig2.push_back(67807.7589940669);
	orig2.push_back(67851.4951239448);
	orig2.push_back(67600.3391510036);
	orig2.push_back(67594.2048305186);
	orig2.push_back(67630.5843137649);
	orig2.push_back(68934.2439368008);
	orig2.push_back(67702.2977947311);
	orig2.push_back(67684.7219085344);
	orig2.push_back(67728.4578702289);
	orig2.push_back(67477.3043976881);
	orig2.push_back(67471.1687279973);
	orig2.push_back(67507.5498072903);
	orig2.push_back(5830.66394408082);
	orig2.push_back(4598.71005228594);
	orig2.push_back(4581.15360165197);
	orig2.push_back(4624.70452633227);
	orig2.push_back(4373.69609467867);
	orig2.push_back(4367.55705089064);
	orig2.push_back(4403.94020021867);
	orig2.push_back(5707.62681332166);
	orig2.push_back(4475.67292114054);
	orig2.push_back(4458.11660339282);
	orig2.push_back(4501.66735988975);
	orig2.push_back(4250.66142863659);
	orig2.push_back(4244.52103564269);
	orig2.push_back(4280.9057810174);
	orig2.push_back(72590.0677075904);
	orig2.push_back(71358.1214857454);
	orig2.push_back(71340.5446151128);
	orig2.push_back(71384.2912090587);
	orig2.push_back(71133.1354437025);
	orig2.push_back(71126.9977570641);
	orig2.push_back(71163.3815383737);
	orig2.push_back(72467.0304853562);
	orig2.push_back(71235.0842631249);
	orig2.push_back(71217.5075253786);
	orig2.push_back(71261.2539511411);
	orig2.push_back(71010.1006861854);
	orig2.push_back(71003.9616503411);
	orig2.push_back(71040.3470276974);

	vector<double> orig3;
	orig3.push_back(5738.72964026776);
	orig3.push_back(4506.77714543238);
	orig3.push_back(4489.21437709082);
	orig3.push_back(4532.83527180686);
	orig3.push_back(4281.8274930884);
	orig3.push_back(4275.71469593015);
	orig3.push_back(4312.08304006539);
	orig3.push_back(5615.692504122);
	orig3.push_back(4383.74000890039);
	orig3.push_back(4366.17737344507);
	orig3.push_back(4409.79809997775);
	orig3.push_back(4158.79282165973);
	orig3.push_back(4152.67867529561);
	orig3.push_back(4189.04861547753);
	orig3.push_back(69057.2811548334);
	orig3.push_back(67825.3350131499);
	orig3.push_back(67807.7589940669);
	orig3.push_back(67851.4951239448);
	orig3.push_back(67600.3391510036);
	orig3.push_back(67594.2048305186);
	orig3.push_back(67630.5843137649);
	orig3.push_back(68934.2439368008);
	orig3.push_back(67702.2977947311);
	orig3.push_back(67684.7219085344);
	orig3.push_back(67728.4578702289);
	orig3.push_back(67477.3043976881);
	orig3.push_back(67471.1687279973);
	orig3.push_back(67507.5498072903);
	orig3.push_back(5830.66394408085);
	orig3.push_back(4598.71005228597);
	orig3.push_back(4581.15360165201);
	orig3.push_back(4624.7045263323);
	orig3.push_back(4373.6960946787);
	orig3.push_back(4367.55705089067);
	orig3.push_back(4403.94020021871);
	orig3.push_back(5707.62681332169);
	orig3.push_back(4475.67292114058);
	orig3.push_back(4458.11660339286);
	orig3.push_back(4501.66735988979);
	orig3.push_back(4250.66142863662);
	orig3.push_back(4244.52103564272);
	orig3.push_back(4280.90578101744);
	orig3.push_back(72590.0677075905);
	orig3.push_back(71358.1214857454);
	orig3.push_back(71340.5446151128);
	orig3.push_back(71384.2912090587);
	orig3.push_back(71133.1354437025);
	orig3.push_back(71126.9977570641);
	orig3.push_back(71163.3815383738);
	orig3.push_back(72467.0304853563);
	orig3.push_back(71235.0842631249);
	orig3.push_back(71217.5075253786);
	orig3.push_back(71261.2539511411);
	orig3.push_back(71010.1006861854);
	orig3.push_back(71003.9616503411);
	orig3.push_back(71040.3470276975);

	vector<double> orig4;
	orig4.push_back(5738.72964026776);
	orig4.push_back(4506.77714543238);
	orig4.push_back(4489.21437709082);
	orig4.push_back(4532.83527180686);
	orig4.push_back(4281.8274930884);
	orig4.push_back(4275.71469593015);
	orig4.push_back(4312.08304006539);
	orig4.push_back(5615.692504122);
	orig4.push_back(4383.74000890039);
	orig4.push_back(4366.17737344507);
	orig4.push_back(4409.79809997775);
	orig4.push_back(4158.79282165973);
	orig4.push_back(4152.67867529561);
	orig4.push_back(4189.04861547753);
	orig4.push_back(69057.2811548334);
	orig4.push_back(67825.3350131499);
	orig4.push_back(67807.7589940669);
	orig4.push_back(67851.4951239448);
	orig4.push_back(67600.3391510036);
	orig4.push_back(67594.2048305187);
	orig4.push_back(67630.5843137649);
	orig4.push_back(68934.2439368009);
	orig4.push_back(67702.2977947311);
	orig4.push_back(67684.7219085344);
	orig4.push_back(67728.4578702289);
	orig4.push_back(67477.3043976881);
	orig4.push_back(67471.1687279973);
	orig4.push_back(67507.5498072903);
	orig4.push_back(5830.66394408085);
	orig4.push_back(4598.71005228597);
	orig4.push_back(4581.15360165201);
	orig4.push_back(4624.7045263323);
	orig4.push_back(4373.6960946787);
	orig4.push_back(4367.55705089067);
	orig4.push_back(4403.94020021871);
	orig4.push_back(5707.62681332169);
	orig4.push_back(4475.67292114058);
	orig4.push_back(4458.11660339286);
	orig4.push_back(4501.66735988979);
	orig4.push_back(4250.66142863663);
	orig4.push_back(4244.52103564272);
	orig4.push_back(4280.90578101744);
	orig4.push_back(72590.0677075905);
	orig4.push_back(71358.1214857454);
	orig4.push_back(71340.5446151128);
	orig4.push_back(71384.2912090587);
	orig4.push_back(71133.1354437025);
	orig4.push_back(71126.9977570641);
	orig4.push_back(71163.3815383738);
	orig4.push_back(72467.0304853563);
	orig4.push_back(71235.0842631249);
	orig4.push_back(71217.5075253786);
	orig4.push_back(71261.2539511411);
	orig4.push_back(71010.1006861854);
	orig4.push_back(71003.9616503411);
	orig4.push_back(71040.3470276974);

	if (_t1.size() != orig1.size()) {
		return false;
	}

	bool out = true;
	for (unsigned int i=0; i<_t1.size(); i++) {
		if (fabs(_t1[i] - orig1[i]) < _epsilon) {
			cout << "1." << i << " OK" << endl;
		} else {
			cout << "1." << i << " NOT OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t2.size(); i++) {
		if (fabs(_t2[i] - orig2[i]) < _epsilon) {
			cout << "2." << i << " OK" << endl;
		} else {
			cout << "2." << i << " NOT OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t3.size(); i++) {
		if (fabs(_t3[i] - orig3[i]) < _epsilon) {
			cout << "3." << i << " OK" << endl;
		} else {
			cout << "3." << i << " NOT OK" << endl;
			out = false;
		}
	}
	for (unsigned int i=0; i<_t4.size(); i++) {
		if (fabs(_t4[i] - orig4[i]) < _epsilon) {
			cout << "4." << i << " OK" << endl;
		} else {
			cout << "4." << i << " NOT OK" << endl;
			out = false;
		}
	}
	return out;
}

