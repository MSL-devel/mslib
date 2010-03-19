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





#include "OptionParser.h"
#include "System.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "PairwiseEnergyCalculator.h"

#include "testData.h"

using namespace MSL;

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

int main() {

	writePdbFile();


	System initialSystem;
	initialSystem.readPdb("/tmp/symmetricTrimer.pdb");

	stringstream seq;
	for (uint c = 0; c< initialSystem.size();c++){
		Chain ch = initialSystem.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<" "<<ch.getPosition(0).getResidueNumber()<<": ";

		for (uint p = 0 ; p < ch.size();p++){
			Position pos = ch.getPosition(p);

			cout << "Position: "<<pos.getResidueNumber()<<endl;
			string chainId = pos.getChainId();
			int resNum     = pos.getResidueNumber();

			int index = -1;

			if (resNum == 9 || resNum == 6){
				index = 1;
			}



			if (index != -1){
				seq << " [ VAL THR ]";
			} else {

				seq << " "<<pos.getCurrentIdentity().getResidueName();
			}

		}
		seq << "\n";
		
	}


	PolymerSequence pseq(seq.str());
	

	System sys;
	string topFile = "/library/charmmTopPar/top_all22_prot.inp";
	string parFile = "/library/charmmTopPar/par_all22_prot.inp";
	cout << "Use toppar " << topFile << ", " << parFile << endl;
	CharmmSystemBuilder CSB(topFile,parFile);

	// Check for type of energy calculation...
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(sys,pseq);


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


	SystemRotamerLoader sysRot(sys, "/library/rotlib/balanced/rotlib-balanced-200.txt");
		
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

		

	return 1;
}
