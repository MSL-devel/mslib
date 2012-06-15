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

#include "LinearProgrammingOptimization.h"
#include "MonteCarloOptimization.h"
#include "DeadEndElimination.h"
#include "testData.h"
using namespace std;

using namespace MSL;


int main() {

	// Write a test energy table from testData.h (/tmp/testEnergy.txt)
	writeEnergyTable(energyTable,"/tmp/testEnergy.txt");
	
	
	cout << "************* TEST DEAD END ELIMINATION ON SMALL ENERGY TABLE *************"<<endl;
	DeadEndElimination dee;
	dee.readEnergyTable("/tmp/testEnergy.txt");
	dee.runSimpleGoldsteinSingles();
	dee.printMe();
	cout << "Eliminated Rotamers: "<<dee.getEliminatedCounter()<<endl;


	cout << "************* TEST LINEAR PROGRAMMING OPTIMIZATION ON SMALL ENERGY TABLE *************"<<endl;
	LinearProgrammingOptimization lp;

	lp.readEnergyTable("/tmp/testEnergy.txt");
	lp.analyzeEnergyTable();	
	lp.createLP();
	lp.solveLP();

	lp.printMe();
	cout << "Optimized Energy: "<<lp.getTotalEnergy()<<endl;


	cout << "************* TEST LINEAR PROGRAMMING OPTIMIZATION ON SMALL ENERGY TABLE WITH INPUT MASK *************"<<endl;

	vector<vector<bool> > inputMask;
	inputMask.resize(lp.getNumPositions());
	for (uint i = 0; i < inputMask.size();i++){
		inputMask[i].resize(lp.getNumRotamers(i));
		for (uint j = 0; j < inputMask[i].size();j++){
			inputMask[i][j] = true;
		}
	}

	// Put some of the known global minimum rotamers as false.
	inputMask[0][5] = false;


	LinearProgrammingOptimization lp2;
	lp2.setVerbose(true);
	lp2.readEnergyTable("/tmp/testEnergy.txt");

	lp2.setInputRotamerMasks(inputMask);
	//lp2.preFilterEnergyTable();

	lp2.analyzeEnergyTable();	
	lp2.createLP();
	lp2.solveLP();

	lp2.printMe();
	cout << "Optimized Energy: "<<lp2.getTotalEnergy()<<endl;


	cout << "************* TEST MONTE CARLO OPTIMIZATION ON SMALL ENERGY TABLE *************"<<endl;	

	MonteCarloOptimization mc;
	cout << "Read Energy Table"<<endl;
	mc.readEnergyTable("/tmp/testEnergy.txt");
	cout << "RunMC"<<endl;
//	mc.setAnnealSchedule(MonteCarloOptimization::EXP_TEMP_ANNEAL,100,1);
//	mc.setNumberOfCycles(10000);
	mc.setInitializationState(MonteCarloOptimization::LOWESTSELF);
	mc.runMC(100,1,10000,EXPONENTIAL,1000,1000,0.01);
	mc.printMe();
	//mc.printLowest() // IMPLEMENT THIS!
	cout << "Optimized Energy: "<<mc.getStateEnergy()<<endl;
	cout << "Sampled Energies: "<<endl;
	mc.printSampledConfigurations();


	cout << "************* TEST MONTE CARLO OPTIMIZATION ON SMALL ENERGY TABLE WITH INPUT MASKS *************"<<endl;	

	// Monte Carlo with input masks
	MonteCarloOptimization mc2;
	cout << "Read Energy Table"<<endl;
	mc2.readEnergyTable("/tmp/testEnergy.txt");

	inputMask.clear();
	inputMask.resize(mc2.getNumPositions());
	for (uint i = 0; i < inputMask.size();i++){
		inputMask[i].resize(mc2.getNumRotamers(i));
		for (uint j = 0; j < inputMask[i].size();j++){
			inputMask[i][j] = true;
		}
	}

	// Put some of the known global minimum rotamers as false.
	inputMask[0][5] = false;


	
	mc2.setInputRotamerMasks(inputMask);

	//mc2.setAnnealSchedule(MonteCarloOptimization::EXP_TEMP_ANNEAL,100,1);
	//mc2.setNumberOfCycles(10000);
	mc2.setInitializationState(MonteCarloOptimization::LOWESTSELF);

	cout << "RunMC2"<<endl;
	mc2.runMC(100,1,10000,EXPONENTIAL,1000,1000,0.01);

	// Print Results
	mc2.printMe();
	//mc2.printLowest() // IMPLEMENT THIS!
	cout << "Optimized Energy: "<<mc2.getStateEnergy()<<endl;
	cout << "Sampled Energies: "<<endl;
	mc2.printSampledConfigurations();

	


	/*
	  TEST DEE ON LARGE ENERGY TABLE...
	 */


	cout << "************* TEST DEADEND ELIMINIATION OPTIMIZATION ON LARGE ENERGY TABLE *************"<<endl;	
	string bigEneTable;
	//bigEneTable = "/data/dwkulp/DFIRE/DFIRE/bigRots/out-MMFF/oligomerSelfPair.tbl.bin";
	bigEneTable = "/home/dwkulp/oligomerSelfPair.tbl.bin";


	
	///data/dwkulp/DFIRE/DFIRE/bigRots/out-MMFF/
	DeadEndElimination deeBig;
	deeBig.readEnergyTable(bigEneTable);
	deeBig.runSimpleGoldsteinSingles();
	deeBig.printMe();
	cout << "Eliminated Rotamers: "<<deeBig.getEliminatedCounter()<<" out of "<<deeBig.getTotalNumberRotamers()<<endl;


	inputMask = deeBig.getMask();

	
	cout << "************* TEST MONTE CARLO OPTIMIZATION ON LARGE ENERGY TABLE AFTER DEADEND ELIMINATION  *************"<<endl;	
	MonteCarloOptimization mc3;
	cout << "Read Energy Table"<<endl;
	mc3.readEnergyTable(bigEneTable);
	mc3.setInputRotamerMasks(inputMask);

	//mc3.setAnnealSchedule(MonteCarloOptimization::EXP_TEMP_ANNEAL,100,1);
	//mc3.setNumberOfCycles(10000);
	mc3.setInitializationState(MonteCarloOptimization::LOWESTSELF);

	cout << "RunMC3"<<endl;
	mc3.runMC(100,1,10000,EXPONENTIAL,1000,1000,0.01);

	// Print Results
	mc3.printMe();
	//mc2.printLowest() // IMPLEMENT THIS!
	cout << "Optimized Energy: "<<mc3.getStateEnergy()<<endl;
	cout << "Sampled Energies: "<<endl;
	mc3.printSampledConfigurations();



	cout << "************* TEST LINEAR PROGRAMMING OPTIMIZATION ON LARGE ENERGY TABLE AFTER DEADEND ELIMINATION  *************"<<endl;	
	cout << "RunLP3"<<endl;
	LinearProgrammingOptimization lp3;
	lp3.readEnergyTable(bigEneTable);

	lp3.setInputRotamerMasks(inputMask);
	lp3.analyzeEnergyTable();	
	lp3.createLP();
	lp3.solveLP();

	lp3.printMe();
	cout << "Optimized Energy: "<<lp3.getTotalEnergy()<<endl;



	/*
	cout << "************* TEST LINEAR PROGRAMMING OPTIMIZATION ON SMALL ENERGY TABLE WITH POSITIONS THAT ARE LINKED  *************"<<endl;	
	LinearProgrammingOptimization lp4;
	lp4.readEnergyTable("/tmp/testEnergy.txt");
	lp4.analyzeEnergyTable();	
	lp4.createLP();
	lp4.solveLP();

	lp4.printMe();
	cout << "Optimized Energy: "<<lp4.getTotalEnergy()<<endl;
	*/



	cout << "************* TEST MONTE CARLO OPTIMIZATION ON SMALL ENERGY TABLE WITH POSITIONS THAT ARE LINKED *************"<<endl;	
	MonteCarloOptimization mc4;
	cout << "Read Energy Table"<<endl;
	mc4.readEnergyTable("/tmp/testEnergy.txt");


	//mc4.setAnnealSchedule(MonteCarloOptimization::EXP_TEMP_ANNEAL,100,1);
	//mc4.setNumberOfCycles(10000);
	mc4.setInitializationState(MonteCarloOptimization::LOWESTSELF);

	cout << "RunMC4"<<endl;
	mc4.runMC(100,1,10000,EXPONENTIAL,1000,1000,0.01);

	// Print Results
	mc4.printMe();
	//mc2.printLowest() // IMPLEMENT THIS!
	cout << "Optimized Energy: "<<mc4.getStateEnergy()<<endl;
	cout << "Sampled Energies: "<<endl;
	mc4.printSampledConfigurations();
	
}
