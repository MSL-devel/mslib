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
#include "OptionParser.h"
#include "energyOptimizations.h"


using namespace std;

using namespace MSL;

/*

  Process:
    *  Read-in input PDB

    *  Build   PolymerSequence.

    *  Build new system from PolymerSequence (this creates bonded energy terms too)

    *  Apply coordinates from structure from PDB
	
    *  Build rotamers

    *  Create energy table

 */

int main(int argc, char *argv[]){

	// Option Parser
	cout << "Setup Options"<<endl;
	EnergyTableOptions opt = setupEnergyTableOptions(argc,argv);

	// Create a system from the structural input options
	System sys;
	createSystem(opt.structOpt, sys);

	// create energy table specified by user...
	cout << "Calculating energy table..."<<endl;
	OnTheFlyManager otfm(&sys, opt.structOpt.parfile);
	CharmmEnergyCalculator * calculator = otfm.getCharmmEnergyCalculator();

	// Set parameters for energy calc
	calculator->setVdwRescalingFactor(opt.vdwScale);
	calculator->setDielectricConstant(opt.dielectric);
	calculator->setUseRdielectric(opt.distanceDependentElectrostatics);

	otfm.calculateEnergyTable(sys);

	// get pair table, print out.
	vector<vector<double> > selfEnergy = otfm.getSelfTable();
	vector<vector<double> > templateEnergy = otfm.getTemplateTable();
	vector<vector<vector<vector<double> > > > pairEnergy = otfm.getPairTable();
	
	ofstream eout;
	eout.open(opt.energyTableName.c_str());
	double fixedE = 0.0;
	int variableIndexI = 0;
	for (uint i = 0; i< selfEnergy.size();i++){
		if (selfEnergy[i].size() > 1){
			fprintf(stdout,"%1s %4d is index %6d\n",sys.getPosition(i).getChainId().c_str(),sys.getPosition(i).getResidueNumber(),variableIndexI);
		}
		for (uint j = 0; j < selfEnergy[i].size();j++){
			if (selfEnergy[i].size() > 1){
		
				char t[40];
				sprintf(t, "%6d %6d %8.3f\n", variableIndexI,j,(selfEnergy[i][j]+templateEnergy[i][j]));
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
					sprintf(t, "%6d %6d %6d %6d %8.3f\n", variableIndexI,j,variableIndexJ,jj,pairEnergy[i][j][ii][jj]);
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


	
	cout << "Done."<<endl;

}







