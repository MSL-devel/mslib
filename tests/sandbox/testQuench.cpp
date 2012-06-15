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
#include <sstream>
#include <vector>
#include <ostream>
#include <fstream>
#include <cmath>



#include "testQuench.h"
#include "Quench.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"
#include "SystemRotamerLoader.h"
#include "PairwiseEnergyCalculator.h"
#include "AtomicPairwiseEnergy.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"

using namespace std;

using namespace MSL;

#include "SysEnv.h"
static SysEnv SYSENV;


int main(int argc, char *argv[]){

	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	Quench quencher(opt.topfile, opt.parfile, opt.rotlib);


	// Set the number of rotamers you want for large and small side chains, respectively
	quencher.setVariableNumberRotamers(50,5);


	//System sys = quencher.runQuench(initialSystem);



	// Set which positions you want to repack (optional)
	vector<int> variablePositions;
	variablePositions.push_back(0);
	variablePositions.push_back(1);
	variablePositions.push_back(2);
	variablePositions.push_back(3);
	variablePositions.push_back(4);
	// Could instead use: variablePositions.push_back(initialSystem.getPositionIndex("A", "73");

	System sys = quencher.runQuench(initialSystem, variablePositions);
	// For all positions, use: System sys = quencher.runQuench(initialSystem);

	stringstream ss;
	ss << "/tmp/currentConformation.pdb";
	string filename = ss.str();
	cout << "Write pdb " << filename << endl;
	PDBWriter writer;
	writer.open(filename);
	if (!writer.write(sys.getAtomPointers())) {
		cerr << "Problem writing " << filename << endl;
	}
	writer.close();

}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "testQuench --pdb PDB\n";
		exit(0);
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		cerr << "WARNING no topfile specified, using default "<<SYSENV.getEnv("MSL_CHARMM_TOP")<<endl;
		opt.topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		cerr << "WARNING no parfile specified, using default "<<SYSENV.getEnv("MSL_CHARMM_PAR")<<endl;
		opt.parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
	}

	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "WARNING no rotlib specified, using default "<<SYSENV.getEnv("MSL_ROTLIB")<<endl;
		opt.rotlib = SYSENV.getEnv("MSL_ROTLIB");
	}

	return opt;
}
