/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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



#include "runQuench.h"
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


int main(int argc, char *argv[]){

	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	Quench quencher(opt.topfile, opt.parfile, opt.rotlib);

	// Set the number of rotamers you want for large and small side chains, respectively
	quencher.setVariableNumberRotamers(opt.largeRotNum,opt.smallRotNum);

	// Set which positions you want to repack (optional)
	if (opt.positions.size() > 0) {	
		vector<int> variablePositions;
		for (uint i = 0; i < opt.positions.size(); i++) {
			vector<string> pos = MslTools::tokenize(opt.positions[i],"_");
			variablePositions.push_back(initialSystem.getPositionIndex(pos[0], pos[1]));
		}
		System sys = quencher.runQuench(initialSystem,variablePositions);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtoms())) {
			cerr << "Problem writing " << opt.outfile << endl;
		}
		writer.close();
	}
	else {
		System sys = quencher.runQuench(initialSystem);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtoms())) {
			cerr << "Problem writing " << opt.outfile << endl;
		}
		writer.close();
	}
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
		cout << "runQuench --pdb PDB [--topfile TOPFILE --parfile PARFILE --rotlib ROTLIB --outfile OUTPDB --positions A_73 A_74 A_75 --smallRotNum NUM --largeRotNum NUM]\n";
		exit(0);
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		cerr << "WARNING no topfile specified, using default /library/charmmTopPar/top_all22_prot.inp\n";
		opt.topfile = "/library/charmmTopPar/top_all22_prot.inp";
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		cerr << "WARNING no parfile specified, using default /library/charmmTopPar/par_all22_prot.inp\n";
		opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
	}

	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "WARNING no rotlib specified, using default /library/rotlib/balanced/rotlib-balanced-200.txt\n";
		opt.rotlib = "/library/rotlib/balanced/rotlib-balanced-200.txt";
	}

	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		cerr << "WARNING no outfile specified, using default /tmp/currentConformation.pdb\n";
		opt.outfile = "/tmp/currentConformation.pdb";
	}

	opt.positions = OP.getStringVectorJoinAll("positions");

	opt.smallRotNum = OP.getInt("largeRotNum");
	if (OP.fail()){
		cerr << "WARNING largeRotNum not specified, using default 50\n";
		opt.largeRotNum = 50;
	}

	opt.smallRotNum = OP.getInt("smallRotNum");
	if (OP.fail()){
		cerr << "WARNING smallRotNum not specified, using default 5\n";
		opt.smallRotNum = 5;
	}

	if (opt.largeRotNum < 0 || opt.smallRotNum < 0) { cerr << "Need a positive rotamer number" << endl; exit(15); }

	return opt;
}
