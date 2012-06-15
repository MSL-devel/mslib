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



#include "runKBQuench.h"
#include "Quench.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "System.h"

using namespace std;

using namespace MSL;



int main(int argc, char *argv[]){

	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	TwoBodyDistanceDependentPotentialTable tbd;
	tbd.readPotentialTable(opt.potfile);

	System initialSystem;
	initialSystem.readPdb(opt.pdb);

	Quench quencher(opt.topfile, opt.parfile, opt.rotlib);

	// Set the number of rotamers you want for large and small side chains, respectively
	quencher.setVariableNumberRotamers(opt.largeRotNum,opt.smallRotNum);

	// Set which positions you want to repack (optional)
	if (opt.positions.size() > 0) {	
		vector<int> variablePositions;
		for (uint i = 0; i < opt.positions.size(); i++) {
			//vector<string> pos = MslTools::tokenize(opt.positions[i],"_");
			variablePositions.push_back(initialSystem.getPositionIndex(opt.positions[i]));
		}
		System sys = quencher.runQuench(initialSystem,variablePositions,tbd);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtomPointers())) {
			cerr << "Problem writing " << opt.outfile << endl;
		}
		writer.close();
	}
	else {
		System sys = quencher.runQuench(initialSystem,tbd);
		cout << "Write pdb " << opt.outfile << endl;
		PDBWriter writer;
		writer.open(opt.outfile);
		if (!writer.write(sys.getAtomPointers())) {
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
		cout << "runKBQuench --pdb PDB [--potfile KBFILE --topfile TOPFILE --parfile PARFILE --rotlib ROTLIB --outfile OUTPDB --positions A_73 A_74 A_75 --smallRotNum NUM --largeRotNum NUM]\n";
		exit(0);
	}

	opt.pdb  = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.potfile = OP.getString("potfile");
	if (OP.fail()){
		cerr << "WARNING no potfile specified, using default /snap/cluster/jedonald/projects/KB_Spring2009/convert2MSL/DPHER_dfalpha1.61_rcut14.5_lcalpha0.5.equate.tbd\n";
		opt.potfile = "/snap/cluster/jedonald/projects/KB_Spring2009/convert2MSL/DPHER_dfalpha1.61_rcut14.5_lcalpha0.5.equate.tbd";
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
