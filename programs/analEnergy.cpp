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

#include "EnergeticAnalysis.h"
#include "CharmmSystemBuilder.h"
#include "PolymerSequence.h"
#include "OptionParser.h"
#include "analEnergy.h"

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	System sys;
	sys.readPdb(opt.pdb);

	PolymerSequence pseq(sys);

	CharmmSystemBuilder CSB(opt.topfile,opt.parfile);

	System outSys;
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(outSys,pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	int numAssignedAtoms = outSys.assignCoordinates(sys.getAtoms(),false);
	fprintf(stdout,"Number of assigned atoms: %d",numAssignedAtoms);

	// Build the all atoms without coordinates (not in initial PDB)
	outSys.buildAllAtoms();

	outSys.writePdb("/tmp/preEA.pdb");
	EnergeticAnalysis ea;

	for (uint i = 0; i < opt.positions.size();i++){

		vector<string> toks = MslTools::tokenize(opt.positions[i],"_");
		if (toks.size() != 2){
			cerr << "Position specification no good: "<<opt.positions[i]<<endl;
			continue;
		}
		if (!outSys.exists(toks[0],MslTools::toInt(toks[1]))){
			cerr << "Position: "<<opt.positions[i]<<" does not exist!"<<endl;
			continue;
		}
		int pos = outSys.getPositionIndex(toks[0],MslTools::toInt(toks[1]));
		cout << "Analyze "<<outSys.getResidue(pos).toString()<<endl;
		ea.analyzePosition(outSys, pos);
	}
}

Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;
	
	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: getSphericalCoordinates conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.configfile = OP.getString("configfile");
	
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configfile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.positions = OP.getStringVector("positions");
	if (OP.fail()){
		cerr << "ERRROR 1111 no chain\n";
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

	return opt;
}

