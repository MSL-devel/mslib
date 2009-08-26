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
#include "generateCrystalLattice.h"
#include "OptionParser.h"
#include "CrystalLattice.h"
#include "MslTools.h"

using namespace std;

int main(int argc, char *argv[]){

	// Option Parser
	Options opt = setupOptions(argc,argv);


	// Crystal Lattice Object
	CrystalLattice cl(opt.pdb);
	
	cl.generateCrystal();
	cl.writeCrystalUnits(opt.outfile,true,opt.singleFile,opt.renameChainsExcept);

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
		cout << "Usage: generateCrystalLattice conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << "#outfile FILE\n";
		cout << "#singleFile 1\n";
		cout << "#renameChains Z\n";
		cout << endl;
		exit(0);
	}

	opt.configFile = OP.getString("config");
	
	if (opt.configFile != "") {
		OP.readFile(opt.configFile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.outfile = OP.getString("out");
	if (OP.fail()){


		opt.outfile = MslTools::getFileName(opt.pdb);

		cout << "WARNING 1111 no out file specified, "<<opt.outfile<<" will be used as a prefix"<<endl;

	}

	opt.singleFile = OP.getBool("singleFile");
	if (OP.fail()){
		opt.singleFile = false;
	}
	opt.nmrStyleFile = OP.getBool("nmrStyleFile");
	if (OP.fail()){
		opt.nmrStyleFile = false;
		opt.singleFile   = true;
	}

	opt.renameChainsExcept = OP.getString("renameChainsExcept");
	if (OP.fail()){
		opt.renameChainsExcept = "";
	}

	return opt;
}
