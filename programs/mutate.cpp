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

// MSL Includes
#include "System.h"
#include "OptionParser.h"
#include "MslTools.h"
#include "PolymerSequence.h"
#include "Position.h"
#include "Residue.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "mutate.h"

// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace MSL;
using namespace std;



int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	cout << "Read in PDB structure"<<endl;
	System sys;
	sys.readPdb(opt.pdb);

	// Get one of the residues 
	Position &pos = sys.getPosition(opt.position);

	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);

	// Get backbone atoms  
	AtomPointerVector backboneAtoms = pdbTop.getBackboneAtoms(pos.getCurrentIdentity());
	
	// *********  MUTATE *************** //
	stringstream newPos;
	newPos <<  opt.position <<","<<opt.newRes;
	AtomContainer newAtoms = pdbTop.getResidue(newPos.str(),backboneAtoms,opt.numRot);

	pos.addIdentity(newAtoms.getAtomPointers(),opt.newRes);
	pos.setActiveIdentity(opt.newRes);

	for (uint i = 0; i < pos.getTotalNumberOfRotamers();i++){
		pos.setActiveRotamer(i);
		char tmp[100];
		sprintf(tmp,"%s-%03d.pdb",opt.outpdb.c_str(),i);
	        string outfile(tmp);
		sys.writePdb(outfile);
	}

}

Options setupOptions(int theArgc, char * theArgv[]){

    // Create the options
    Options opt;
    

    // Parse the options
    OptionParser OP;
    OP.readArgv(theArgc, theArgv);
    OP.setRequired(opt.required);    
    OP.setAllowed(opt.optional);
    //    OP.setShortOptionEquivalent(opt.equivalent);
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"


    if (OP.countOptions() == 0){
        cout << "Usage:" << endl;
        cout << endl;
        cout << "mutate CONF\n";
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

    if (OP.getBool("help")){

	    cout << "# PDB "<<endl;	
	    cout << "pdb foo.pdb"<<endl<<endl;
	    cout << "# Output pdb" <<endl;
	    cout << "outpdb /tmp/out"<<endl<<endl;
	    cout << "# Rotamer library"<<endl;
	    cout << "rotlib /library/rotlib/balanced/rotlib-balanced-200.txt"<<endl<<endl;
	    cout << "# Position"<<endl;
	    cout << "position B,2"<<endl;
	    cout << "# New residue"<<endl;
	    cout << "newRes PHE"<<endl<<endl;
	    cout << "# Number of rotamers"<<endl;
	    cout << "numRot 10"<<endl<<endl;
	    exit(1);
    }


    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	    cerr << "ERROR 1111 no pdb file"<<endl;
    }

    opt.rotlib = OP.getString("rotlib");
    if (OP.fail()){
	    cerr << "ERROR 1111 no rotlib file"<<endl;
	    exit(1111);
    }

    opt.position = OP.getString("position");
    if (OP.fail()){
	    cerr << "ERROR 1111 no position"<<endl;
	    exit(1111);
    }
    opt.newRes = OP.getString("newRes");
    if (OP.fail()){
	    cerr << "ERROR 1111 no new residue (newRes)"<<endl;
	    exit(1111);
    }

    opt.outpdb = OP.getString("outpdb");
    if (OP.fail()){
	    opt.outpdb = "/tmp/out";
	    cerr << "WARNING no outpdb file specifed will write to: "<<opt.outpdb<<endl;
    }

    opt.numRot = OP.getInt("numRot");
    if (OP.fail()){
	    cerr << "numRot defaults to 10"<<endl;
	    opt.numRot = 10;
    }
    opt.debug = OP.getBool("debug");
    if (OP.fail()){
	    opt.debug = false;
    }


    
    return opt;
}
