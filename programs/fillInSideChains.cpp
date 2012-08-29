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
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "PolymerSequence.h"
#include "Position.h"
#include "Residue.h"
#include "MslTools.h"
#include "CharmmEnergyCalculator.h"
#include "fillInSideChains.h"

// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace MSL;
using namespace std;



int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	cout << "Read in PDB structure"<<endl;
	System initSys;
	initSys.readPdb(opt.pdb);

	initSys.writePdb("/tmp/init.pdb");

	// Read through and mark positions that don't have full side chains?
	CharmmTopologyReader CTR(opt.topfile);
	CTR.read();

	vector<int> residuesToFill;
	for (uint r = 0; r < initSys.positionSize();r++){

		Residue &res = initSys.getResidue(r);

		// Default to HSD for HIS
		if (res.getResidueName() == "HIS") res.setResidueName("HSD");

		if (!CTR.residueExists(res.getResidueName())){
			cerr << "ERROR 2222 Residue: "<<res.toString()<<" ; residue type does not exist in topology file: "<<opt.topfile<<endl;
			exit(2222);
		}
		

		vector<string> topologyAtoms = CTR.getResidue(res.getResidueName()).getAllTopoAtomNames();

		for (uint t = 0; t < topologyAtoms.size();t++){
			
			// Skip over hydrogen atoms..
			if (topologyAtoms[t].substr(0,1) == "H") continue;

			// If ILE and CD, look for CD1
			if (res.getResidueName() == "ILE" && topologyAtoms[t] == "CD" && !res.atomExists(topologyAtoms[t]) && res.atomExists("CD1")){
				//res.getAtom("CD1").setName("CD");
				//res.updateAtomMap(res.getAtom("CD1"));
				continue;
			}


			if (!res.atomExists(topologyAtoms[t])){
				
				residuesToFill.push_back(r);
				cout << "*** Residue to fill: "<<res.toString()<<" due to atom from topology not found in pdb: "<<topologyAtoms[t]<<endl;
				break;
			}
		}

	}


	PolymerSequence seq(initSys);

	cout << "Build Charmm System"<<endl;
	System sys;

	CharmmSystemBuilder CSB(sys,opt.topfile,opt.parfile);
	CSB.setBuildNonBondedInteractions(false);
	CSB.buildSystem(seq);



	int numAssignedAtoms = sys.assignCoordinates(initSys.getAtomPointers());
	fprintf(stdout, "\tNumber of assigned atoms: %8d\n",numAssignedAtoms);
	if (numAssignedAtoms == 0){
		cerr << "ERROR 2222 zero assigned atoms, means after re-building from sequence we can not match chain,residue numbers with original coordinates.\n";
		exit(2222);
	}

	sys.buildAllAtoms();


	if (opt.debug) {
		string filename = "/tmp/initialBuild.pdb";
		cout << "Write initial build pdb " << filename << endl;
		PDBWriter writer;
		writer.open(filename);
		if (!writer.write(sys.getAtomPointers())) {
			cerr << "Problem writing " << filename << endl;
		}
		writer.close();
	}


	cout << "Read rotamer library " << opt.rotlib << " and load rotamers"<<endl;	

	SystemRotamerLoader sysRot(sys, opt.rotlib);


	// For each residue to fill add 100 rotamers
	for (uint r = 0; r < residuesToFill.size();r++){

		Residue  &res = sys.getResidue(residuesToFill[r]);
		Position *pos = res.getParentPosition();

		//sysRot.loadRotamers(pos, "BALANCED-200",res.getResidueName(),0,99);
		sysRot.loadRotamers(pos,res.getResidueName(),0,99,"");

	}

	
	// Quencher-type
	CharmmEnergyCalculator calculator(opt.parfile);
	for (uint p = 0; p < sys.positionSize();p++){
		Position &pos = sys.getPosition(p);

		if (pos.getTotalNumberOfRotamers() > 1){


			double minEnergy = MslTools::doubleMax;
			int    minIndex  = -1;
			for (uint c = 0; c < pos.getTotalNumberOfRotamers();c++){
				double self = calculator.calculateSelfEnergy(sys,p,c);
				double temp = calculator.calculateTemplateEnergy(sys,p,c);	


				if (temp+self < minEnergy){
					minEnergy = temp+self;
					minIndex  = c;
				}

			}

			if (minIndex == -1){
				cerr << "ERROR 8754 on position "<<pos.getCurrentIdentity().toString()<< " no low energy rotamer ? .. do nothing and continue."<<endl;
				continue;
			}

			// Set to lowest self+template energy rotamer.
			cout << "*** Setting "<<pos.getCurrentIdentity().toString()<<" to rotamer "<<minIndex<<" with self+template energy = "<<minEnergy<<endl;
			pos.setActiveRotamer(minIndex);

		}
	}


	
	sys.writePdb(opt.outpdb);
	
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
        cout << "fillInSideChains CONF\n";
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
	    cout << "outpdb /tmp/out.pdb"<<endl<<endl;
	    cout << "# Rotamer library"<<endl;
	    cout << "rotlib /library/rotlib/balanced/rotlib-balanced-200.txt"<<endl<<endl;
	    cout << "# CHARMM Topology"<<endl;
	    cout << "topfile /library/charmmTopPar/top_all27_prot_lipid.inp"<<endl<<endl;
	    cout << "# CHARMM Parameter"<<endl;
	    cout << "parfile /library/charmmTopPar/par_all27_prot_lipid.inp"<<endl<<endl;
	    exit(1);
    }


    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	    cerr << "ERROR 1111 no pdb file"<<endl;
    }

    opt.rotlib = OP.getString("rotlib");
    if (OP.fail()){
	    cerr << "ERROR 1111 no rotlib file"<<endl;
	    
    }

    opt.topfile = OP.getString("topfile");
    if (OP.fail()){
	    cerr << "ERROR 1111 no topfile file"<<endl;
    }
    opt.parfile = OP.getString("parfile");
    if (OP.fail()){
	    cerr << "ERROR 1111 no parfile file"<<endl;
    }

    opt.outpdb = OP.getString("outpdb");
    if (OP.fail()){
	    opt.outpdb = "/tmp/out.pdb";
	    cerr << "WARNING no outpdb file specifed will write to: "<<opt.outpdb<<endl;
    }

    opt.debug = OP.getBool("debug");
    if (OP.fail()){
	    opt.debug = false;
    }


    
    return opt;
}
