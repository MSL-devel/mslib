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


// MSL Includes
#include "System.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "PhiPsiStatistics.h"
#include "getPseudoTorsion.h"

// STL Includes
#include<iostream>
using namespace std;



int main(int argc, char *argv[]){
    
    // Option Parser
    Options opt = setupOptions(argc,argv);

    // Read PDB
    System sys;
    sys.readPdb(opt.pdb);

    if (!sys.exists(opt.chain, opt.resi-1)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi-1<<" does not exist"<<endl;
	    exit(2333);
    }

    if (!sys.exists(opt.chain, opt.resi)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi<<" does not exist"<<endl;
	    exit(2333);
    }

    if (!sys.exists(opt.chain, opt.resi+1)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi+1<<" does not exist"<<endl;
	    exit(2333);
    }
    if (!sys.exists(opt.chain, opt.resi+2)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi+2<<" does not exist"<<endl;
	    exit(2333);
    }
    if (!sys.exists(opt.chain, opt.resi+3)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi+3<<" does not exist"<<endl;
	    exit(2333);
    }

    if (!sys.exists(opt.chain, opt.resi+4)){
	    cerr << "ERROR residue "<<opt.chain<<","<<opt.resi+4<<" does not exist"<<endl;
	    exit(2333);
    }

    // r0 and r5 are only there for phi/psi values. others r1-r4 compute pseudoTorsion
    Residue &r0 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi-1));
    Residue &r1 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi));
    Residue &r2 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi+1));
    Residue &r3 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi+2));
    Residue &r4 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi+3));
    Residue &r5 = sys.getResidue(sys.getPositionIndex(opt.chain,opt.resi+4));

    if (!(r1.exists("CA") && r2.exists("CA") && r3.exists("CA") && r4.exists("CA"))){
	    cerr << "ERROR a 'CA' atom does not exist on one of the residues."<<endl;
	    exit(2334);
    }

    


    double torsion = r1("CA").dihedral(r2("CA"),r3("CA"),r4("CA"));

    fprintf(stdout, "%10s %25s %s %3d %8.3f",MslTools::getFileName(opt.pdb).c_str(),opt.name.c_str(),r1.getChainId().c_str(),r1.getResidueNumber(),torsion);

    double phi     = PhiPsiStatistics::getPhi(r0,r1);
    double psi     = PhiPsiStatistics::getPsi(r1,r2);
    fprintf(stdout, " %8.3f %8.3f",phi,psi);


    phi     = PhiPsiStatistics::getPhi(r1,r2);
    psi     = PhiPsiStatistics::getPsi(r2,r3);
    fprintf(stdout, " %8.3f %8.3f",phi,psi);


    phi     = PhiPsiStatistics::getPhi(r2,r3);
    psi     = PhiPsiStatistics::getPsi(r3,r4);
    fprintf(stdout, " %8.3f %8.3f",phi,psi);

    phi     = PhiPsiStatistics::getPhi(r3,r4);
    psi     = PhiPsiStatistics::getPsi(r4,r5);
    fprintf(stdout, " %8.3f %8.3f\n",phi,psi);

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
		cout << "getDihedrals --pdb PDB --selection SEL [ --debug ]\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.chain = OP.getString("chain");
	if (OP.fail()){
		cerr << "ERROR 1111 chain not specified.\n";
		exit(1111);
	}
	opt.resi = OP.getInt("resi");
	if (OP.fail()){
		cerr << "ERROR 1111 resi not specified.\n";
		exit(1111);
	}
	opt.name = OP.getString("name");
	if (OP.fail()){
		cerr << "ERROR 1111 selection not specified.\n";
		exit(1111);
	}


	return opt;
}
