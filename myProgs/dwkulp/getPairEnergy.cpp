/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
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
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "Timer.h"
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "EnergySet.h"
#include "SysEnv.h"
#include "getPairEnergy.h"


using namespace std;

using namespace MSL;

// Global objects
Timer t;
double startTime = 0.0;
static SysEnv SYSENV;

int main(int argc, char *argv[]){

	// Option Parser
	Options opt = setupOptions(argc, argv);

	// Create a system from the structural input options
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topfile,opt.parfile);


        //CSB.setBuildNonBondedInteractions(false);
	CSB.setDielectricConstant(opt.dielectric);
	CSB.setUseRdielectric(opt.distanceDependentElectrostatics);
	if(!CSB.buildSystemFromPDB(opt.pdb)) {
		cerr << "Unable to build" << opt.pdb << endl;
		exit(0);
	}

	// Add Side Chains
	sys.buildAllAtoms();
	CSB.updateNonBonded(opt.ctonnb,opt.ctofnb,opt.cutnb);

	sys.calcEnergy();
	EnergySet *e = sys.getEnergySet();

	HydrogenBondBuilder HBB(sys,opt.hbondfile);
	HBB.buildInteractions(opt.cuthb);  

	//cout << sys.getEnergySummary();

	Position &pos1 = sys.getPosition(opt.position1);
	for (uint p = 0; p < pos1.atomSize();p++){
	  pos1.getAtom(p).setSelectionFlag("pos1",true);
	}

	Position &pos2 = sys.getPosition(opt.position2);
	for (uint p = 0; p < pos2.atomSize();p++){
	  pos2.getAtom(p).setSelectionFlag("pos2",true);
	}
	
	map<string,double> weights = e->getWeightMap();
	map<string,double> energies;
	e->setAllTermsInactive();
	for (map<string,double>::iterator it = weights.begin();it != weights.end();it++){
	  string term = it->first;
	  e->setTermActive(term,true);
	  energies[term] = e->calcEnergy("pos1","pos2");
	  e->setTermActive(term,false);
	}
	for (map<string,double>::iterator it = energies.begin();it != energies.end();it++){
	  fprintf(stdout,"%10s %-15s %10.4f\n",MslTools::getFileName(opt.pdb).c_str(),it->first.c_str(), it->second);
	}


}



Options setupOptions(int theArgc, char * theArgv[]){
    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile optionn
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: getPairEnergy " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdb PDB\n";
	cout << "position1 POSITION_DEF\n";
	cout << "position2 POSITION_DEF\n";
	cout << endl;
	exit(0);
    }

    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdb specified."<<endl;
	exit(1111);
    }
    opt.position1 = OP.getString("position1");
    if (OP.fail()){
	cerr << "ERROR 1111 no position1 specified."<<endl;
	exit(1111);
    }
    opt.position2 = OP.getString("position2");
    if (OP.fail()){
	cerr << "ERROR 1111 no position2 specified."<<endl;
	exit(1111);
    }
    opt.topfile = OP.getString("topfile");
    if (OP.fail()){
      string envVar = "MSL_CHARMM_TOP";
      if (SYSENV.isDefined(envVar) && MslTools::fileExists(SYSENV.getEnv(envVar))){
        opt.topfile = SYSENV.getEnv(envVar);
	cerr << "WARNING 1111 no topfile specified using: "<<opt.topfile<<endl;
      } else {
	cerr << "ERROR 2222 no topfile defined and default("<<SYSENV.getEnv(envVar)<<") does not exist."<<endl;
	exit(2222);
      }
    }

    opt.parfile = OP.getString("parfile");
    if (OP.fail()){
      string envVar = "MSL_CHARMM_PAR";
      if (SYSENV.isDefined(envVar) && MslTools::fileExists(SYSENV.getEnv(envVar))){
        opt.parfile = SYSENV.getEnv(envVar);
	cerr << "WARNING 1111 no parfile specified using: "<<opt.parfile<<endl;
      } else {
	cerr << "ERROR 2222 no parfile defined and default("<<SYSENV.getEnv(envVar)<<") does not exist."<<endl;
	exit(2222);
      }
    }

    opt.hbondfile = OP.getString("hbondfile");
    if (OP.fail()){
      string envVar = "MSL_HBOND_CA_PAR";
      if (SYSENV.isDefined(envVar) && MslTools::fileExists(SYSENV.getEnv(envVar))){
        opt.hbondfile = SYSENV.getEnv(envVar);
	cerr << "WARNING 1111 no hbondfile specified using: "<<opt.hbondfile<<endl;
      } else {
	cerr << "ERROR 2222 no hbondfile defined and default("<<SYSENV.getEnv(envVar)<<") does not exist."<<endl;
	exit(2222);
      }
    }


    opt.dielectric = OP.getDouble("dielectric");
    if (OP.fail()){
        cerr << "WARNING no dielectric specified, using 1 as default.\n";
	opt.dielectric = 1.0;

    }
    opt.distanceDependentElectrostatics = OP.getBool("distanceDielectric");
    if (OP.fail()){
        cerr << "WARNING no distanceDielectric specified, using 'true' as default.\n";
        opt.distanceDependentElectrostatics = true;
    }

    opt.ctonnb = OP.getDouble("cuton");
    if (OP.fail()){
      opt.ctonnb =  -1.0;
    }
    opt.ctofnb = OP.getDouble("cutoff");
    if (OP.fail()){
      opt.ctofnb =  -1.0;
    }
    opt.cutnb = OP.getDouble("cutnb");
    if (OP.fail()){
      opt.cutnb =  -1.0;
    }
    opt.cuthb = OP.getDouble("cuthb");
    if(OP.fail()) {
      opt.cuthb = 10.0;
    }
    return opt;
}
