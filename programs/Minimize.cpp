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
//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes
#include "Minimize.h"
#include "System.h"
#include "OptionParser.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "GSLMinimizer.h"
#include "EnergySet.h"
#include "Timer.h"

using namespace std;
using namespace MSL;
using namespace MSL::MslTools;


int main(int argc, char *argv[]){
    
	// Option Parser
        Options opt = setupOptions(argc,argv);

	// Read PDB
	System initSys;
	initSys.readPdb(opt.pdb);

	// Create a sequence
	PolymerSequence pseq(initSys);

	CharmmSystemBuilder CSB(opt.topfile,opt.parfile);
	CSB.setDielectricConstant(opt.dielectric);
	CSB.setUseRdielectric(opt.distanceDependentElectrostatics);
	System sys;
	CSB.buildSystem(sys,pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).
	int numAssignedAtoms = sys.assignCoordinates(initSys.getAtomPointers());

	// Build the all atoms without coordinates (not in initial PDB)
	sys.buildAllAtoms();


	//CharmmEnergy::instance()->setDielectricConstant(opt.dielectric);
	//CharmmEnergy::instance()->setUseRdielectric(opt.distanceDependentElectrostatics);

        // Set atoms to minimize
	for (uint i =0 ; i < sys.getAtomPointers().size();i++){
		sys.getAtomPointers()(i).setSelectionFlag("all",true);
		sys.getAtomPointers()(i).setMinimizationIndex(i+1);
	}

	Timer t;
	double startTime = t.getWallTime();
	EnergySet *es = sys.getEnergySet();
        es->setTermActive("CHARMM_VDW",true);
        es->setTermActive("CHARMM_ELEC",false);
        es->setTermActive("CHARMM_BOND",false);
        es->setTermActive("CHARMM_ANGL",false);
        es->setTermActive("CHARMM_DIHE",false);
        es->setTermActive("CHARMM_U-BR",false);
	es->setTermActive("CHARMM_IMPR",false);

	GSLMinimizer min(es,&sys.getAtomPointers());
	//min.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	min.setMinimizeAlgorithm(GSLMinimizer::BFGS);
	min.setMaxIterations(opt.steps);
	min.Minimize();
	fprintf(stdout, "Time: %8.6f\n",(t.getWallTime()-startTime));

	/*
	for (uint r = 0; r < sys.residueSize();r++){

		Residue &res = sys.getResidue(r);

		int index = 0;
		for (uint i = 0; i < res.getAtomPointers().size();i++){
			res.getAtomPointers()(i).setMinimizationIndex(index+1);
			res.getAtomPointers()(i).setSelectionFlag("active",true);
			index++;
		}


		// Setup a minimizer object
		GSLMinimizer min(es,&res.getAtomPointers());
		min.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
		min.setMaxIterations(opt.steps);
		min.Minimize();

//		min.setMinimizeAlgorithm(GSLMinimizer::CG_FLETCHER_REEVES);
//		min.setMaxIterations(opt.steps);
//		min.Minimize();


		//char fname[80];
		//sprintf(fname,"%s_mini%03d.pdb",MslTools::getFileName(opt.pdb).c_str(),r);
		//sys.writePdb(fname);

		for (uint i = 0; i < res.getAtomPointers().size();i++){
			res.getAtomPointers()(i).setMinimizationIndex(-1);
			res.getAtomPointers()(i).setSelectionFlag("all",true);
		}
	}
	*/
	fprintf(stdout,"%s_mini.pdb %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),sys.getEnergySet()->calcEnergy());
	char fname[80];
	sprintf(fname,"%s_mini.pdb",MslTools::getFileName(opt.pdb).c_str());
	sys.writePdb(fname);


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
		cout << "Usage: Minimize conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.configFile = OP.getString("configfile");
	
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

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		opt.parfile =  "/library/charmmTopPar/par_all27_prot_lipid.inp";	
		cerr << "WARNING 1111 no parfile specified using: "<<opt.parfile<<endl;
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		opt.topfile =  "/library/charmmTopPar/top_all27_prot_lipid.inp";
		cerr << "WARNING 1111 no topfile specified using: "<<opt.topfile<<endl;
	}

	opt.steps = OP.getInt("steps");
	if (OP.fail()){
		opt.steps = 50;
	}

	opt.dielectric = OP.getDouble("dielectric");
	if (OP.fail()){
		cerr << "WARNING no dielectric specified, using 80 as default.\n";
		opt.dielectric = 80;

	}
	opt.distanceDependentElectrostatics = OP.getBool("distanceDielectric");
	if (OP.fail()){
		cerr << "WARNING no distanceDielectric specified, using 'false' as default.\n";
		opt.distanceDependentElectrostatics = false;
	}

	return opt;
}
