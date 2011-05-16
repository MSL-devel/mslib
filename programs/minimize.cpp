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
#include "AtomSelection.h"
#include "GSLMinimizer.h"
#include "EnergySet.h"
#include "Timer.h"

using namespace std;
using namespace MSL;
using namespace MSL::MslTools;

void usage() {
	cout << "Usage: Minimize conf" << endl;
	cout << endl;
	cout << "Required Parameters" << endl;
	cout << "<pdb> PDB" << endl;
	cout << endl;
	cout << "Optional Parameters" << endl;
	cout << "<topfile> CHARMM_TOPOLOGY_FILE" << endl;
	cout << "<parfile> CHARMM_PARAMETER_FILE" << endl;
	cout << "<method> Minimization_algorithm_to_use" << endl;
	cout << "<steps> no_of_steps (default 50) " << endl;
	cout << "<stepsize>  Initial trial stepsize (default 0.1) " << endl;
	cout << "<tolerance>  (default 0.01) " << endl;
	cout << "<dielectric> dielectric_constant_value(default 1)" << endl;
	cout << "<distanceDielectric> true or false (default false)" << endl;
	cout << "<selection> selection_string_for_spring_constrained_atoms" << endl;
	cout << "<springconstant> specified_during_constrained_minimization" << endl;
	cout << "<cuton> non-bonded_interactions" << endl;
	cout << "<cutoff> non-bonded_interactions" << endl;
	cout << "<cutnb> non_bonded_interactions" << endl;

}

int main(int argc, char *argv[]){
    
	// Option Parser
        Options opt = setupOptions(argc,argv);

	// Read PDB
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topfile,opt.parfile);
	if(opt.cutnb != -1.0) {
		CSB.setBuildNonBondedInteractions(false);
	}
	CSB.setDielectricConstant(opt.dielectric);
	CSB.setUseRdielectric(opt.distanceDependentElectrostatics);
	if(!CSB.buildSystemFromPDB(opt.pdb)) {
		cerr << "Unable to build System form Pdb " << opt.pdb << endl;
		exit(0);
	}

	sys.buildAllAtoms();
	if(opt.cutnb != -1.0) {
		CSB.updateNonBonded(opt.ctonnb,opt.ctofnb,opt.cutnb);
	}

        // Set atoms to minimize
	for (uint i =0 ; i < sys.getAtomPointers().size();i++){
		sys.getAtomPointers()(i).setMinimizationIndex(i+1);
	}


	EnergySet *es = sys.getEnergySet();
	es->setAllTermsActive();
        //es->setTermActive("CHARMM_VDW",false);
        //es->setTermActive("CHARMM_ELEC",false);
        //es->setTermActive("CHARMM_BOND",false);
        //es->setTermActive("CHARMM_ANGL",false);
        //es->setTermActive("CHARMM_DIHE",false);
        //es->setTermActive("CHARMM_U-BR",false);
	//es->setTermActive("CHARMM_IMPR",false);

	GSLMinimizer min(*es,sys.getAtomPointers());
	
	if(opt.method == "NELDERMEAD1") {
		min.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD1);
	} else if (opt.method == "NELDERMEAD2") {
		min.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD2);
	} else if (opt.method == "NELDERMEAD2RAND") {
		min.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD2RAND);
	} else if (opt.method == "STEEPEST_DESCENT") {
		min.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	} else if (opt.method == "BFGS") {
		min.setMinimizeAlgorithm(GSLMinimizer::BFGS);
	} else if (opt.method == "CG_POLAK_RIBIERE") {
		min.setMinimizeAlgorithm(GSLMinimizer::CG_POLAK_RIBIERE);
	} else if (opt.method == "CG_FLETCHER_REEVES") {
		min.setMinimizeAlgorithm(GSLMinimizer::CG_FLETCHER_REEVES);
	} else {
		cerr << "Unknown method " << opt.method << endl;
		cerr << "Should be one of: " << endl;
		cerr << " - NELDERMEAD1 " << endl;
		cerr << " - NELDERMEAD2 " << endl;
		cerr << " - NELDERMEAD2RAND " << endl;
		cerr << " - STEEPEST_DESCENT " << endl;
		cerr << " - BFGS " << endl;
		cerr << " - CG_POLAK_RIBIERE " << endl;
		cerr << " - CG_FLETCHER_REEVES " << endl;
		exit(0);
	}

	min.setStepSize(opt.stepSize);
	min.setTolerance(opt.tolerance);
	min.setMaxIterations(opt.steps);

	AtomSelection sel(sys.getAtomPointers());
	AtomPointerVector springAtoms;
	if(opt.selection != "") {
		springAtoms = sel.select(opt.selection);
		cout << "Spring Atoms " << springAtoms.size() << endl;
		min.constrainAtoms(springAtoms,opt.springConstant);
	}

	sys.calcEnergy();
	cout << sys.getEnergySummary() << endl;
	time_t start,end;
	time(&start);
	min.Minimize();
	time(&end);
	if(opt.selection != "") {
		cout << sys.getEnergySummary() << endl;
		min.removeConstraints();
	}

	string outfile = MslTools::getFileName(opt.pdb) + "_mini.pdb";
	sys.writePdb(MslTools::getFileName(opt.pdb) + "_mini.pdb");
	cout << opt.method << " Time: " << difftime(end,start) << " seconds" << endl;;
	cout << "Final Energy: " << sys.calcEnergy() << endl;
	cout << "Written Minimized structure to " << outfile << endl;
	cout << sys.getEnergySummary() << endl;


}


Options setupOptions(int theArgc, char * theArgv[]){

	// Create the options
	Options opt;
	
	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	vector<vector<string> > dep(2);
	dep[0].push_back("selection");
	dep[0].push_back("springconstant");
	dep[1].push_back("springconstant");
	dep[1].push_back("selection");
	OP.setDependentOn(dep);
	
	vector<vector<string> > interd(1);
	interd[0].push_back("cuton");
	interd[0].push_back("cutoff");
	interd[0].push_back("cutnb");
	OP.setInterDependent(interd);

	if (OP.countOptions() == 0){
		usage();
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
	if(!OP.checkOptions()) {
		cout << OP.getErrors() << endl;
		exit(1111);
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

	opt.stepSize = OP.getDouble("stepsize");
	if (OP.fail()){
		cerr << "WARNING no stepsize specified, using 0.1 as default.\n";
		opt.stepSize = 0.1;
	}

	opt.tolerance = OP.getDouble("tolerance");
	if (OP.fail()){
		cerr << "WARNING no stepsize specified, using 0.01 as default.\n";
		opt.tolerance = 0.01;
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
	opt.method = MslTools::toUpper(OP.getString("method"));
	if (OP.fail()){
		opt.method =  "BFGS";
		cerr << "WARNING 1111 no method specified. Using: "<< opt.method <<endl;
	}
	opt.selection = MslTools::toUpper(OP.getString("selection"));
	if (OP.fail()){
		opt.selection =  "";
	}
	opt.springConstant = OP.getDouble("springconstant");
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
	return opt;
}
