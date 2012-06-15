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
//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes
#include "minimize.h"
#include "System.h"
#include "SysEnv.h"
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
static SysEnv SYSENV;

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
	cout << "<cycles> no_of_cycles (default 1 - applicable only for constraint minimization) " << endl;
	cout << "<tolerance>  (default 0.01) " << endl;
	cout << "<dielectric> dielectric_constant_value(default 1)" << endl;
	cout << "<distanceDielectric> true or false (default false)" << endl;
	cout << "<constrainedatoms> selection_string_for_spring_constrained_atoms" << endl;
	cout << "<fixedatoms> selection_string_for_fixed_atoms" << endl;
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

	GSLMinimizer min(sys);
	
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
		cerr << " - NELDERMEAD2 (GSL v.1.14 or above)" << endl;
		cerr << " - NELDERMEAD2RAND (GSL v.1.14 or above)" << endl;
		cerr << " - STEEPEST_DESCENT " << endl;
		cerr << " - BFGS " << endl;
		cerr << " - CG_POLAK_RIBIERE " << endl;
		cerr << " - CG_FLETCHER_REEVES " << endl;
		exit(0);
	}

	min.setStepSize(opt.stepSize);
	min.setTolerance(opt.tolerance);
	min.setMaxIterations(opt.steps);

	if(opt.fixedAtoms != "") {
		cout << "Fixing atoms in selection " << opt.fixedAtoms << endl;
		min.fixAtoms(opt.fixedAtoms);
	}
	bool constrainedMinimization = false;
	if(opt.constrainedAtoms != "") {
		constrainedMinimization = true;
		cout << "Constraining atoms in selection " << opt.constrainedAtoms << endl;
		min.setConstraintForce(opt.constrainedAtoms,opt.springConstant);
	}

	sys.calcEnergy();
	cout << "Initial Energy .... " << endl;
	cout << sys.getEnergySummary() << endl;
	time_t start,end;
	time(&start);

	min.minimize();
	if(constrainedMinimization) {
		for(int i = 1; i < opt.cycles; i++) {
			cout << "After Cycle " << i << " .... " << endl;
			cout << sys.getEnergySummary() << endl;
			min.resetConstraints();
			min.minimize();
		}
	}
	time(&end);

	cout << "Final Energy .... " << endl;
	cout << sys.getEnergySummary() << endl;

	if(constrainedMinimization) { 
		min.removeConstraints();
	}

	sys.writePdb(opt.outfile);
	cout << opt.method << " Minimization Time: " << difftime(end,start) << " seconds" << endl;
	cout << "Written Minimized structure to " << opt.outfile << endl;

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
	dep[0].push_back("constrainedatoms");
	dep[0].push_back("springconstant");
	dep[1].push_back("springconstant");
	dep[1].push_back("constrainedatoms");
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
		string env = "MSL_CHARMM_PAR";
		if(SYSENV.isDefined(env)) {
			opt.parfile =  SYSENV.getEnv(env);	
		} else {
			opt.parfile = "/library/charmmTopPar/par_all22_prot.inp";
		}
		cerr << "WARNING 1111 no parfile specified using: "<<opt.parfile<<endl;
	}

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
		string env = "MSL_CHARMM_TOP";
		if(SYSENV.isDefined(env)) {
			opt.topfile =  SYSENV.getEnv(env);	
		} else {
			opt.topfile = "/library/charmmTopPar/top_all22_prot.inp";
		}

		cerr << "WARNING 1111 no topfile specified using: "<<opt.topfile<<endl;
	}
	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		opt.outfile = MslTools::getFileName(opt.pdb) + "_mini.pdb";
		cerr << "WARNING 1111 no outfile specified using: "<< opt.outfile<<endl;
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

	opt.cycles = OP.getInt("cycles");
	if (OP.fail()){
		opt.cycles = 1;
	}

	opt.tolerance = OP.getDouble("tolerance");
	if (OP.fail()){
		cerr << "WARNING no tolerance specified, using 0.01 as default.\n";
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
	opt.constrainedAtoms = MslTools::toUpper(OP.getString("constrainedatoms"));
	if (OP.fail()){
		opt.constrainedAtoms =  "";
	}
	opt.fixedAtoms = MslTools::toUpper(OP.getString("fixedatoms"));
	if (OP.fail()){
		opt.fixedAtoms =  "";
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
