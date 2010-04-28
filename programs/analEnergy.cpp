/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
 Sabareesh Subramaniam, Ben Mueller

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
#include "ResidueSelection.h"
#include "analEnergy.h"
#include "release.h"
#include "AtomicPairwiseEnergy.h"
#include "AtomSelection.h"

using namespace MSL;
using namespace std;


int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	System outSys;
//	CharmmSystemBuilder CSB(outSys,opt.topfile,opt.parfile);
//	CSB.buildSystemFromPDB(opt.pdb);

	

	System sys;
	sys.readPdb(opt.pdb);
	for (uint i = 0; i < sys.positionSize();i++){
	  if (sys.getResidue(i).getResidueName() == "HIS"){
	    sys.getResidue(i).setResidueName("HSD");
	  }

	}
	PolymerSequence pseq(sys);

	CharmmSystemBuilder CSB(outSys,opt.topfile,opt.parfile,"");

	fprintf(stdout, "Building system");
	fflush(stdout);

	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).


	int numAssignedAtoms = outSys.assignCoordinates(sys.getAtomPointers(),false);
	fprintf(stdout,"Number of assigned atoms: %d",numAssignedAtoms);
	fflush(stdout);
	// Build the all atoms without coordinates (not in initial PDB)
	outSys.buildAllAtoms();

	CSB.updateNonBonded(8.0, 9.0, 12.0);

	//outSys.writePdb("/tmp/preEA.pdb");


	// If no positions or selection to analyze, then print out total energy breakdown!
	cout << "Calculate total energy\n";
	EnergySet *es =outSys.getEnergySet();

	if (opt.selection2 != ""){
	  AtomSelection sel(outSys.getAtomPointers());
	  AtomPointerVector ats1 = sel.select(opt.selection);
	  AtomPointerVector ats2 = sel.select(opt.selection2);

	  AtomicPairwiseEnergy ape(opt.parfile);
	  ape.setNonBondedCutoffs(8,12); // 0->8 is full 8-12 is switched, >12 is 0.

	  map<string,double> ene =   ape.calculatePairwiseNonBondedEnergy(outSys,ats1,ats2);

	  map<string,double>::iterator it;
	  for (it = ene.begin();it  != ene.end();it++){
	    fprintf(stdout,"%30s %-10s: %-15s %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),opt.prefix.c_str(),it->first.c_str(),it->second);
	  }

	  /*
		es->setAllTermsInactive();
		es->setTermActive("CHARMM_VDW");		
		double vdw = es->calcEnergy(opt.selection,opt.selection2);

		es->setAllTermsInactive();
	        es->setTermActive("CHARMM_ELEC");
	    	double elec = es->calcEnergy(opt.selection,opt.selection2);

		fprintf(stdout,"CHARMM VDW ELEC: %8.3f %8.3f\n",vdw,elec);  
	  */
	  
		
	}

	//cout << outSys.calcEnergy() << endl;
	//cout << outSys.getEnergySummary();	


	EnergeticAnalysis ea;
	ea.setParameterFile(opt.parfile);

	if (opt.positions.size() > 0) {
	  for (uint i = 0; i < opt.positions.size();i++){
		int pos = outSys.getPositionIndex(opt.positions[i]); // takes "A_37" "A_37A" "A,37" or "A 37"
		cout << "Analyze "<<outSys.getResidue(pos).toString()<<endl;
		ea.analyzePosition(outSys, pos);
	  }
	}

	if (opt.selection != "" && opt.selection2 == ""){
	  ResidueSelection resSel(outSys);
	  vector<Residue *> list = resSel.select(opt.selection);
	  fprintf(stdout, "Selected %d residues with '%s'\n",(int)list.size(),opt.selection.c_str());
	  fflush(stdout);
	  for (uint i = 0; i < list.size();i++){
	    int pos = outSys.getPositionIndex(list[i]->getPositionId());
	    cout << "Analyze "<<list[i]->toString()<<endl;
	    ea.analyzePosition(outSys, pos);
	  }
	  list.clear();
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
		cout << "Usage: analEnergy conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << "positions CHAIN_RESNUM\n";
		cout << "select SEL \n";
		cout << "# topfile /library/charmmTopPar/top_all22_prot.inp";
		cout << "# parfile /library/charmmTopPar/par_all22_prot.inp";
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
		cerr << "WARNING 1111 no positions to analyze\n";
	}	
	
	opt.selection = OP.getString("select");
	if (OP.fail()){
	  cerr << "WARNING 1111 no selections to analyze\n";
	}

	opt.selection2 = OP.getString("select2");
	
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

	opt.prefix = OP.getString("prefix");
	if (OP.fail()){
	  opt.prefix = "ENE";
	}

	return opt;
}

