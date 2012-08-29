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
#include "analEnergy.h"

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	// sys contains the structure found in the pdb
	System sys;
	sys.readPdb(opt.pdb);
	PolymerSequence pseq(sys);

	// outSys will be used to selectively calculate the energy terms
	System outSys;
	CharmmSystemBuilder CSB(outSys,opt.topfile,opt.parfile);

	cout << "Building system" << endl;

	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Assign the true coordiantes to the atoms in outSys
	int numAssignedAtoms = outSys.assignCoordinates(sys.getAtomPointers(),false);
	cout << "Number of assigned atoms: " << numAssignedAtoms << endl;

	// Build the all atoms without coordinates that are not in initial PDB
	outSys.buildAllAtoms();

	CSB.updateNonBonded(opt.cuton, opt.cutoff, opt.cutnb);

	// Print out total energy breakdown
	cout << "Calculate total energy: ";
	EnergySet *es =outSys.getEnergySet();
	cout << es->calcEnergy() << endl;

	// If both selections are listed, give pairwise interaction between them
	if ((opt.selection != "") && (opt.selection2 != "")) {
		AtomSelection sel(outSys.getAtomPointers());
		AtomPointerVector ats1 = sel.select(opt.selection);
		AtomPointerVector ats2 = sel.select(opt.selection2);

		CharmmEnergyCalculator calculator(opt.parfile);
		calculator.setNonBondedCutoffs(opt.cuton,opt.cutoff);

		map<string,double> ene =   calculator.calculatePairwiseNonBondedEnergy(ats1,ats2);

		map<string,double>::iterator it;
		cout << "Energy between selections in " << MslTools::getFileName(opt.pdb) << endl;
		for (it = ene.begin();it  != ene.end();it++){
			fprintf(stdout,"%-15s %8.3f\n",it->first.c_str(),it->second);
		}
	}

	EnergeticAnalysis ea;
	ea.setParameterFile(opt.parfile);
	ea.setPymolOutput(opt.pymolOutput);

	// Analyze specific positions
	if (opt.positions.size() > 0) {
	  for (uint i = 0; i < opt.positions.size();i++){
		int pos = outSys.getPositionIndex(opt.positions[i]); // takes "A_37" "A_37A" "A,37" or "A 37"
		cout << "Analyze "<<outSys.getResidue(pos).toString()<<endl;
		ea.analyzePosition(outSys, pos);
	  }
	}

	// If only one selection is given, analyze positions within this selection
	if (opt.selection != "" && opt.selection2 == ""){
	  ResidueSelection resSel(outSys);
	  vector<Residue *> list = resSel.select(opt.selection);
	  cout << "Selected " << list.size() << " residues from " << opt.selection << endl;
	  for (uint i = 0; i < list.size();i++){
	    int pos = outSys.getPositionIndex(list[i]->getPositionId());
	    cout << "Analyze "<<list[i]->toString()<<endl;
	    ea.analyzePosition(outSys, pos);
	  }
	  list.clear();
	}
}

