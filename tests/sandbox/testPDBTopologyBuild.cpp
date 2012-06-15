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


#include <iostream>

#include "System.h"
#include "PDBTopologyBuilder.h"
#include "PDBWriter.h"
#include "AtomSelection.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;


int main() {

	System sys;

	PolymerSequence seq("\
A: ALA ILE VAL ILE\n\
B: ARG HIS THR GLY");


	cout << seq << endl;

	cout << "Build from PDB topology" << endl;
	PDBTopologyBuilder builder(sys, "/library/charmmTopPar/top_pdb2.3_noH.inp");
	builder.buildSystem(seq);
	sys.printIcTable();

	/*
	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}
	*/
	if (!sys.seed()) {
		cerr << "Cannot seed the system" << endl;
	}
	sys.buildAtoms();
	AtomSelection sel(sys.getAllAtomPointers());
	sel.select("chainB, chain B");

	Transforms tr;
	tr.translate(sel.getSelection("chainB"), CartesianPoint(10,5,5));
	string filename = "/tmp/buildFromPDBTopology.pdb";

	if (!sys.writePdb(filename)) {
		cerr << "Cannot write output file " << filename << endl;
		exit(1);
	}

	cout << "Written pdb file " << filename << endl;
	cout << endl;

	cout << sys.getAtomPointers();

	exit(0);
	cout << "=====================================================" << endl;
	cout << "Create a new System from the PDB we previously saved with the buildSystemFromPDB function" << endl;
	
	System sys2;
	PDBTopologyBuilder builder2(sys2, "/library/charmmTopPar/top_pdb2.3_noH.inp");
	builder2.buildSystemFromPDB("/tmp/buildFromPDBTopology.pdb");

	AtomSelection sel2(sys2.getAtomPointers());
	AtomPointerVector noCrd = sel2.select("noCrd, HASCOOR 0");

	if (noCrd.size() != 0) {
		cout << "Error building from PDB, there are " << noCrd.size() << " atoms without coordinates out of " << sys2.atomSize() << endl;
	} else { 
		cout << "System build OK from PDB " << filename << endl;
	}

	filename = "/tmp/buildFromPDBTopology2.pdb";

	if (!sys2.writePdb(filename)) {
		cerr << "Cannot write output file " << filename << endl;
		exit(1);
	}

	cout << "Written pdb file " << filename << endl;
	cout << endl;

	return 0;
}
