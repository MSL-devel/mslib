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

#include <iostream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"

using namespace std;

int main() {

	System sys;

	PolymerSequence seq("\
A: ALA-ACE ILE VAL ILE\n\
B: ARG HSD THR GLY");

/*
	PolymerSequence seq("\
A: GLY GLY\n\
B: ARG HSD THR GLY");
*/

	cout << seq << endl;

	//CharmmSystemBuilder CSB("/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp");
	CharmmSystemBuilder CSB("/library/charmmTopPar/top_all22_prot.inp", "/library/charmmTopPar/par_all22_prot.inp");
	CSB.buildSystem(sys, seq);
	sys.printIcTable();

	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
	}
	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
	}
	sys.buildAtoms();
	string filename = "/tmp/buildFromCharmmTopology.pdb";
	PDBWriter writer(filename);
    writer.open();
	writer.write(sys.getAtoms());
	writer.close();

	cout << "Written pdb file " << filename << endl;
	cout << endl;

	AtomPointerVector atoms = sys.getAtoms();
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		cout << **k << endl;
	}

	return 0;
}
