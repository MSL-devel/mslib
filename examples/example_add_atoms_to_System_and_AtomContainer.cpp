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
#include <cstdlib>

#include "AtomContainer.h"
#include "System.h"
#include "MslTools.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to read and write PDBs with the
 *  AtomContainer
 *******************************************************************/

int main(int argc, char *argv[]) {

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to add atoms to the System and the AtomContainer (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;




	cout << "Add a tri-peptide ALA, ILE, ALA (19 atoms) to the AtomContainer with the addAtom function" << endl;
	AtomContainer container;
	container.addAtom("A,1,ALA,N",    2.143,  1.328,  0.000, "N");
	container.addAtom("A,1,ALA,CA",   1.539,  0.000,  0.000, "C");
	container.addAtom("A,1,ALA,CB",   2.095, -0.791,  1.207, "C");
	container.addAtom("A,1,ALA,C",    0.000,  0.000,  0.000, "C");
	container.addAtom("A,1,ALA,O",   -0.661,  1.037,  0.000, "O");
	container.addAtom("A,2,ILE,N",   -0.612, -1.210, -0.000, "N");
	container.addAtom("A,2,ILE,CA",  -2.052, -1.462, -0.000, "C");
	container.addAtom("A,2,ILE,CB",  -2.790, -0.764, -1.194, "C");
	container.addAtom("A,2,ILE,CG2", -2.141, -1.025, -2.572, "C");
	container.addAtom("A,2,ILE,CG1", -4.319, -1.018, -1.213, "C");
	container.addAtom("A,2,ILE,CD1", -5.054, -0.338, -2.380, "C");
	container.addAtom("A,2,ILE,C",   -2.221, -2.971, -0.000, "C");
	container.addAtom("A,2,ILE,O",   -1.239, -3.712, -0.000, "O");
	container.addAtom("A,3,ALA,N",   -3.474, -3.466, -0.000, "N");
	container.addAtom("A,3,ALA,CA",  -3.791, -4.877, -0.000, "C");
	container.addAtom("A,3,ALA,CB",  -3.085, -5.538,  1.207, "C");
	container.addAtom("A,3,ALA,C",   -5.297, -5.192, -0.000, "C");
	container.addAtom("A,3,ALA,OT1", -6.104, -4.223,  0.000, "O");
	container.addAtom("A,3,ALA,OT2", -5.649, -6.401, -0.000, "O");

	cout << endl;

	// all atoms can be conveniently be printed using the << operator
	cout << "Print the AtomContainer (note, it prints all atoms)" << endl;
	cout << container << endl;
	cout << endl;
	cout << "=============================" << endl;

	// cycle among atoms with the [] operator
	cout << "Cycle over the atoms of the AtomContainer" << endl;
	for (unsigned int i=0; i<container.size(); i++) {
		cout << "Atom " << i << " is " << container[i] << endl;
	}
	cout << endl;
	cout << "=============================" << endl;

	// write the coordinates to PDB
	cout << "Write the coordinates to file: /tmp/example0001_out.pdb" << endl;
	if (!container.writePdb("/tmp/example0001_out.pdb")) {
		cout << "Error writing pdb file" << endl;
	} else {
		cout << "OK" << endl;
	}

	cout << endl;
	cout << "=============================" << endl;

	cout << "Do the same thing with the System (note, this is not the most efficient way)" << endl;
	System sys;
	sys.addAtom("A,1,ALA,N",    2.143,  1.328,  0.000, "N");
	sys.addAtom("A,1,ALA,CA",   1.539,  0.000,  0.000, "C");
	sys.addAtom("A,1,ALA,CB",   2.095, -0.791,  1.207, "C");
	sys.addAtom("A,1,ALA,C",    0.000,  0.000,  0.000, "C");
	sys.addAtom("A,1,ALA,O",   -0.661,  1.037,  0.000, "O");
	sys.addAtom("A,2,ILE,N",   -0.612, -1.210, -0.000, "N");
	sys.addAtom("A,2,ILE,CA",  -2.052, -1.462, -0.000, "C");
	sys.addAtom("A,2,ILE,CB",  -2.790, -0.764, -1.194, "C");
	sys.addAtom("A,2,ILE,CG2", -2.141, -1.025, -2.572, "C");
	sys.addAtom("A,2,ILE,CG1", -4.319, -1.018, -1.213, "C");
	sys.addAtom("A,2,ILE,CD1", -5.054, -0.338, -2.380, "C");
	sys.addAtom("A,2,ILE,C",   -2.221, -2.971, -0.000, "C");
	sys.addAtom("A,2,ILE,O",   -1.239, -3.712, -0.000, "O");
	sys.addAtom("A,3,ALA,N",   -3.474, -3.466, -0.000, "N");
	sys.addAtom("A,3,ALA,CA",  -3.791, -4.877, -0.000, "C");
	sys.addAtom("A,3,ALA,CB",  -3.085, -5.538,  1.207, "C");
	sys.addAtom("A,3,ALA,C",   -5.297, -5.192, -0.000, "C");
	sys.addAtom("A,3,ALA,OT1", -6.104, -4.223,  0.000, "O");
	sys.addAtom("A,3,ALA,OT2", -5.649, -6.401, -0.000, "O");

	cout << endl;

	// all atoms can be conveniently be printed using the << operator
	cout << "Print the System (note, it prints the amino acid sequence)" << endl;
	cout << sys << endl;
	cout << endl;
	cout << "=============================" << endl;

	cout << sys.getAtomPointers() << endl;

	// cycle among atoms with the [] operator
	cout << "Cycle over the atoms of the System" << endl;
	for (unsigned int i=0; i<sys.atomSize(); i++) {
		cout << "Atom " << i << " is " << sys[i] << endl;
	}
	cout << endl;
	cout << "=============================" << endl;

	// write the coordinates to PDB
	cout << "Write the coordinates to file: /tmp/example0002_out.pdb" << endl;
	if (!sys.writePdb("/tmp/example0002_out.pdb")) {
		cout << "Error writing pdb file" << endl;
	} else {
		cout << "OK" << endl;
	}


	cout << endl;
	cout << "=============================" << endl;

	// if all atoms are passed at once (as an AtomPointerVector) creating the molecule in
	// the System will have significant less overhead.
	cout << "Add atoms to the System more efficiently by passing an AtomVpointerVector with all atoms at once" << endl;
	AtomPointerVector atoms;
	atoms.push_back(new Atom("A,1,ALA,N",    2.143,  1.328,  0.000, "N"));
	atoms.push_back(new Atom("A,1,ALA,CA",   1.539,  0.000,  0.000, "C"));
	atoms.push_back(new Atom("A,1,ALA,CB",   2.095, -0.791,  1.207, "C"));
	atoms.push_back(new Atom("A,1,ALA,C",    0.000,  0.000,  0.000, "C"));
	atoms.push_back(new Atom("A,1,ALA,O",   -0.661,  1.037,  0.000, "O"));
	atoms.push_back(new Atom("A,2,ILE,N",   -0.612, -1.210, -0.000, "N"));
	atoms.push_back(new Atom("A,2,ILE,CA",  -2.052, -1.462, -0.000, "C"));
	atoms.push_back(new Atom("A,2,ILE,CB",  -2.790, -0.764, -1.194, "C"));
	atoms.push_back(new Atom("A,2,ILE,CG2", -2.141, -1.025, -2.572, "C"));
	atoms.push_back(new Atom("A,2,ILE,CG1", -4.319, -1.018, -1.213, "C"));
	atoms.push_back(new Atom("A,2,ILE,CD1", -5.054, -0.338, -2.380, "C"));
	atoms.push_back(new Atom("A,2,ILE,C",   -2.221, -2.971, -0.000, "C"));
	atoms.push_back(new Atom("A,2,ILE,O",   -1.239, -3.712, -0.000, "O"));
	atoms.push_back(new Atom("A,3,ALA,N",   -3.474, -3.466, -0.000, "N"));
	atoms.push_back(new Atom("A,3,ALA,CA",  -3.791, -4.877, -0.000, "C"));
	atoms.push_back(new Atom("A,3,ALA,CB",  -3.085, -5.538,  1.207, "C"));
	atoms.push_back(new Atom("A,3,ALA,C",   -5.297, -5.192, -0.000, "C"));
	atoms.push_back(new Atom("A,3,ALA,OT1", -6.104, -4.223,  0.000, "O"));
	atoms.push_back(new Atom("A,3,ALA,OT2", -5.649, -6.401, -0.000, "O"));
	System sys2;
	sys2.addAtoms(atoms);

	// garbage collection: delete all pointers
	for (unsigned int i=0; i<atoms.size(); i++) {
		delete atoms[i];
	}
	atoms.clear();

	cout << endl;

	// all atoms can be conveniently be printed using the << operator
	cout << "Print the System (note, it prints the amino acid sequence)" << endl;
	cout << sys2 << endl;
	cout << endl;
	cout << "=============================" << endl;

	cout << sys2.getAtomPointers() << endl;

	// cycle among atoms with the [] operator
	cout << "Cycle over the atoms of the System" << endl;
	for (unsigned int i=0; i<sys2.atomSize(); i++) {
		cout << "Atom " << i << " is " << sys2[i] << endl;
	}
	cout << endl;
	cout << "=============================" << endl;

	// write the coordinates to PDB
	cout << "Write the coordinates to file: /tmp/example0003_out.pdb" << endl;
	if (!sys2.writePdb("/tmp/example0003_out.pdb")) {
		cout << "Error writing pdb file" << endl;
	} else {
		cout << "OK" << endl;
	}


	return 0;
}
