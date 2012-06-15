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


#include "PDBReader.h"
#include "PDBWriter.h"
#include "CartesianPoint.h"
#include "AtomSelection.h"
#include "System.h"
#include "testData.h"

using namespace MSL;
using namespace std;

int main() {

	// Write a test PDB file from testData.h (/tmp/testPdb.pdb)
	writePdbFile();


	/*
	cout << "\n=== Test Read /tmp/testPdb.pdb into AtomPointerVector ===\n\n";

	AtomPointerVector av;
	PDBReader rAv;
	//rAv.open(argv[1]);
	rAv.open("/tmp/testPdb.pdb");
	rAv.read();
	av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();

	for (AtomPointerVector::iterator itAv = av.begin(); itAv != av.end() ; itAv++){
		cout << (*itAv)->toString()<<endl;
	}

	cout << "Make selection for residue 1, create a residue from the subset atom vector and print the atoms from the residue" << endl;
	AtomSelection sel(av);
	
	AtomPointerVector subset = sel.select("res1, resi 1");
	cout << "Selected residue 1, atom vector with size " << subset.size() << endl;
	for (AtomPointerVector::iterator itAv = subset.begin(); itAv != subset.end() ; itAv++){
		cout << (*itAv)->toString()<<endl;
	}

	
	cout << "Print from the residue:" << endl;
	Residue Ala1(subset, "ALA", 1, "");

	for (AtomPointerVector::iterator k=Ala1.getAtomPointers().begin(); k!= Ala1.getAtomPointers().end(); k++) {
		cout << **k << endl;
	}

	// ======================================================

	cout << "Make selection for chain A, create a chain from the subset atom vector and print the atoms from the chain" << endl;
	
	subset = sel.select("chainA, chain A");
	cout << "Selected chain A, atom vector with size " << subset.size() << endl;
	for (AtomPointerVector::iterator itAv = subset.begin(); itAv != subset.end() ; itAv++){
		cout << (*itAv)->toString()<<endl;
	}

	
	cout << "Print from the chain:" << endl;
	Chain chainA(subset, "A");

	for (AtomPointerVector::iterator k=chainA.getAtomPointers().begin(); k!= chainA.getAtomPointers().end(); k++) {
		cout << **k << endl;
	}

	// ======================================================

	cout << "Create a system from the atom vector" << endl;
	
	System sys(av);
	cout << "The system has " << sys.atomSize() << " atoms" << endl;

	for (AtomPointerVector::iterator k=sys.getAtomPointers().begin(); k!= sys.getAtomPointers().end(); k++) {
		cout << **k << endl;
	}

	
	cout << " * Does chain A exists? ";
	if (sys.exists("A")) {
		cout << "YES.  Get the chain..." << endl;
		Chain & found = sys.getLastFoundChain();
		cout << "ChainID of found chain = " << found.getChainId() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does chain B exists? ";
	if (sys.exists("B")) {
		cout << "YES.  Get the chain..." << endl;
		Chain & found = sys.getLastFoundChain();
		cout << "ChainID of found chain = " << found.getChainId() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does chain C exists? ";
	if (sys.exists("C")) {
		cout << "YES.  Get the chain..." << endl;
		Chain & found = sys.getLastFoundChain();
		cout << "ChainID of found chain = " << found.getChainId() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 4 A exists? ";
	if (sys.exists("A", 4)) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 5 A exists? ";
	if (sys.exists("A", 5)) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 4 C exists? ";
	if (sys.exists("C", 4)) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 4 A exists (using a string resNum)? ";
	if (sys.exists("A", "4")) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 4A A exists (using a string resNum)? ";
	if (sys.exists("A", "4A")) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does residue 3A B exists (using a string resNum)? ";
	if (sys.exists("B", "3A")) {
		cout << "YES.  Get the position..." << endl;
		Position & foundPos = sys.getLastFoundPosition();
		cout << "chain resnum/icode of found position = " << foundPos.getChainId() << " " << foundPos.getResidueNumber() << "/\"" << foundPos.getResidueIcode() << "\"" << endl;
		cout << "... get the residue..." << endl;
		Residue & foundRes = sys.getLastFoundResidue();
		cout << "chain resnum/icode name of found residue = " << foundRes.getChainId() << " " << foundRes.getResidueNumber() << "/\"" << foundRes.getResidueIcode() << "\" " << foundRes.getResidueName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does atom 4 A CA exists? ";
	if (sys.exists("A", 4, "CA")) {
		cout << "YES.  Get the atom..." << endl;
		Atom & found = sys.getLastFoundAtom();
		cout << "chain resnum/icode name of found atom = " << found.getChainId() << " " << found.getResidueNumber() << "/\"" << found.getResidueIcode() << "\" " << found.getName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " * Does atom 4 A CD1 exists? ";
	if (sys.exists("A", 4, "CD1")) {
		cout << "YES.  Get the atom..." << endl;
		Atom & found = sys.getLastFoundAtom();
		cout << "chain resnum/icode name of found atom = " << found.getChainId() << " " << found.getResidueNumber() << "/\"" << found.getResidueIcode() << "\" " << found.getName() << endl;
	} else {
		cout << "NO" << endl;
	}
	cout << " = = = =" << endl;

	cout << "Request and print the CB of Ala 2 in chain B" << endl;
	cout << sys("B")(2)("CA") << endl; // using int for resNum
	cout << sys("B")("2")("CA") << endl; // using string for resNum

	// Writer tests
	//
	//  TEST:
	//     Can PDBWriter write a PDB from an AtomPointerVector?
	//
	cout << "\n=== Write the PDB files from the System ===\n\n";
	
	cout << " #1 AAAA AALA with LEU B 3 in rotamer 0 (AAAA-AAL0A.pdb)" << endl;
	AtomPointerVector sysAtoms = sys.getAtomPointers();

	PDBWriter w("AAAA-AAL0A.pdb");
    w.open();
	w.write(sysAtoms);
	w.close();

	cout << " #2 Change identity for position A 2 ALA->LEU, with LEU B 3 still in rotamer 0 (ALAA-AAL0A.pdb)" << endl;
	cout << "Position 2 A has size " << sys("A").getPosition(2).size() << " and is currently in the " << sys("A")(2).getResidueName() << " state" << endl;
	if (sys("A").getPosition(2).identityExists("LEU")) {
		cout << "The identity LEU is available at position 2 A" << endl;
	} else {
		cout << "The identity LEU is NOT available at position 2 A" << endl;
		exit(1);
	}
	sys("A").getPosition(2).setActiveIdentity("LEU");
	//sys("A").getPosition(2).setActiveIdentity(1); // this would also do it
	cout << "SWITCH IDENTITY... Position 2 A is now in the " << sys("A")(2).getResidueName() << " state" << endl;
	cout << "* Write PDB file ALAA-AAL0A.pdb" << endl;
	sysAtoms = sys.getAtomPointers();

	w.open("ALAA-AAL0A.pdb");
	w.write(sysAtoms);
	w.close();

	cout << " #3 Change rotamer (alt conf) for position B 3 (ALAA-AAL1A.pdb)" << endl;
	cout << "Position 3 B has " << sys("B")(3).getNumberOfAltConformations() << " conformations ";
	cout << "(Chi 1 = " << sys("B")(3)("N").dihedral(sys("B")(3)("CA"), sys("B")(3)("CB"), sys("B")(3)("CG"));
	cout << "; Chi 2 = " << sys("B")(3)("CA").dihedral(sys("B")(3)("CB"), sys("B")(3)("CG"), sys("B")(3)("CD1")) << ")" << endl;

	sys("B")(3).setActiveConformation(1);
	cout << "SWITCH CONFORMATION TO #1... ";
	cout << "(Chi 1 = " << sys("B")(3)("N").dihedral(sys("B")(3)("CA"), sys("B")(3)("CB"), sys("B")(3)("CG"));
	cout << "; Chi 2 = " << sys("B")(3)("CA").dihedral(sys("B")(3)("CB"), sys("B")(3)("CG"), sys("B")(3)("CD1")) << ")" << endl;
	cout << " * Write PDB file ALAA-AAL1A.pdb" << endl;
	sysAtoms = sys.getAtomPointers();
	w.open("ALAA-AAL1A.pdb");
	w.write(sysAtoms);
	w.close();

	sys("B")(3).setActiveConformation(2);
	cout << "SWITCH CONFORMATION TO #2... ";
	cout << "(Chi 1 = " << sys("B")(3)("N").dihedral(sys("B")(3)("CA"), sys("B")(3)("CB"), sys("B")(3)("CG"));
	cout << "; Chi 2 = " << sys("B")(3)("CA").dihedral(sys("B")(3)("CB"), sys("B")(3)("CG"), sys("B")(3)("CD1")) << ")" << endl;
	cout << " * Write PDB file ALAA-AAL2A.pdb" << endl;
	sysAtoms = sys.getAtomPointers();
	w.open("ALAA-AAL2A.pdb");
	w.write(sysAtoms);
	w.close();



	// Read from string
	PDBReader strRead;
	strRead.read(pdbtext);
	cout << "String-read PDB: "<<endl;
	for (uint i = 0; i < strRead.getAtomPointers().size();i++){
		cout << strRead.getAtomPointers()[0];
	}
	*/
	return 0;

};


