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


#include "CRDReader.h"
#include "CRDWriter.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "AtomPointerVector.h"
#include "testData.h"

using namespace MSL;
using namespace std;

int main() {

	// Write a test PDB file from testData.h (/tmp/testPdb.pdb)
	
	PDBReader pRead;
	if(!pRead.read(pdbForTestCrd)) {
		cerr << "Unable to read pdbForTestCrd" << endl;
		exit(0);
	}

	AtomPointerVector atoms = pRead.getAtomPointers();
	cout << atoms << endl;

	CRDReader cRead;
	
	if(!cRead.read(testCrdFile)) {
		cerr << "Unable to read crdFile" << endl;
		exit(0);
	}

	cout << cRead.getAtomPointers();
	cRead.assignCoordinates(atoms);
		
	PDBWriter w("/tmp/testPDB.pdb");
	if(!w.open() || !w.write(atoms)) {
		cerr << "Unable to write /tmp/testPDB.pdb" << endl;
		exit(0);
	}


	CRDWriter w2("/tmp/testCRD.crd");
	if(!w2.open()) {
		cerr << "Unable to open /tmp/testCRD.pdb" << endl;
		exit(0);
	}

	w2.addRemark("This is a test file created by testCRDIO");
	if(!w2.write(atoms)) {
		cerr << "Unable to write /tmp/testCRD.pdb" << endl;
		exit(0);
	}


	return 0;

};



