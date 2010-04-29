/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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


#include "CRDReader.h"
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
	w.open();
	if(!w.write(atoms)) {
		cerr << "Unable to write /tmp/testPDB.pdb" << endl;
		exit(0);
	}


	return 0;

};



