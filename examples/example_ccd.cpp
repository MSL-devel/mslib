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

#include "System.h"
#include "MslTools.h"
#include "CCD.h"
#include "SysEnv.h"


using namespace std;
using namespace MSL;

static SysEnv SYSENV;

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_ccd <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to use CCD in MSL (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	string pdbFile = "example0008_caOnly.pdb";
	pdbFile = (string)argv[1] + "/" + pdbFile;

	// Read pdb file into a System object
	System sys;
	sys.readPdb(pdbFile);

	CCD sampleCCD;

	// Do local sampling inside CCD object
	sampleCCD.localSample(sys.getAtomPointers(),10,10);

	cout << "Create new system"<<endl;
	System newSystem;
	newSystem.addAtoms(sampleCCD.getAtomPointers());
	cout << "Write  new system"<<endl;

	newSystem.writePdb("/tmp/ccdEnsemble.pdb",true);

}
