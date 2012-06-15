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
#include <string>

#include "AtomSelection.h"
#include "IcEntry.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "System.h"

using namespace MSL;
using namespace std;


/*************************************************************
 *  
 *    THIS TEST PROGRAM NEEDS TO BE DEVELOPED!!!
 *
 *
 *************************************************************/

int main(){


	string pdbtext = "\
ATOM      1  N   ALA A   1      -0.557   1.474 -12.560  1.00  0.00              \n\
ATOM      2  CA  ALA A   1      -1.686   1.490 -11.633  1.00  0.00              \n\
ATOM      3  CB  ALA A   1      -3.019   1.370 -12.406  1.00  0.00              \n\
ATOM      4  C   ALA A   1      -1.553   0.374 -10.582  1.00  0.00              \n\
ATOM      5  O   ALA A   1      -0.603  -0.407 -10.563  1.00  0.00              \n\
ATOM      6  N   ALA A   2      -2.542   0.281  -9.659  1.00  0.00              \n\
ATOM      7  CA  ALA A   2      -2.634  -0.689  -8.571  1.00  0.00              \n\
ATOM      8  CB  ALA A   2      -1.341  -0.604  -7.727  1.00  0.00              \n\
ATOM      9  C   ALA A   2      -3.884  -0.536  -7.686  1.00  0.00              \n\
TER      10      ALA A   2                                                      \n\
END                                                                             \n";

	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Write a pdb file and read it into an atom vector     |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	
	// write the input test PDB file
	ofstream pdb_fs;
	pdb_fs.open("/tmp/testPdb.pdb");
	pdb_fs << pdbtext;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file /tmp/testPdb.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file /tmp/testPdb.pdb" << endl;
	}

	pdb_fs.close();

	PDBReader rAv;
	//rAv.open(argv[1]);
	rAv.open("/tmp/testPdb.pdb");
	rAv.read();
	AtomPointerVector av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();
	System sys(av);
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}


	/********************************************************
	 *  Create a moledule, a ALA ALA dipeptide, param 19 atoms
	 *
	 *   RESI ALA          0.00
	 *   GROUP   
	 *   ATOM N    NH1    -0.47  !     |
	 *   ATOM HN   H       0.31  !  HN-N
	 *   ATOM CA   CT1     0.07  !     |     HB1
	 *   ATOM HA   HB      0.09  !     |    /
	 *   GROUP                   !  HA-CA--CB-HB2
	 *   ATOM CB   CT3    -0.27  !     |    \
	 *   ATOM HB1  HA      0.09  !     |     HB3
	 *   ATOM HB2  HA      0.09  !   O=C
	 *   ATOM HB3  HA      0.09  !     |
	 *   GROUP                   !
	 *   ATOM C    C       0.51
	 *   ATOM O    O      -0.51
	 *   BOND CB CA  N  HN  N  CA  
	 *   BOND C  CA  C  +N  CA HA  CB HB1  CB HB2  CB HB3 
	 *   DOUBLE O  C 
	 *   IMPR N -C CA HN  C CA +N O   
	 *   DONOR HN N   
	 *   ACCEPTOR O C   
	 *   IC -C   CA   *N   HN    1.3551 126.4900  180.0000 115.4200  0.9996
	 *   IC -C   N    CA   C     1.3551 126.4900  180.0000 114.4400  1.5390
	 *   IC N    CA   C    +N    1.4592 114.4400  180.0000 116.8400  1.3558
	 *   IC +N   CA   *C   O     1.3558 116.8400  180.0000 122.5200  1.2297
	 *   IC CA   C    +N   +CA   1.5390 116.8400  180.0000 126.7700  1.4613
	 *   IC N    C    *CA  CB    1.4592 114.4400  123.2300 111.0900  1.5461
	 *   IC N    C    *CA  HA    1.4592 114.4400 -120.4500 106.3900  1.0840
	 *   IC C    CA   CB   HB1   1.5390 111.0900  177.2500 109.6000  1.1109
	 *   IC HB1  CA   *CB  HB2   1.1109 109.6000  119.1300 111.0500  1.1119
	 *   IC HB1  CA   *CB  HB3   1.1109 109.6000 -119.5800 111.6100  1.1114
	 *
	 ********************************************************/

	sys.addIcEntry("A 0 C",  "A 1 N",  "A 1 CA", "A 1 C",  1.3551, 126.4900, 180.0000, 114.4400, 1.5390);
	sys.addIcEntry("A 1 N",  "A 1 CA", "A 1 C",  "A 2 N",  1.4592, 114.4400, 180.0000, 116.8400, 1.3558);
	sys.addIcEntry("A 2 N",  "A 1 CA", "A 1 C",  "A 1 O",  1.3558, 116.8400, 180.0000, 122.5200, 1.2297, true);
	sys.addIcEntry("A 1 CA", "A 1 C",  "A 2 N",  "A 2 CA", 1.5390, 116.8400, 180.0000, 126.7700, 1.4613);
	sys.addIcEntry("A 1 N",  "A 1 C",  "A 1 CA", "A 1 CB", 1.4592, 114.4400, 123.2300, 111.0900, 1.5461, true);
	sys.addIcEntry("A 1 C",  "A 2 N",  "A 2 CA", "A 2 C",  1.3551, 126.4900, 180.0000, 114.4400, 1.5390);
	sys.addIcEntry("A 2 N",  "A 2 CA", "A 2 C",  "A 3 N",  1.4592, 114.4400, 180.0000, 116.8400, 1.3558);
	sys.addIcEntry("A 2 N",  "A 2 C",  "A 2 CA", "A 2 CB", 1.4592, 114.4400, 123.2300, 111.0900, 1.5461, true);

	sys.fillIcFromCoor();
	sys.printIcTable();
	sys.wipeAllCoordinates();

	cout << "Test seeding with A 1 C, A 1 CA, A 1 N" << endl;
	sys.seed("A 1 C", "A 1 CA", "A 1 N");
	sys.buildAllAtoms();
	for (AtomPointerVector::iterator k = sys.getAtomPointers().begin(); k != sys.getAtomPointers().end() ; k++){
		cout << *(*k) << endl;
	}
	string filename = "/tmp/builtAtoms-A_1_C-A_1_CA-A_1_N.pdb";
	PDBWriter writer(filename);
    writer.open();
	writer.write(sys.getAtomPointers());
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;

	System sys2(sys);
}
