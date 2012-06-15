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

using namespace MSL;
using namespace std;


int main() {


	string pdbtext = "\
ATOM      1 1H   ALA A   1      -0.079  -2.341  -3.811  1.00  0.00           H  \n\
ATOM      2 2H   ALA A   1       0.350  -0.712  -3.601  1.00  0.00           H  \n\
ATOM      3 3H   ALA A   1      -0.809  -1.422  -2.584  1.00  0.00           H  \n\
ATOM      4  N   ALA A   1       0.072  -1.586  -3.112  1.00  0.00           N  \n\
ATOM      5  CA  ALA A   1       1.121  -1.981  -2.192  1.00  0.00           C  \n\
ATOM      6  HA  ALA A   1       0.834  -2.871  -1.644  1.00  0.00           H  \n\
ATOM      7  CB  ALA A   1       2.417  -2.319  -2.965  1.00  0.00           C  \n\
ATOM      8 1HB  ALA A   1       2.218  -3.144  -3.682  1.00  0.00           H  \n\
ATOM      9 2HB  ALA A   1       2.775  -1.441  -3.546  1.00  0.00           H  \n\
ATOM     10 3HB  ALA A   1       3.228  -2.646  -2.281  1.00  0.00           H  \n\
ATOM     11  C   ALA A   1       1.370  -0.901  -1.154  1.00  0.00           C  \n\
ATOM     12  O   ALA A   1       1.549  -1.196   0.026  1.00  0.00           O  \n\
ATOM     13  N   ALA A   2       1.384   0.372  -1.594  1.00  0.00           N  \n\
ATOM     14  H   ALA A   2       1.256   0.631  -2.551  1.00  0.00           H  \n\
ATOM     15  CA  ALA A   2       1.592   1.470  -0.676  1.00  0.00           C  \n\
ATOM     16  HA  ALA A   2       2.518   1.340  -0.128  1.00  0.00           H  \n\
ATOM     17  CB  ALA A   2       1.703   2.804  -1.449  1.00  0.00           C  \n\
ATOM     18 1HB  ALA A   2       2.550   2.751  -2.166  1.00  0.00           H  \n\
ATOM     19 2HB  ALA A   2       0.776   3.007  -2.030  1.00  0.00           H  \n\
ATOM     20 3HB  ALA A   2       1.886   3.660  -0.764  1.00  0.00           H  \n\
ATOM     21  C   ALA A   2       0.485   1.530   0.362  1.00  0.00           C  \n\
ATOM     22  O   ALA A   2       0.745   1.757   1.542  1.00  0.00           O  \n\
ATOM     23  N   LEU A   2       1.384   0.372  -1.594  1.00  0.00           N  \n\
ATOM     24  H   LEU A   2       1.258   0.650  -2.544  1.00  0.00           H  \n\
ATOM     25  CA  LEU A   2       1.592   1.470  -0.676  1.00  0.00           C  \n\
ATOM     26  HA  LEU A   2       2.517   1.270  -0.155  1.00  0.00           H  \n\
ATOM     27  CB  LEU A   2       1.660   2.816  -1.442  1.00  0.00           C  \n\
ATOM     28 1HB  LEU A   2       0.723   2.906  -2.030  1.00  0.00           H  \n\
ATOM     29 2HB  LEU A   2       1.645   3.649  -0.710  1.00  0.00           H  \n\
ATOM     30  CG  LEU A   2       2.861   2.976  -2.409  1.00  0.00           C  \n\
ATOM     31  HG  LEU A   2       2.925   2.109  -3.115  1.00  0.00           H  \n\
ATOM     32  CD1 LEU A   2       2.671   4.240  -3.268  1.00  0.00           C  \n\
ATOM     33 1HD1 LEU A   2       2.607   5.142  -2.625  1.00  0.00           H  \n\
ATOM     34 2HD1 LEU A   2       3.524   4.369  -3.966  1.00  0.00           H  \n\
ATOM     35 3HD1 LEU A   2       1.738   4.172  -3.864  1.00  0.00           H  \n\
ATOM     36  CD2 LEU A   2       4.184   3.048  -1.624  1.00  0.00           C  \n\
ATOM     37 1HD2 LEU A   2       4.332   2.121  -1.031  1.00  0.00           H  \n\
ATOM     38 2HD2 LEU A   2       5.057   3.158  -2.300  1.00  0.00           H  \n\
ATOM     39 3HD2 LEU A   2       4.174   3.910  -0.924  1.00  0.00           H  \n\
ATOM     40  C   LEU A   2       0.485   1.530   0.362  1.00  0.00           C  \n\
ATOM     41  O   LEU A   2       0.745   1.757   1.542  1.00  0.00           O  \n\
ATOM     42  N   ALA A   3      -0.771   1.325  -0.078  1.00  0.00           N  \n\
ATOM     43  H   ALA A   3      -1.005   1.154  -1.035  1.00  0.00           H  \n\
ATOM     44  CA  ALA A   3      -1.889   1.341   0.841  1.00  0.00           C  \n\
ATOM     45  HA  ALA A   3      -1.919   2.276   1.389  1.00  0.00           H  \n\
ATOM     46  CB  ALA A   3      -3.223   1.221   0.068  1.00  0.00           C  \n\
ATOM     47 1HB  ALA A   3      -3.315   2.066  -0.648  1.00  0.00           H  \n\
ATOM     48 2HB  ALA A   3      -3.264   0.274  -0.513  1.00  0.00           H  \n\
ATOM     49 3HB  ALA A   3      -4.097   1.255   0.753  1.00  0.00           H  \n\
ATOM     50  C   ALA A   3      -1.758   0.240   1.878  1.00  0.00           C  \n\
ATOM     51  O   ALA A   3      -2.026   0.457   3.058  1.00  0.00           O  \n\
ATOM     52  N   ALA A   4      -1.339  -0.962   1.438  1.00  0.00           N  \n\
ATOM     53  H   ALA A   4      -1.130  -1.163   0.481  1.00  0.00           H  \n\
ATOM     54  CA  ALA A   4      -1.163  -2.066   2.357  1.00  0.00           C  \n\
ATOM     55  HA  ALA A   4      -2.078  -2.257   2.905  1.00  0.00           H  \n\
ATOM     56  CB  ALA A   4      -0.816  -3.359   1.585  1.00  0.00           C  \n\
ATOM     57 1HB  ALA A   4      -1.631  -3.596   0.869  1.00  0.00           H  \n\
ATOM     58 2HB  ALA A   4       0.125  -3.237   1.004  1.00  0.00           H  \n\
ATOM     59 3HB  ALA A   4      -0.698  -4.226   2.270  1.00  0.00           H  \n\
ATOM     60  C   ALA A   4      -0.101  -1.747   3.394  1.00  0.00           C  \n\
ATOM     61  OXT ALA A   4       0.467  -0.624   3.333  1.00  0.00           O  \n\
ATOM     62  O   ALA A   4       0.159  -2.622   4.263  1.00  0.00           O  \n\
TER      63                                                                     \n\
ATOM     64 1H   ALA B   1       9.921  -2.341  -3.811  1.00  0.00           H  \n\
ATOM     65 2H   ALA B   1      10.350  -0.712  -3.601  1.00  0.00           H  \n\
ATOM     66 3H   ALA B   1       9.191  -1.422  -2.584  1.00  0.00           H  \n\
ATOM     67  N   ALA B   1      10.072  -1.586  -3.112  1.00  0.00           N  \n\
ATOM     68  CA  ALA B   1      11.121  -1.981  -2.192  1.00  0.00           C  \n\
ATOM     69  HA  ALA B   1      10.834  -2.871  -1.644  1.00  0.00           H  \n\
ATOM     70  CB  ALA B   1      12.417  -2.319  -2.965  1.00  0.00           C  \n\
ATOM     71 1HB  ALA B   1      12.218  -3.144  -3.682  1.00  0.00           H  \n\
ATOM     72 2HB  ALA B   1      12.775  -1.441  -3.546  1.00  0.00           H  \n\
ATOM     73 3HB  ALA B   1      13.228  -2.646  -2.281  1.00  0.00           H  \n\
ATOM     74  C   ALA B   1      11.370  -0.901  -1.154  1.00  0.00           C  \n\
ATOM     75  O   ALA B   1      11.549  -1.196   0.026  1.00  0.00           O  \n\
ATOM     76  N   ALA B   2      11.384   0.372  -1.594  1.00  0.00           N  \n\
ATOM     77  H   ALA B   2      11.256   0.631  -2.551  1.00  0.00           H  \n\
ATOM     78  CA  ALA B   2      11.592   1.470  -0.676  1.00  0.00           C  \n\
ATOM     79  HA  ALA B   2      12.518   1.340  -0.128  1.00  0.00           H  \n\
ATOM     80  CB  ALA B   2      11.703   2.804  -1.449  1.00  0.00           C  \n\
ATOM     81 1HB  ALA B   2      12.550   2.751  -2.166  1.00  0.00           H  \n\
ATOM     82 2HB  ALA B   2      10.776   3.007  -2.030  1.00  0.00           H  \n\
ATOM     83 3HB  ALA B   2      11.886   3.660  -0.764  1.00  0.00           H  \n\
ATOM     84  C   ALA B   2      10.485   1.530   0.362  1.00  0.00           C  \n\
ATOM     85  O   ALA B   2      10.745   1.757   1.542  1.00  0.00           O  \n\
ATOM     86  N  ALEU B   3       9.229   1.325  -0.078  1.00  0.00           N  \n\
ATOM     87  H  ALEU B   3       8.976   1.153  -1.028  1.00  0.00           H  \n\
ATOM     88  CA ALEU B   3       8.111   1.341   0.841  1.00  0.00           C  \n\
ATOM     89  HA ALEU B   3       8.149   2.286   1.362  1.00  0.00           H  \n\
ATOM     90  CB ALEU B   3       6.773   1.178   0.075  1.00  0.00           C  \n\
ATOM     91 1HB ALEU B   3       6.845   0.239  -0.513  1.00  0.00           H  \n\
ATOM     92 2HB ALEU B   3       5.955   1.019   0.808  1.00  0.00           H  \n\
ATOM     93  CG ALEU B   3       6.409   2.334  -0.891  1.00  0.00           C  \n\
ATOM     94  HG ALEU B   3       7.252   2.546  -1.597  1.00  0.00           H  \n\
ATOM     95  CD1ALEU B   3       5.197   1.931  -1.750  1.00  0.00           C  \n\
ATOM     96 1HD1ALEU B   3       4.319   1.712  -1.107  1.00  0.00           H  \n\
ATOM     97 2HD1ALEU B   3       4.922   2.749  -2.448  1.00  0.00           H  \n\
ATOM     98 3HD1ALEU B   3       5.424   1.023  -2.347  1.00  0.00           H  \n\
ATOM     99  CD2ALEU B   3       6.111   3.624  -0.105  1.00  0.00           C  \n\
ATOM    100 1HD2ALEU B   3       6.999   3.928   0.488  1.00  0.00           H  \n\
ATOM    101 2HD2ALEU B   3       5.853   4.466  -0.781  1.00  0.00           H  \n\
ATOM    102 3HD2ALEU B   3       5.264   3.466   0.595  1.00  0.00           H  \n\
ATOM    103  C  ALEU B   3       8.242   0.240   1.878  1.00  0.00           C  \n\
ATOM    104  O  ALEU B   3       7.974   0.457   3.058  1.00  0.00           O  \n\
ATOM    105  N  BLEU B   3       9.229   1.325  -0.078  1.00  0.00           N  \n\
ATOM    106  H  BLEU B   3       8.976   1.153  -1.028  1.00  0.00           H  \n\
ATOM    107  CA BLEU B   3       8.111   1.341   0.841  1.00  0.00           C  \n\
ATOM    108  HA BLEU B   3       8.149   2.286   1.362  1.00  0.00           H  \n\
ATOM    109  CB BLEU B   3       6.800   1.212   0.024  1.00  0.00           C  \n\
ATOM    110 1HB BLEU B   3       6.756   2.083  -0.663  1.00  0.00           H  \n\
ATOM    111 2HB BLEU B   3       6.866   0.313  -0.624  1.00  0.00           H  \n\
ATOM    112  CG BLEU B   3       5.487   1.176   0.845  1.00  0.00           C  \n\
ATOM    113  HG BLEU B   3       5.510   0.354   1.606  1.00  0.00           H  \n\
ATOM    114  CD1BLEU B   3       5.286   2.494   1.616  1.00  0.00           C  \n\
ATOM    115 1HD1BLEU B   3       5.233   3.352   0.914  1.00  0.00           H  \n\
ATOM    116 2HD1BLEU B   3       4.344   2.464   2.203  1.00  0.00           H  \n\
ATOM    117 3HD1BLEU B   3       6.130   2.671   2.316  1.00  0.00           H  \n\
ATOM    118  CD2BLEU B   3       4.297   0.918  -0.098  1.00  0.00           C  \n\
ATOM    119 1HD2BLEU B   3       4.424  -0.053  -0.621  1.00  0.00           H  \n\
ATOM    120 2HD2BLEU B   3       3.335   0.879   0.455  1.00  0.00           H  \n\
ATOM    121 3HD2BLEU B   3       4.228   1.718  -0.864  1.00  0.00           H  \n\
ATOM    122  C  BLEU B   3       8.242   0.240   1.878  1.00  0.00           C  \n\
ATOM    123  O  BLEU B   3       7.974   0.457   3.058  1.00  0.00           O  \n\
ATOM    124  N  CLEU B   3       9.229   1.325  -0.078  1.00  0.00           N  \n\
ATOM    125  H  CLEU B   3       8.976   1.153  -1.028  1.00  0.00           H  \n\
ATOM    126  CA CLEU B   3       8.111   1.341   0.841  1.00  0.00           C  \n\
ATOM    127  HA CLEU B   3       8.149   2.286   1.362  1.00  0.00           H  \n\
ATOM    128  CB CLEU B   3       6.771   1.179   0.079  1.00  0.00           C  \n\
ATOM    129 1HB CLEU B   3       6.765   0.160  -0.360  1.00  0.00           H  \n\
ATOM    130 2HB CLEU B   3       5.935   1.204   0.809  1.00  0.00           H  \n\
ATOM    131  CG CLEU B   3       6.518   2.202  -1.058  1.00  0.00           C  \n\
ATOM    132  HG CLEU B   3       7.319   2.136  -1.837  1.00  0.00           H  \n\
ATOM    133  CD1CLEU B   3       5.185   1.882  -1.759  1.00  0.00           C  \n\
ATOM    134 1HD1CLEU B   3       4.343   1.934  -1.038  1.00  0.00           H  \n\
ATOM    135 2HD1CLEU B   3       4.989   2.608  -2.576  1.00  0.00           H  \n\
ATOM    136 3HD1CLEU B   3       5.206   0.861  -2.195  1.00  0.00           H  \n\
ATOM    137  CD2CLEU B   3       6.508   3.640  -0.506  1.00  0.00           C  \n\
ATOM    138 1HD2CLEU B   3       7.484   3.877  -0.032  1.00  0.00           H  \n\
ATOM    139 2HD2CLEU B   3       6.333   4.386  -1.308  1.00  0.00           H  \n\
ATOM    140 3HD2CLEU B   3       5.713   3.757   0.260  1.00  0.00           H  \n\
ATOM    141  C  CLEU B   3       8.242   0.240   1.878  1.00  0.00           C  \n\
ATOM    142  O  CLEU B   3       7.974   0.457   3.058  1.00  0.00           O  \n\
ATOM    143  N   ALA B   4       8.661  -0.962   1.438  1.00  0.00           N  \n\
ATOM    144  H   ALA B   4       8.870  -1.163   0.481  1.00  0.00           H  \n\
ATOM    145  CA  ALA B   4       8.837  -2.066   2.357  1.00  0.00           C  \n\
ATOM    146  HA  ALA B   4       7.922  -2.257   2.905  1.00  0.00           H  \n\
ATOM    147  CB  ALA B   4       9.184  -3.359   1.585  1.00  0.00           C  \n\
ATOM    148 1HB  ALA B   4       8.369  -3.596   0.869  1.00  0.00           H  \n\
ATOM    149 2HB  ALA B   4      10.125  -3.237   1.004  1.00  0.00           H  \n\
ATOM    150 3HB  ALA B   4       9.302  -4.226   2.270  1.00  0.00           H  \n\
ATOM    151  C   ALA B   4       9.899  -1.747   3.394  1.00  0.00           C  \n\
ATOM    152  OXT ALA B   4      10.467  -0.624   3.333  1.00  0.00           O  \n\
ATOM    153  O   ALA B   4      10.159  -2.622   4.263  1.00  0.00           O  \n\
TER     154                                                                     \n\
END                                                                             \n";

	
	// write the input test PDB file
	ofstream pdb_fs;
	pdb_fs.open("testPdb.pdb");
	pdb_fs << pdbtext;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file testPdb.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file testPdb.pdb" << endl;
	}

	pdb_fs.close();

	cout << "\n=== Test Read testPdb.pdb into AtomPointerVector ===\n\n";

	AtomPointerVector av;
	PDBReader rAv;
	//rAv.open(argv[1]);
	rAv.open("testPdb.pdb");
	rAv.read();
	av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();

	//for (AtomPointerVector::iterator itAv = av.begin(); itAv != av.end() ; itAv++){
	//	cout << (*itAv)->toString()<<endl;
	//}

	cout << "Create a system from the atom vector" << endl;
	
	System sys(av);
	cout << "The system has " << sys.atomSize() << " atoms" << endl;

	for (AtomPointerVector::iterator k=sys.getAtomPointers().begin(); k!= sys.getAtomPointers().end(); k++) {
		cout << **k << endl;
	}

	cout << endl;
	cout << "==========================" << endl;
	cout << "Print the positions using an iterator over the getPositions() function" << endl;
	
	for (vector<Position*>::iterator k=sys.getPositions().begin(); k!=sys.getPositions().end(); k++) {
		cout << (*k)->getChainId() << " " << (*k)->getResidueNumber() << (*k)->getResidueIcode() << " " << (*k)->getResidueName() << endl;
	}

	cout << endl;
	cout << "==========================" << endl;
	cout << "Print the positions using an indes over the getPosition(size_t _n) function" << endl;
	
	for (unsigned int i=0; i<sys.positionSize(); i++) {
		cout << sys.getPosition(i).getChainId() << " " << sys.getPosition(i).getResidueNumber() << sys.getPosition(i).getResidueIcode() << " " << sys.getPosition(i).getResidueName() << endl;
	}

	cout << endl;
	cout << "==========================" << endl;
	cout << "Print the residues using an indes over the getResidue(size_t _n) function" << endl;
	
	for (unsigned int i=0; i<sys.positionSize(); i++) {
		cout << sys.getResidue(i).getChainId() << " " << sys.getResidue(i).getResidueNumber() << sys.getResidue(i).getResidueIcode() << " " << sys.getResidue(i).getResidueName() << endl;
	}

	return 0;

};


