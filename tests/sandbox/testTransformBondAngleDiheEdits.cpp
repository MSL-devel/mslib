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
#include "AtomBondBuilder.h"
#include "System.h"
#include "Transforms.h"
#include "CharmmSystemBuilder.h"


using namespace std;

using namespace MSL;


int main(int argc,char *argv[]) {

	string inputPdb = "/tmp/testEdit-before.pdb";

	string pdbContent = "ATOM      1  N   LEU A   1      32.638  63.213 139.911  1.00 26.64           N  \n\
ATOM      2  CA  LEU A   1      33.480  64.135 140.669  1.00 27.35           C  \n\
ATOM      3  C   LEU A   1      32.799  64.558 141.964  1.00 26.91           C  \n\
ATOM      4  O   LEU A   1      33.465  64.790 142.974  1.00 27.54           O  \n\
ATOM      5  CB  LEU A   1      33.811  65.389 139.856  1.00 28.81           C  \n\
ATOM      6  CG  LEU A   1      35.081  65.375 139.008  1.00 31.78           C  \n\
ATOM      7  CD1 LEU A   1      35.400  66.798 138.584  1.00 32.80           C  \n\
ATOM      8  CD2 LEU A   1      36.243  64.817 139.808  1.00 32.23           C  \n\
ATOM      9  N   ALA A   2      31.473  64.676 141.933  1.00 24.98           N  \n\
ATOM     10  CA  ALA A   2      30.735  65.061 143.128  1.00 24.58           C  \n\
ATOM     11  C   ALA A   2      30.966  63.974 144.181  1.00 24.15           C  \n\
ATOM     12  O   ALA A   2      30.871  64.230 145.387  1.00 21.91           O  \n\
ATOM     13  CB  ALA A   2      29.247  65.203 142.814  1.00 24.56           C  \n\
ATOM     14  N   GLY A   3      31.271  62.764 143.706  1.00 22.86           N  \n\
ATOM     15  CA  GLY A   3      31.545  61.643 144.591  1.00 20.65           C  \n\
ATOM     16  C   GLY A   3      32.724  61.940 145.505  1.00 22.02           C  \n\
ATOM     17  O   GLY A   3      32.975  61.224 146.476  1.00 20.63           O  \n\
ATOM     18  N   THR A   4      33.458  63.004 145.193  1.00 22.28           N  \n\
ATOM     19  CA  THR A   4      34.586  63.413 146.014  1.00 22.02           C  \n\
ATOM     20  C   THR A   4      34.078  63.768 147.407  1.00 21.29           C  \n\
ATOM     21  O   THR A   4      34.687  63.407 148.413  1.00 23.44           O  \n\
ATOM     22  CB  THR A   4      35.289  64.684 145.453  1.00 23.12           C  \n\
ATOM     23  OG1 THR A   4      35.978  64.378 144.234  1.00 21.17           O  \n\
ATOM     24  CG2 THR A   4      36.288  65.230 146.470  1.00 23.53           C  \n\
ATOM     25  N   PHE A   5      32.945  64.461 147.456  1.00 19.45           N  \n\
ATOM     26  CA  PHE A   5      32.380  64.920 148.718  1.00 19.78           C  \n\
ATOM     27  C   PHE A   5      31.454  63.979 149.458  1.00 20.88           C  \n\
ATOM     28  O   PHE A   5      31.527  63.879 150.682  1.00 22.08           O  \n\
ATOM     29  CB  PHE A   5      31.664  66.260 148.511  1.00 19.24           C  \n\
ATOM     30  CG  PHE A   5      32.508  67.288 147.812  1.00 18.36           C  \n\
ATOM     31  CD1 PHE A   5      32.450  67.429 146.426  1.00 18.98           C  \n\
ATOM     32  CD2 PHE A   5      33.407  68.068 148.530  1.00 19.77           C  \n\
ATOM     33  CE1 PHE A   5      33.278  68.330 145.755  1.00 17.70           C  \n\
ATOM     34  CE2 PHE A   5      34.248  68.977 147.876  1.00 22.07           C  \n\
ATOM     35  CZ  PHE A   5      34.181  69.107 146.479  1.00 21.53           C  \n\
ATOM     36  N   SER A   6      30.577  63.292 148.740  1.00 21.04           N  \n\
ATOM     37  CA  SER A   6      29.639  62.398 149.400  1.00 21.94           C  \n\
ATOM     38  C   SER A   6      29.407  61.107 148.628  1.00 22.78           C  \n\
ATOM     39  O   SER A   6      30.070  60.836 147.626  1.00 22.74           O  \n\
ATOM     40  CB  SER A   6      28.308  63.120 149.638  1.00 22.76           C  \n\
ATOM     41  OG  SER A   6      27.715  63.525 148.415  1.00 24.76           O  \n\
ATOM     42  N   THR A   7      28.443  60.321 149.088  1.00 21.12           N  \n\
ATOM     43  CA  THR A   7      28.174  59.031 148.472  1.00 22.27           C  \n\
ATOM     44  C   THR A   7      27.049  58.950 147.454  1.00 21.33           C  \n\
ATOM     45  O   THR A   7      26.051  59.661 147.525  1.00 21.53           O  \n\
ATOM     46  CB  THR A   7      27.917  57.972 149.552  1.00 21.27           C  \n\
ATOM     47  OG1 THR A   7      26.708  58.292 150.250  1.00 23.02           O  \n\
ATOM     48  CG2 THR A   7      29.071  57.942 150.545  1.00 22.14           C  \n\
ATOM     49  N   TYR A   8      27.244  58.052 146.498  1.00 21.60           N  \n\
ATOM     50  CA  TYR A   8      26.289  57.803 145.439  1.00 19.87           C  \n\
ATOM     51  C   TYR A   8      26.239  56.301 145.233  1.00 20.85           C  \n\
ATOM     52  O   TYR A   8      27.225  55.604 145.475  1.00 22.86           O  \n\
ATOM     53  CB  TYR A   8      26.716  58.549 144.180  1.00 18.71           C  \n\
ATOM     54  CG  TYR A   8      26.649  60.040 144.401  1.00 18.92           C  \n\
ATOM     55  CD1 TYR A   8      27.788  60.771 144.732  1.00 18.01           C  \n\
ATOM     56  CD2 TYR A   8      25.421  60.709 144.375  1.00 20.42           C  \n\
ATOM     57  CE1 TYR A   8      27.708  62.137 145.040  1.00 18.83           C  \n\
ATOM     58  CE2 TYR A   8      25.331  62.069 144.682  1.00 19.12           C  \n\
ATOM     59  CZ  TYR A   8      26.474  62.771 145.014  1.00 18.75           C  \n\
ATOM     60  OH  TYR A   8      26.383  64.102 145.345  1.00 20.73           O  \n\
ATOM     61  N   PRO A   9      25.084  55.776 144.803  1.00 19.77           N  \n\
ATOM     62  CA  PRO A   9      24.908  54.340 144.580  1.00 20.03           C  \n\
ATOM     63  C   PRO A   9      25.613  53.749 143.373  1.00 22.53           C  \n\
ATOM     64  O   PRO A   9      25.816  54.419 142.356  1.00 22.35           O  \n\
ATOM     65  CB  PRO A   9      23.397  54.204 144.463  1.00 20.99           C  \n\
ATOM     66  CG  PRO A   9      23.035  55.462 143.719  1.00 21.16           C  \n\
ATOM     67  CD  PRO A   9      23.854  56.514 144.454  1.00 20.56           C  \n\
ATOM     68  N   ASN A  10      25.983  52.479 143.503  1.00 23.54           N  \n\
ATOM     69  CA  ASN A  10      26.621  51.751 142.422  1.00 23.65           C  \n\
ATOM     70  C   ASN A  10      25.596  51.725 141.283  1.00 25.49           C  \n\
ATOM     71  O   ASN A  10      24.392  51.625 141.520  1.00 25.75           O  \n\
ATOM     72  CB  ASN A  10      26.949  50.327 142.879  1.00 21.81           C  \n\
ATOM     73  CG  ASN A  10      27.546  49.479 141.773  1.00 22.12           C  \n\
ATOM     74  OD1 ASN A  10      26.845  49.058 140.852  1.00 22.07           O  \n\
ATOM     75  ND2 ASN A  10      28.853  49.228 141.853  1.00 21.18           N  \n\
ATOM     76  N   PRO A  11      26.060  51.824 140.034  1.00 26.62           N  \n\
ATOM     77  CA  PRO A  11      25.167  51.811 138.869  1.00 28.98           C  \n\
ATOM     78  C   PRO A  11      24.186  50.633 138.816  1.00 30.06           C  \n\
ATOM     79  O   PRO A  11      23.065  50.777 138.342  1.00 31.59           O  \n\
ATOM     80  CB  PRO A  11      26.140  51.792 137.693  1.00 28.46           C  \n\
ATOM     81  CG  PRO A  11      27.335  52.515 138.237  1.00 27.24           C  \n\
ATOM     82  CD  PRO A  11      27.466  51.952 139.620  1.00 26.05           C  \n\
ATOM     83  N   HIS A  12      24.600  49.477 139.317  1.00 31.56           N  \n\
ATOM     84  CA  HIS A  12      23.749  48.293 139.268  1.00 32.11           C  \n\
ATOM     85  C   HIS A  12      22.709  48.112 140.372  1.00 32.10           C  \n\
ATOM     86  O   HIS A  12      22.015  47.099 140.394  1.00 33.73           O  \n\
ATOM     87  CB  HIS A  12      24.614  47.032 139.195  1.00 33.31           C  \n\
ATOM     88  CG  HIS A  12      25.599  47.035 138.067  1.00 36.62           C  \n\
ATOM     89  ND1 HIS A  12      26.891  47.497 138.207  1.00 38.93           N  \n\
ATOM     90  CD2 HIS A  12      25.484  46.626 136.781  1.00 38.44           C  \n\
ATOM     91  CE1 HIS A  12      27.530  47.368 137.058  1.00 39.08           C  \n\
ATOM     92  NE2 HIS A  12      26.700  46.841 136.176  1.00 40.23           N  \n\
ATOM     93  N   ILE A  13      22.588  49.062 141.292  1.00 32.16           N  \n\
ATOM     94  CA  ILE A  13      21.589  48.910 142.346  1.00 32.50           C  \n\
ATOM     95  C   ILE A  13      20.682  50.123 142.453  1.00 32.60           C  \n\
ATOM     96  O   ILE A  13      21.077  51.228 142.100  1.00 32.90           O  \n\
ATOM     97  CB  ILE A  13      22.235  48.634 143.723  1.00 33.40           C  \n\
ATOM     98  CG1 ILE A  13      23.018  49.852 144.203  1.00 30.64           C  \n\
ATOM     99  CG2 ILE A  13      23.150  47.414 143.627  1.00 35.17           C  \n\
ATOM    100  CD1 ILE A  13      23.484  49.719 145.628  1.00 31.41           C  \n\
ATOM    101  N   ASN A  14      19.462  49.907 142.939  1.00 32.76           N  \n\
ATOM    102  CA  ASN A  14      18.486  50.981 143.067  1.00 32.73           C  \n\
ATOM    103  C   ASN A  14      18.333  51.507 144.490  1.00 33.06           C  \n\
ATOM    104  O   ASN A  14      19.050  51.089 145.402  1.00 34.45           O  \n\
ATOM    105  CB  ASN A  14      17.117  50.522 142.546  1.00 32.10           C  \n\
ATOM    106  CG  ASN A  14      16.552  49.350 143.332  1.00 33.49           C  \n\
ATOM    107  OD1 ASN A  14      16.749  49.246 144.546  1.00 33.12           O  \n\
ATOM    108  ND2 ASN A  14      15.828  48.471 142.644  1.00 33.37           N  \n\
ATOM    109  N   PHE A  15      17.379  52.421 144.660  1.00 31.73           N  \n\
ATOM    110  CA  PHE A  15      17.086  53.061 145.941  1.00 30.58           C  \n\
ATOM    111  C   PHE A  15      16.926  52.103 147.118  1.00 30.34           C  \n\
ATOM    112  O   PHE A  15      17.644  52.206 148.108  1.00 30.58           O  \n\
ATOM    113  CB  PHE A  15      15.813  53.907 145.816  1.00 30.27           C  \n\
ATOM    114  CG  PHE A  15      15.593  54.843 146.966  1.00 29.36           C  \n\
ATOM    115  CD1 PHE A  15      16.246  56.069 147.015  1.00 30.40           C  \n\
ATOM    116  CD2 PHE A  15      14.750  54.494 148.011  1.00 30.97           C  \n\
ATOM    117  CE1 PHE A  15      16.068  56.940 148.089  1.00 30.02           C  \n\
ATOM    118  CE2 PHE A  15      14.560  55.357 149.097  1.00 32.66           C  \n\
ATOM    119  CZ  PHE A  15      15.224  56.586 149.132  1.00 31.71           C  \n\
ATOM    120  N   VAL A  16      15.970  51.183 147.018  1.00 30.12           N  \n\
ATOM    121  CA  VAL A  16      15.705  50.223 148.086  1.00 29.21           C  \n\
ATOM    122  C   VAL A  16      16.918  49.347 148.392  1.00 28.79           C  \n\
ATOM    123  O   VAL A  16      17.214  49.058 149.551  1.00 28.79           O  \n\
ATOM    124  CB  VAL A  16      14.497  49.321 147.730  1.00 29.48           C  \n\
ATOM    125  CG1 VAL A  16      14.303  48.250 148.790  1.00 29.33           C  \n\
ATOM    126  CG2 VAL A  16      13.241  50.171 147.613  1.00 31.28           C  \n\
ATOM    127  N   GLN A  17      17.620  48.925 147.351  1.00 28.97           N  \n\
ATOM    128  CA  GLN A  17      18.798  48.088 147.527  1.00 29.74           C  \n\
ATOM    129  C   GLN A  17      19.902  48.871 148.239  1.00 28.62           C  \n\
ATOM    130  O   GLN A  17      20.581  48.341 149.119  1.00 28.00           O  \n\
ATOM    131  CB  GLN A  17      19.270  47.583 146.165  1.00 32.50           C  \n\
ATOM    132  CG  GLN A  17      18.161  46.820 145.440  1.00 41.70           C  \n\
ATOM    133  CD  GLN A  17      18.509  46.433 144.017  1.00 45.11           C  \n\
ATOM    134  OE1 GLN A  17      18.728  47.289 143.153  1.00 43.79           O  \n\
ATOM    135  NE2 GLN A  17      18.555  45.131 143.764  1.00 48.60           N  \n\
TER     136      GLN A  17                                                      \n\
END                                                                             \n";

	ofstream pdb_fs;
	pdb_fs.open(inputPdb.c_str());
	if (pdb_fs.fail()) {
		cerr << "Error writing test pdb " << inputPdb << endl;
		exit(1);
	}
	pdb_fs << pdbContent;
	pdb_fs.close();

	System sys;
	if (!sys.readPdb(inputPdb)) {
		cerr << "Cannot read pdb file " << inputPdb << endl;
		exit(1);
	}

	AtomPointerVector av = sys.getAtomPointers();

	AtomBondBuilder abb;
	abb.buildConnections(av);
	
	cout << endl;
	cout << "===========================================" << endl;
	cout << endl;
	Atom * pAtomN4 = NULL;
	if (sys.atomExists("A, 4, N")) {
		pAtomN4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomN4 << endl;
	} else {
		cout << "Atom A 4 N not found" << endl;
		exit(1);
	}
	Atom * pAtomCA4 = NULL;
	if (sys.atomExists("A, 4, CA")) {
		pAtomCA4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCA4 << endl;
	} else {
		cout << "Atom A 4 CA not found" << endl;
		exit(1);
	}
	Atom * pAtomC4 = NULL;
	if (sys.atomExists("A, 4, C")) {
		pAtomC4 = &(sys.getLastFoundAtom());
	//	cout << *pAtomC4 << endl;
	} else {
		cout << "Atom A 4 C not found" << endl;
		exit(1);
	}

	// PHE 5
	Atom * pAtomN5 = NULL;
	if (sys.atomExists("A, 5, N")) {
		pAtomN5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomN5 << endl;
	} else {
		cout << "Atom A 5 N not found" << endl;
		exit(1);
	}
	Atom * pAtomCA5 = NULL;
	if (sys.atomExists("A, 5, CA")) {
		pAtomCA5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCA5 << endl;
	} else {
		cout << "Atom A 5 CA not found" << endl;
		exit(1);
	}
	Atom * pAtomCB5 = NULL;
	if (sys.atomExists("A, 5, CB")) {
		pAtomCB5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCB5 << endl;
	} else {
		cout << "Atom A 5 CB not found" << endl;
		exit(1);
	}
	Atom * pAtomCG5 = NULL;
	if (sys.atomExists("A, 5, CG")) {
		pAtomCG5 = &(sys.getLastFoundAtom());
	//	cout << *pAtomCG5 << endl;
	} else {
		cout << "Atom A 5 CG not found" << endl;
		exit(1);
	}

	// save the original coordinates
	cout << "Save the coordinates of the original pdb " << inputPdb << endl;
	av.saveCoor("orig");

	Transforms T;

	cout << " ==============================" << endl;
	cout << "Set the PHE 5 (CA, CB) bond to 1.65" << endl;
	cout << "  Bond before " << pAtomCA5->distance(*pAtomCB5) << endl;
	T.setBondDistance(*pAtomCA5, *pAtomCB5, 1.65);
	cout << "  Bond after " << pAtomCA5->distance(*pAtomCB5) << endl;
	string outPdb = "/tmp/testEdit-bond.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 (N, CA, CB) angle to 112.47" << endl;
	cout << "  Angle before " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	T.setBondAngle(*pAtomN5, *pAtomCA5, *pAtomCB5, 112.47);
	cout << "  Angle after " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	outPdb = "/tmp/testEdit-angle.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 chi1 (N, CA, CB, CG) to 178.34" << endl;
	cout << "  Chi1 before " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	T.setDihedral(*pAtomN5, *pAtomCA5, *pAtomCB5, *pAtomCG5, 178.34);
	cout << "  Chi1 after " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	outPdb = "/tmp/testEdit-chi1.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the THR 4 psi (N, CA, C, N+1) to 175.67" << endl;
	cout << "  Psi before " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	T.setDihedral(*pAtomN4, *pAtomCA4, *pAtomC4, *pAtomN5, 175.67);
	cout << "  Psi after " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	outPdb = "/tmp/testEdit-psi.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	T.setNaturalMovements(true);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 (CA, CB) bond to 1.65 with a natural movement" << endl;
	cout << "  Bond before " << pAtomCA5->distance(*pAtomCB5) << endl;
	T.setBondDistance(*pAtomCA5, *pAtomCB5, 1.65);
	cout << "  Bond after " << pAtomCA5->distance(*pAtomCB5) << endl;
	outPdb = "/tmp/testEdit-bond-nat.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 (N, CA, CB) angle to 112.47 with a natural movement" << endl;
	cout << "  Angle before " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	T.setBondAngle(*pAtomN5, *pAtomCA5, *pAtomCB5, 112.47);
	cout << "  Angle after " << pAtomN5->angle(*pAtomCA5, *pAtomCB5) << endl;
	outPdb = "/tmp/testEdit-angle-nat.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the PHE 5 chi1 (N, CA, CB, CG) to 178.34 with a natural movement" << endl;
	cout << "  Chi1 before " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	T.setDihedral(*pAtomN5, *pAtomCA5, *pAtomCB5, *pAtomCG5, 178.34);
	cout << "  Chi1 after " << pAtomN5->dihedral(*pAtomCA5, *pAtomCB5, *pAtomCG5) << endl;
	outPdb = "/tmp/testEdit-chi1-nat.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

	cout << " ==============================" << endl;
	av.applySavedCoor("orig");
	cout << "Set the THR 4 psi (N, CA, C, N+1) to 175.67 with a natural movement" << endl;
	cout << "  Psi before " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	T.setDihedral(*pAtomN4, *pAtomCA4, *pAtomC4, *pAtomN5, 175.67);
	cout << "  Psi after " << pAtomN4->dihedral(*pAtomCA4, *pAtomC4, *pAtomN5) << endl;
	outPdb = "/tmp/testEdit-psi-nat.pdb";
	cout << "write " << outPdb << endl;
	sys.writePdb(outPdb);

}

