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

#include "MslTools.h"
#include "System.h"
#include "System.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"
#include "AtomSelection.h"
#include "Transforms.h"


using namespace MSL;
using namespace std;

#include "SysEnv.h"
static SysEnv SYSENV;

int main() {

	string inputPdb = "/tmp/dummy3.pdb";

	string pdbContent = "ATOM      1 DUM1 DM3 A   1       0.000   0.000   0.000  1.00  0.00              \n\
ATOM      2 DUM2 DM3 A   1       1.000   0.000   0.000  1.00  0.00              \n\
ATOM      3 DUM3 DM3 A   1       2.000   0.000   0.000  1.00  0.00              \n\
ATOM      4 DUM4 DM3 A   1       1.000   0.000   0.000  1.00  0.00              \n\
TER       5                                                                     \n\
END                                                                             \n";


	ofstream pdb_fs;
	pdb_fs.open(inputPdb.c_str());
	if (pdb_fs.fail()) {
		cerr << "Error writing test pdb " << inputPdb << endl;
		exit(1);
	}
	pdb_fs << pdbContent;
	pdb_fs.close();

	PDBReader reader;
	if (!reader.open(inputPdb) || !reader.read()) {
		cerr << "Cannot read PDB file " << inputPdb << endl;
		exit(1);
	}
	PolymerSequence seq(reader.getAtomPointers());

	cout << seq << endl;

	System sys;

	/************************************************************************
	 *  - call the Charmm system builder with the topology and parameter files
	 *  - read the solvation file
	 *  - set NOT to build the non-bonded interaction at first
	 *  - build the system
	 *  - 
	 ************************************************************************/
	string topfile = SYSENV.getEnv("MSL_DIR")+"/toppar/charmm/charmm/toph19_dum.inp";
	string parfile =  SYSENV.getEnv("MSL_DIR")+"/toppar/charmm/param19_eef1.1.inp";
	string solvfile = SYSENV.getEnv("MSL_DIR")+"/toppar/charmm/solvpar.inp";
	CharmmSystemBuilder CSB(sys, topfile, parfile);
	CSB.readSolvation(solvfile);
	//CharmmSystemBuilder CSB(sys, "/exports/home/asenes/fromThemis/dummyTest/toph19_eef1.1.dum.inp", "/exports/home/asenes/fromThemis/dummyTest/param19_eef1.1.nowildcards.inp");
	//CSB.readSolvation("/exports/home/asenes/fromThemis/dummyTest/solvpar_dum.inp");
	CSB.setBuildNonBondedInteractions(false);
	CSB.setElec14factor(0.4);
	CSB.setUseRdielectric(true);
	CSB.buildSystem(seq);
	sys.assignCoordinates(reader.getAtomPointers());
	CSB.updateNonBonded(997.0, 998.0, 999.0);
	sys.getEnergySet()->setAllTermsInactive();
	sys.getEnergySet()->setTermActive("CHARMM_EEF1");
	sys.getEnergySet()->setTermActive("CHARMM_EF1R");
	cout << sys.getAtomPointers() << endl;

	AtomSelection sel(sys.getAtomPointers());
	AtomPointerVector apv = sel.select("atoms14, name DUM1 or name DUM3");
	cout << apv.size() << endl;

	for (double x=1.0; x<=3.0; x+=0.1) {
		sys[3].setCoor(x, 0.0, 0.0);
		cout << sys[3] << endl;
		cout << sys.calcEnergy() << endl;
		cout << sys.getEnergySummary();
	}
	/*
	map<string, vector<Interaction*> > * pTerms = sys.getEnergySet()->getEnergyTerms();
	
	for (map<string, vector<Interaction*> >::iterator t=pTerms->begin(); t!=pTerms->end(); t++) {
		cout << t->first << endl;
		//if (t->first != "CHARMM_VDW" && t->first != "CHARMM_ANGL") {
		if (t->first != "CHARMM_EEF1" && t->first != "CHARMM_EF1R") {
			continue;
		}
		for (vector<Interaction*>::iterator i=t->second.begin(); i!=t->second.end(); i++) {
			cout << "  " << (*i)->toString() << endl;
		}
	}
	*/

/*
	string inputPdb = "/tmp/dyna.pdb";

	string pdbContent = "ATOM      1  HT1 TYR M   1      -2.987   1.517  14.829  1.00  0.00              \n\
ATOM      2  HT2 TYR M   1      -1.309   1.247  14.576  1.00  0.00              \n\
ATOM      3  HT3 TYR M   1      -2.347  -0.070  14.735  1.00  0.00              \n\
ATOM      4  N   TYR M   1      -2.275   0.905  14.384  1.00  0.00              \n\
ATOM      5  CA  TYR M   1      -2.400   0.904  12.939  1.00  0.00              \n\
ATOM      6  CB  TYR M   1      -2.497   2.351  12.433  1.00  0.00              \n\
ATOM      7  CG  TYR M   1      -2.654   2.526  10.923  1.00  0.00              \n\
ATOM      8  CD1 TYR M   1      -1.683   3.212  10.218  1.00  0.00              \n\
ATOM      9  CE1 TYR M   1      -1.803   3.377   8.853  1.00  0.00              \n\
ATOM     10  CD2 TYR M   1      -3.754   2.011  10.260  1.00  0.00              \n\
ATOM     11  CE2 TYR M   1      -3.877   2.173   8.894  1.00  0.00              \n\
ATOM     12  CZ  TYR M   1      -2.898   2.852   8.199  1.00  0.00              \n\
ATOM     13  OH  TYR M   1      -3.008   3.003   6.829  1.00  0.00              \n\
ATOM     14  HH  TYR M   1      -2.289   3.561   6.525  1.00  0.00              \n\
ATOM     15  C   TYR M   1      -1.168   0.203  12.378  1.00  0.00              \n\
ATOM     16  O   TYR M   1      -0.092   0.264  12.977  1.00  0.00              \n\
ATOM     17  N   GLY M   2      -1.340  -0.484  11.251  1.00  0.00              \n\
ATOM     18  H   GLY M   2      -2.221  -0.484  10.824  1.00  0.00              \n\
ATOM     19  CA  GLY M   2      -0.241  -1.140  10.565  1.00  0.00              \n\
ATOM     20  C   GLY M   2       0.040  -0.428   9.252  1.00  0.00              \n\
ATOM     21  O   GLY M   2      -0.853   0.203   8.679  1.00  0.00              \n\
ATOM     22  N   GLY M   3       1.268  -0.528   8.739  1.00  0.00              \n\
ATOM     23  H   GLY M   3       1.919  -1.098   9.193  1.00  0.00              \n\
ATOM     24  CA  GLY M   3       1.670   0.170   7.524  1.00  0.00              \n\
ATOM     25  C   GLY M   3       1.098  -0.453   6.257  1.00  0.00              \n\
ATOM     26  O   GLY M   3       1.814  -1.070   5.465  1.00  0.00              \n\
ATOM     27  N   PHE M   4      -0.210  -0.297   6.048  1.00  0.00              \n\
ATOM     28  H   PHE M   4      -0.739   0.187   6.721  1.00  0.00              \n\
ATOM     29  CA  PHE M   4      -0.870  -0.849   4.876  1.00  0.00              \n\
ATOM     30  CB  PHE M   4      -2.383  -0.955   5.093  1.00  0.00              \n\
ATOM     31  CG  PHE M   4      -2.798  -1.728   6.339  1.00  0.00              \n\
ATOM     32  CD1 PHE M   4      -3.619  -1.116   7.267  1.00  0.00              \n\
ATOM     33  CD2 PHE M   4      -2.351  -3.019   6.555  1.00  0.00              \n\
ATOM     34  CE1 PHE M   4      -3.980  -1.791   8.414  1.00  0.00              \n\
ATOM     35  CE2 PHE M   4      -2.718  -3.686   7.706  1.00  0.00              \n\
ATOM     36  CZ  PHE M   4      -3.530  -3.074   8.636  1.00  0.00              \n\
ATOM     37  C   PHE M   4      -0.555  -0.005   3.647  1.00  0.00              \n\
ATOM     38  O   PHE M   4      -1.182   1.020   3.368  1.00  0.00              \n\
ATOM     39  N   LEU M   5       0.507  -0.401   2.947  1.00  0.00              \n\
ATOM     40  H   LEU M   5       1.104  -1.060   3.358  1.00  0.00              \n\
ATOM     41  CA  LEU M   5       0.846   0.200   1.667  1.00  0.00              \n\
ATOM     42  CB  LEU M   5       2.322  -0.040   1.342  1.00  0.00              \n\
ATOM     43  CG  LEU M   5       3.388   0.518   2.280  1.00  0.00              \n\
ATOM     44  CD1 LEU M   5       4.758   0.016   1.862  1.00  0.00              \n\
ATOM     45  CD2 LEU M   5       3.358   2.039   2.312  1.00  0.00              \n\
ATOM     46  C   LEU M   5      -0.044  -0.351   0.558  1.00  0.00              \n\
ATOM     47  O   LEU M   5       0.221  -1.418  -0.005  1.00  0.00              \n\
ATOM     48  N   ARG M   6      -1.131   0.369   0.267  1.00  0.00              \n\
ATOM     49  H   ARG M   6      -1.316   1.165   0.812  1.00  0.00              \n\
ATOM     50  CA  ARG M   6      -2.089  -0.008  -0.770  1.00  0.00              \n\
ATOM     51  CB  ARG M   6      -3.199   1.036  -0.918  1.00  0.00              \n\
ATOM     52  CG  ARG M   6      -4.304   1.035   0.139  1.00  0.00              \n\
ATOM     53  CD  ARG M   6      -3.847   1.561   1.491  1.00  0.00              \n\
ATOM     54  NE  ARG M   6      -4.942   1.594   2.444  1.00  0.00              \n\
ATOM     55  HE  ARG M   6      -5.844   1.380   2.124  1.00  0.00              \n\
ATOM     56  CZ  ARG M   6      -4.753   1.898   3.732  1.00  0.00              \n\
ATOM     57  NH1 ARG M   6      -5.805   1.978   4.543  1.00  0.00              \n\
ATOM     58 HH11 ARG M   6      -5.682   2.208   5.509  1.00  0.00              \n\
ATOM     59 HH12 ARG M   6      -6.725   1.813   4.186  1.00  0.00              \n\
ATOM     60  NH2 ARG M   6      -3.531   2.090   4.235  1.00  0.00              \n\
ATOM     61 HH21 ARG M   6      -3.403   2.311   5.201  1.00  0.00              \n\
ATOM     62 HH22 ARG M   6      -2.727   1.961   3.653  1.00  0.00              \n\
ATOM     63  C   ARG M   6      -1.483  -0.273  -2.143  1.00  0.00              \n\
ATOM     64  O   ARG M   6      -0.588   0.438  -2.611  1.00  0.00              \n\
ATOM     65  N   ARG M   7      -1.987  -1.318  -2.800  1.00  0.00              \n\
ATOM     66  H   ARG M   7      -2.740  -1.806  -2.399  1.00  0.00              \n\
ATOM     67  CA  ARG M   7      -1.458  -1.791  -4.071  1.00  0.00              \n\
ATOM     68  CB  ARG M   7      -1.821  -3.270  -4.264  1.00  0.00              \n\
ATOM     69  CG  ARG M   7      -0.881  -4.295  -3.632  1.00  0.00              \n\
ATOM     70  CD  ARG M   7      -0.824  -4.257  -2.111  1.00  0.00              \n\
ATOM     71  NE  ARG M   7       0.099  -5.268  -1.623  1.00  0.00              \n\
ATOM     72  HE  ARG M   7       0.200  -6.092  -2.145  1.00  0.00              \n\
ATOM     73  CZ  ARG M   7       0.794  -5.120  -0.493  1.00  0.00              \n\
ATOM     74  NH1 ARG M   7       1.581  -6.115  -0.090  1.00  0.00              \n\
ATOM     75 HH11 ARG M   7       2.113  -6.040   0.753  1.00  0.00              \n\
ATOM     76 HH12 ARG M   7       1.635  -6.955  -0.631  1.00  0.00              \n\
ATOM     77  NH2 ARG M   7       0.729  -3.996   0.225  1.00  0.00              \n\
ATOM     78 HH21 ARG M   7       1.254  -3.889   1.069  1.00  0.00              \n\
ATOM     79 HH22 ARG M   7       0.164  -3.230  -0.085  1.00  0.00              \n\
ATOM     80  C   ARG M   7      -1.826  -0.964  -5.305  1.00  0.00              \n\
ATOM     81  O   ARG M   7      -2.443  -1.436  -6.262  1.00  0.00              \n\
ATOM     82  N   ILE M   8      -1.438   0.312  -5.298  1.00  0.00              \n\
ATOM     83  H   ILE M   8      -1.077   0.685  -4.460  1.00  0.00              \n\
ATOM     84  CA  ILE M   8      -1.541   1.156  -6.483  1.00  0.00              \n\
ATOM     85  CB  ILE M   8      -1.704   2.649  -6.062  1.00  0.00              \n\
ATOM     86  CG2 ILE M   8      -0.505   3.173  -5.275  1.00  0.00              \n\
ATOM     87  CG1 ILE M   8      -2.013   3.572  -7.238  1.00  0.00              \n\
ATOM     88  CD  ILE M   8      -3.373   3.312  -7.919  1.00  0.00              \n\
ATOM     89  C   ILE M   8      -0.344   0.915  -7.407  1.00  0.00              \n\
ATOM     90  O   ILE M   8      -0.438   1.002  -8.635  1.00  0.00              \n\
ATOM     91  N   ARG M   9       0.773   0.513  -6.786  1.00  0.00              \n\
ATOM     92  H   ARG M   9       0.704   0.325  -5.822  1.00  0.00              \n\
ATOM     93  CA  ARG M   9       2.067   0.323  -7.432  1.00  0.00              \n\
ATOM     94  CB  ARG M   9       3.088  -0.215  -6.418  1.00  0.00              \n\
ATOM     95  CG  ARG M   9       3.121   0.533  -5.085  1.00  0.00              \n\
ATOM     96  CD  ARG M   9       3.958  -0.199  -4.038  1.00  0.00              \n\
ATOM     97  NE  ARG M   9       3.519  -1.569  -3.803  1.00  0.00              \n\
ATOM     98  HE  ARG M   9       3.957  -2.284  -4.312  1.00  0.00              \n\
ATOM     99  CZ  ARG M   9       2.551  -1.896  -2.940  1.00  0.00              \n\
ATOM    100  NH1 ARG M   9       2.257  -3.182  -2.768  1.00  0.00              \n\
ATOM    101 HH11 ARG M   9       1.526  -3.466  -2.148  1.00  0.00              \n\
ATOM    102 HH12 ARG M   9       2.751  -3.880  -3.286  1.00  0.00              \n\
ATOM    103  NH2 ARG M   9       1.862  -0.970  -2.269  1.00  0.00              \n\
ATOM    104 HH21 ARG M   9       1.139  -1.224  -1.627  1.00  0.00              \n\
ATOM    105 HH22 ARG M   9       2.049  -0.000  -2.433  1.00  0.00              \n\
ATOM    106  C   ARG M   9       2.036  -0.551  -8.691  1.00  0.00              \n\
ATOM    107  O   ARG M   9       2.606  -0.113  -9.694  1.00  0.00              \n\
ATOM    108  N   PRO M  10       1.401  -1.742  -8.775  1.00  0.00              \n\
ATOM    109  CD  PRO M  10       0.901  -2.548  -7.661  1.00  0.00              \n\
ATOM    110  CA  PRO M  10       1.245  -2.478 -10.027  1.00  0.00              \n\
ATOM    111  CB  PRO M  10       0.594  -3.781  -9.594  1.00  0.00              \n\
ATOM    112  CG  PRO M  10      -0.137  -3.429  -8.326  1.00  0.00              \n\
ATOM    113  C   PRO M  10       0.425  -1.737 -11.079  1.00  0.00              \n\
ATOM    114  O   PRO M  10       0.795  -1.744 -12.254  1.00  0.00              \n\
ATOM    115  N   LYS M  11      -0.647  -1.033 -10.688  1.00  0.00              \n\
ATOM    116  H   LYS M  11      -0.805  -0.898  -9.728  1.00  0.00              \n\
ATOM    117  CA  LYS M  11      -1.524  -0.358 -11.637  1.00  0.00              \n\
ATOM    118  CB  LYS M  11      -2.808   0.106 -10.941  1.00  0.00              \n\
ATOM    119  CG  LYS M  11      -3.909   0.612 -11.875  1.00  0.00              \n\
ATOM    120  CD  LYS M  11      -4.457  -0.510 -12.755  1.00  0.00              \n\
ATOM    121  CE  LYS M  11      -5.485  -0.001 -13.761  1.00  0.00              \n\
ATOM    122  NZ  LYS M  11      -4.873   0.814 -14.796  1.00  0.00              \n\
ATOM    123  HZ1 LYS M  11      -4.392   1.628 -14.362  1.00  0.00              \n\
ATOM    124  HZ2 LYS M  11      -4.178   0.245 -15.321  1.00  0.00              \n\
ATOM    125  HZ3 LYS M  11      -5.605   1.154 -15.451  1.00  0.00              \n\
ATOM    126  C   LYS M  11      -0.809   0.816 -12.300  1.00  0.00              \n\
ATOM    127  O   LYS M  11      -1.019   1.094 -13.484  1.00  0.00              \n\
ATOM    128  N   LEU M  12       0.076   1.482 -11.547  1.00  0.00              \n\
ATOM    129  H   LEU M  12       0.112   1.277 -10.585  1.00  0.00              \n\
ATOM    130  CA  LEU M  12       0.950   2.514 -12.093  1.00  0.00              \n\
ATOM    131  CB  LEU M  12       1.832   3.107 -10.995  1.00  0.00              \n\
ATOM    132  CG  LEU M  12       1.170   3.878  -9.857  1.00  0.00              \n\
ATOM    133  CD1 LEU M  12       2.194   4.218  -8.788  1.00  0.00              \n\
ATOM    134  CD2 LEU M  12       0.482   5.137 -10.364  1.00  0.00              \n\
ATOM    135  C   LEU M  12       1.821   2.019 -13.245  1.00  0.00              \n\
ATOM    136  O   LEU M  12       2.052   2.754 -14.207  1.00  0.00              \n\
ATOM    137  N   LYS M  13       2.294   0.768 -13.180  1.00  0.00              \n\
ATOM    138  H   LYS M  13       2.150   0.245 -12.360  1.00  0.00              \n\
ATOM    139  CA  LYS M  13       2.998   0.158 -14.300  1.00  0.00              \n\
ATOM    140  CB  LYS M  13       3.837  -1.030 -13.821  1.00  0.00              \n\
ATOM    141  CG  LYS M  13       4.693  -1.655 -14.919  1.00  0.00              \n\
ATOM    142  CD  LYS M  13       5.561  -2.779 -14.374  1.00  0.00              \n\
ATOM    143  CE  LYS M  13       6.340  -3.473 -15.490  1.00  0.00              \n\
ATOM    144  NZ  LYS M  13       7.270  -2.575 -16.151  1.00  0.00              \n\
ATOM    145  HZ1 LYS M  13       6.745  -1.789 -16.584  1.00  0.00              \n\
ATOM    146  HZ2 LYS M  13       7.946  -2.201 -15.455  1.00  0.00              \n\
ATOM    147  HZ3 LYS M  13       7.785  -3.093 -16.891  1.00  0.00              \n\
ATOM    148  C   LYS M  13       2.013  -0.262 -15.390  1.00  0.00              \n\
ATOM    149  O   LYS M  13       2.191   0.093 -16.557  1.00  0.00              \n\
ATOM    150  N   TRP M  14       0.936  -0.972 -15.024  1.00  0.00              \n\
ATOM    151  H   TRP M  14       0.828  -1.182 -14.071  1.00  0.00              \n\
ATOM    152  CA  TRP M  14      -0.064  -1.466 -15.970  1.00  0.00              \n\
ATOM    153  CB  TRP M  14      -1.144  -2.284 -15.256  1.00  0.00              \n\
ATOM    154  CG  TRP M  14      -0.685  -3.605 -14.646  1.00  0.00              \n\
ATOM    155  CD2 TRP M  14      -1.173  -4.185 -13.507  1.00  0.00              \n\
ATOM    156  CE2 TRP M  14      -0.454  -5.334 -13.430  1.00  0.00              \n\
ATOM    157  CE3 TRP M  14      -2.292  -3.966 -12.740  1.00  0.00              \n\
ATOM    158  CD1 TRP M  14       0.298  -4.394 -15.214  1.00  0.00              \n\
ATOM    159  NE1 TRP M  14       0.410  -5.458 -14.443  1.00  0.00              \n\
ATOM    160  HE1 TRP M  14       1.024  -6.207 -14.595  1.00  0.00              \n\
ATOM    161  CZ2 TRP M  14      -0.770  -6.367 -12.584  1.00  0.00              \n\
ATOM    162  CZ3 TRP M  14      -2.624  -4.982 -11.855  1.00  0.00              \n\
ATOM    163  CH2 TRP M  14      -1.880  -6.157 -11.778  1.00  0.00              \n\
ATOM    164  C   TRP M  14      -0.746  -0.403 -16.825  1.00  0.00              \n\
ATOM    165  O   TRP M  14      -1.234  -0.702 -17.917  1.00  0.00              \n\
ATOM    166  N   ASP M  15      -0.801   0.844 -16.342  1.00  0.00              \n\
ATOM    167  H   ASP M  15      -0.552   0.985 -15.400  1.00  0.00              \n\
ATOM    168  CA  ASP M  15      -1.281   1.971 -17.133  1.00  0.00              \n\
ATOM    169  CB  ASP M  15      -1.367   3.225 -16.259  1.00  0.00              \n\
ATOM    170  CG  ASP M  15      -1.959   4.443 -16.962  1.00  0.00              \n\
ATOM    171  OD1 ASP M  15      -3.174   4.483 -17.154  1.00  0.00              \n\
ATOM    172  OD2 ASP M  15      -1.206   5.357 -17.294  1.00  0.00              \n\
ATOM    173  C   ASP M  15      -0.421   2.256 -18.363  1.00  0.00              \n\
ATOM    174  O   ASP M  15      -0.946   2.653 -19.405  1.00  0.00              \n\
ATOM    175  N   ASN M  16       0.897   2.034 -18.292  1.00  0.00              \n\
ATOM    176  H   ASN M  16       1.270   1.584 -17.502  1.00  0.00              \n\
ATOM    177  CA  ASN M  16       1.793   2.346 -19.399  1.00  0.00              \n\
ATOM    178  CB  ASN M  16       3.194   2.649 -18.850  1.00  0.00              \n\
ATOM    179  CG  ASN M  16       4.171   3.215 -19.876  1.00  0.00              \n\
ATOM    180  OD1 ASN M  16       4.402   4.419 -19.940  1.00  0.00              \n\
ATOM    181  ND2 ASN M  16       4.786   2.394 -20.717  1.00  0.00              \n\
ATOM    182 HD21 ASN M  16       4.560   1.439 -20.672  1.00  0.00              \n\
ATOM    183 HD22 ASN M  16       5.429   2.783 -21.342  1.00  0.00              \n\
ATOM    184  C   ASN M  16       1.804   1.208 -20.423  1.00  0.00              \n\
ATOM    185  O   ASN M  16       2.775   0.462 -20.578  1.00  0.00              \n\
ATOM    186  N   GLN M  17       0.684   1.063 -21.134  1.00  0.00              \n\
ATOM    187  H   GLN M  17      -0.071   1.664 -20.936  1.00  0.00              \n\
ATOM    188  CA  GLN M  17       0.559   0.061 -22.186  1.00  0.00              \n\
ATOM    189  CB  GLN M  17      -0.907  -0.099 -22.596  1.00  0.00              \n\
ATOM    190  CG  GLN M  17      -1.800  -0.772 -21.553  1.00  0.00              \n\
ATOM    191  CD  GLN M  17      -1.478  -2.240 -21.303  1.00  0.00              \n\
ATOM    192  OE1 GLN M  17      -1.452  -3.063 -22.215  1.00  0.00              \n\
ATOM    193  NE2 GLN M  17      -1.227  -2.633 -20.063  1.00  0.00              \n\
ATOM    194 HE21 GLN M  17      -1.261  -1.966 -19.348  1.00  0.00              \n\
ATOM    195 HE22 GLN M  17      -1.019  -3.579 -19.923  1.00  0.00              \n\
ATOM    196  C   GLN M  17       1.423   0.347 -23.413  1.00  0.00              \n\
ATOM    197  OT1 GLN M  17       2.042  -0.587 -23.921  1.00  0.00              \n\
ATOM    198  OT2 GLN M  17       1.471   1.497 -23.853  1.00  0.00              \n\
TER     199                                                                     \n\
END                                                                             \n";
	ofstream pdb_fs;
	pdb_fs.open(inputPdb.c_str());
	if (pdb_fs.fail()) {
		cerr << "Error writing test pdb " << inputPdb << endl;
		exit(1);
	}
	pdb_fs << pdbContent;
	pdb_fs.close();

	PDBReader reader;
	if (!reader.open(inputPdb) || !reader.read()) {
		cerr << "Cannot read PDB file " << inputPdb << endl;
		exit(1);
	}
	PolymerSequence seq(reader.getAtomPointers());

	cout << seq << endl;

	System sys;

	/ ************************************************************************
	 *  - call the Charmm system builder with the topology and parameter files
	 *  - read the solvation file
	 *  - set NOT to build the non-bonded interaction at first
	 *  - build the system
	 *  - 
	 ************************************************************************ /
	CharmmSystemBuilder CSB("/library/charmmTopPar/toph19_eef1.1.inp", "/library/charmmTopPar/param19_eef1.1.nowildcards.inp");
	CSB.readSolvation("/library/charmmTopPar/osolvpar.inp");
	CSB.setBuildNonBondedInteractions(false);
	CSB.setElec14factor(0.4);
	CSB.setUseRdielectric(true);
	CSB.buildSystem(sys, seq);
	sys.assignCoordinates(reader.getAtomPointers());
	CSB.updateNonBonded(sys, 997.0, 998.0, 999.0);

	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();
	map<string, vector<Interaction*> > * pTerms = sys.getEnergySet()->getEnergyTerms();
	
	for (map<string, vector<Interaction*> >::iterator t=pTerms->begin(); t!=pTerms->end(); t++) {
		cout << t->first << endl;
		//if (t->first != "CHARMM_VDW" && t->first != "CHARMM_ANGL") {
		if (t->first != "CHARMM_EEF1" && t->first != "CHARMM_EF1R") {
			continue;
		}
		for (vector<Interaction*>::iterator i=t->second.begin(); i!=t->second.end(); i++) {
			cout << "  " << (*i)->toString() << endl;
		}
	}
*/

/*

	PolymerSequence seq("\
A: ALA-ACE ILE VAL ILE\n\
B: ARG HSD THR GLY");

	cout << seq << endl;

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
	AtomSelection sel(sys.getAllAtomPointers());
	sel.select("chainB, chain B");

	Transforms tr;
	tr.translate(sel.getSelection("chainB"), CartesianPoint(10,5,5));
	string filename = "/tmp/buildFromCharmmTopology.pdb";

	if (!sys.writePdb(filename)) {
		cerr << "Cannot write output file " << filename << endl;
		exit(1);
	}

	cout << "Written pdb file " << filename << endl;
	cout << endl;

	cout << sys.getAtomPointers();

	cout << sys.calcEnergy() << endl;
	cout << sys.getEnergySummary();

	return 0;
}
*/
}


