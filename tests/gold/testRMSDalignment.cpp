/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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
#include <iomanip> 

#include "PDBReader.h"
#include "PDBWriter.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"

using namespace std;

using namespace MSL;


int main() {


	bool result = true;
	double epsilon = 1e-8;



	string pdb1text = "\
ATOM      1  N   PHE A  10     -39.268  23.899  15.234  1.00 62.84           N  \n\
ATOM      2  CA  PHE A  10     -40.538  24.408  14.746  1.00 64.72           C  \n\
ATOM      3  C   PHE A  10     -40.461  25.871  14.317  1.00 65.69           C  \n\
ATOM      4  O   PHE A  10     -41.252  26.327  13.494  1.00 66.60           O  \n\
ATOM      5  CB  PHE A  10     -41.600  24.247  15.830  1.00 64.93           C  \n\
ATOM      6  CG  PHE A  10     -42.064  22.833  16.013  1.00 67.56           C  \n\
ATOM      7  CD1 PHE A  10     -42.362  22.344  17.281  1.00 69.28           C  \n\
ATOM      8  CD2 PHE A  10     -42.235  21.994  14.913  1.00 69.13           C  \n\
ATOM      9  CE1 PHE A  10     -42.828  21.038  17.459  1.00 70.07           C  \n\
ATOM     10  CE2 PHE A  10     -42.700  20.689  15.072  1.00 70.17           C  \n\
ATOM     11  CZ  PHE A  10     -42.998  20.208  16.350  1.00 70.70           C  \n\
ATOM     12  N   GLY A  11     -39.507  26.603  14.879  1.00 63.84           N  \n\
ATOM     13  CA  GLY A  11     -39.359  27.998  14.529  1.00 59.41           C  \n\
ATOM     14  C   GLY A  11     -38.645  28.182  13.206  1.00 58.26           C  \n\
ATOM     15  O   GLY A  11     -39.097  28.969  12.371  1.00 57.74           O  \n\
ATOM     16  N   THR A  12     -37.544  27.458  13.006  1.00 55.52           N  \n\
ATOM     17  CA  THR A  12     -36.757  27.570  11.777  1.00 54.94           C  \n\
ATOM     18  C   THR A  12     -37.494  26.923  10.610  1.00 55.40           C  \n\
ATOM     19  O   THR A  12     -37.239  27.217   9.439  1.00 56.00           O  \n\
ATOM     20  CB  THR A  12     -35.355  26.925  11.943  1.00 54.83           C  \n\
ATOM     21  OG1 THR A  12     -34.753  27.411  13.152  1.00 55.97           O  \n\
ATOM     22  CG2 THR A  12     -34.454  27.291  10.766  1.00 48.04           C  \n\
ATOM     23  N   PHE A  13     -38.414  26.035  10.946  1.00 54.33           N  \n\
ATOM     24  CA  PHE A  13     -39.228  25.391   9.944  1.00 53.91           C  \n\
ATOM     25  C   PHE A  13     -40.208  26.477   9.465  1.00 53.42           C  \n\
ATOM     26  O   PHE A  13     -40.407  26.660   8.266  1.00 54.42           O  \n\
ATOM     27  CB  PHE A  13     -39.959  24.206  10.582  1.00 54.50           C  \n\
ATOM     28  CG  PHE A  13     -41.082  23.662   9.757  1.00 57.21           C  \n\
ATOM     29  CD1 PHE A  13     -42.313  24.318   9.711  1.00 57.37           C  \n\
ATOM     30  CD2 PHE A  13     -40.916  22.501   9.020  1.00 58.99           C  \n\
ATOM     31  CE1 PHE A  13     -43.359  23.825   8.943  1.00 58.00           C  \n\
ATOM     32  CE2 PHE A  13     -41.959  21.996   8.245  1.00 60.95           C  \n\
ATOM     33  CZ  PHE A  13     -43.186  22.662   8.208  1.00 59.66           C  \n\
ATOM     34  N   TRP A  14     -40.797  27.205  10.411  1.00 50.52           N  \n\
ATOM     35  CA  TRP A  14     -41.754  28.266  10.099  1.00 50.02           C  \n\
ATOM     36  C   TRP A  14     -41.099  29.315   9.210  1.00 49.53           C  \n\
ATOM     37  O   TRP A  14     -41.658  29.714   8.199  1.00 50.95           O  \n\
ATOM     38  CB  TRP A  14     -42.239  28.933  11.389  1.00 47.35           C  \n\
ATOM     39  CG  TRP A  14     -43.457  29.816  11.246  1.00 46.53           C  \n\
ATOM     40  CD1 TRP A  14     -43.569  31.121  11.634  1.00 45.34           C  \n\
ATOM     41  CD2 TRP A  14     -44.765  29.422  10.798  1.00 45.81           C  \n\
ATOM     42  NE1 TRP A  14     -44.858  31.555  11.469  1.00 47.12           N  \n\
ATOM     43  CE2 TRP A  14     -45.614  30.533  10.959  1.00 45.58           C  \n\
ATOM     44  CE3 TRP A  14     -45.299  28.233  10.282  1.00 45.60           C  \n\
ATOM     45  CZ2 TRP A  14     -46.974  30.493  10.625  1.00 47.45           C  \n\
ATOM     46  CZ3 TRP A  14     -46.649  28.193   9.952  1.00 47.81           C  \n\
ATOM     47  CH2 TRP A  14     -47.473  29.318  10.127  1.00 47.10           C  \n\
ATOM     48  N   LEU A  15     -39.906  29.746   9.607  1.00 48.13           N  \n\
ATOM     49  CA  LEU A  15     -39.146  30.746   8.876  1.00 46.11           C  \n\
ATOM     50  C   LEU A  15     -38.887  30.310   7.432  1.00 46.12           C  \n\
ATOM     51  O   LEU A  15     -39.109  31.072   6.495  1.00 47.84           O  \n\
ATOM     52  CB  LEU A  15     -37.826  31.028   9.614  1.00 42.25           C  \n\
ATOM     53  CG  LEU A  15     -36.848  32.058   9.045  1.00 42.99           C  \n\
ATOM     54  CD1 LEU A  15     -35.848  32.506  10.117  1.00 42.65           C  \n\
ATOM     55  CD2 LEU A  15     -36.124  31.468   7.864  1.00 43.00           C  \n\
END                                                                             \n";

	string pdb2text = "\
ATOM      1  N   PHE A  10      11.910   5.611  31.387  1.00 41.00           N  \n\
ATOM      2  CA  PHE A  10      10.641   5.669  30.671  1.00 41.26           C  \n\
ATOM      3  C   PHE A  10      10.365   7.118  30.284  1.00 38.94           C  \n\
ATOM      4  O   PHE A  10       9.859   7.397  29.195  1.00 40.52           O  \n\
ATOM      5  CB  PHE A  10       9.522   5.121  31.547  1.00 45.26           C  \n\
ATOM      6  CG  PHE A  10       9.681   3.673  31.866  1.00 52.10           C  \n\
ATOM      7  CD1 PHE A  10       9.418   2.702  30.905  1.00 54.89           C  \n\
ATOM      8  CD2 PHE A  10      10.133   3.271  33.118  1.00 55.64           C  \n\
ATOM      9  CE1 PHE A  10       9.605   1.337  31.192  1.00 57.91           C  \n\
ATOM     10  CE2 PHE A  10      10.326   1.909  33.417  1.00 58.52           C  \n\
ATOM     11  CZ  PHE A  10      10.061   0.942  32.452  1.00 57.90           C  \n\
ATOM     12  N   GLY A  11      10.719   8.040  31.173  1.00 34.92           N  \n\
ATOM     13  CA  GLY A  11      10.521   9.443  30.877  1.00 30.92           C  \n\
ATOM     14  C   GLY A  11      11.150   9.782  29.541  1.00 28.80           C  \n\
ATOM     15  O   GLY A  11      10.461  10.157  28.597  1.00 27.51           O  \n\
ATOM     16  N   THR A  12      12.463   9.617  29.459  1.00 26.69           N  \n\
ATOM     17  CA  THR A  12      13.212   9.915  28.248  1.00 24.33           C  \n\
ATOM     18  C   THR A  12      12.783   9.060  27.060  1.00 24.25           C  \n\
ATOM     19  O   THR A  12      12.690   9.555  25.930  1.00 23.87           O  \n\
ATOM     20  CB  THR A  12      14.709   9.697  28.482  1.00 25.08           C  \n\
ATOM     21  OG1 THR A  12      15.127  10.467  29.614  1.00 26.79           O  \n\
ATOM     22  CG2 THR A  12      15.505  10.107  27.265  1.00 25.00           C  \n\
ATOM     23  N   PHE A  13      12.562   7.768  27.308  1.00 23.67           N  \n\
ATOM     24  CA  PHE A  13      12.131   6.852  26.258  1.00 22.81           C  \n\
ATOM     25  C   PHE A  13      10.882   7.480  25.651  1.00 22.97           C  \n\
ATOM     26  O   PHE A  13      10.746   7.575  24.427  1.00 24.20           O  \n\
ATOM     27  CB  PHE A  13      11.823   5.468  26.850  1.00 22.28           C  \n\
ATOM     28  CG  PHE A  13      11.063   4.563  25.916  1.00 23.44           C  \n\
ATOM     29  CD1 PHE A  13       9.689   4.755  25.703  1.00 23.48           C  \n\
ATOM     30  CD2 PHE A  13      11.721   3.540  25.217  1.00 22.58           C  \n\
ATOM     31  CE1 PHE A  13       8.974   3.951  24.811  1.00 24.23           C  \n\
ATOM     32  CE2 PHE A  13      11.018   2.726  24.321  1.00 24.12           C  \n\
ATOM     33  CZ  PHE A  13       9.634   2.932  24.114  1.00 24.83           C  \n\
ATOM     34  N   TRP A  14       9.987   7.931  26.525  1.00 21.71           N  \n\
ATOM     35  CA  TRP A  14       8.763   8.584  26.106  1.00 19.63           C  \n\
ATOM     36  C   TRP A  14       9.085   9.858  25.338  1.00 18.46           C  \n\
ATOM     37  O   TRP A  14       8.508  10.129  24.288  1.00 19.33           O  \n\
ATOM     38  CB  TRP A  14       7.927   8.936  27.318  1.00 20.83           C  \n\
ATOM     39  CG  TRP A  14       6.547   9.364  26.983  1.00 23.31           C  \n\
ATOM     40  CD1 TRP A  14       5.997  10.592  27.187  1.00 23.86           C  \n\
ATOM     41  CD2 TRP A  14       5.502   8.534  26.467  1.00 23.98           C  \n\
ATOM     42  NE1 TRP A  14       4.667  10.577  26.844  1.00 25.40           N  \n\
ATOM     43  CE2 TRP A  14       4.339   9.323  26.402  1.00 25.41           C  \n\
ATOM     44  CE3 TRP A  14       5.435   7.196  26.067  1.00 24.31           C  \n\
ATOM     45  CZ2 TRP A  14       3.121   8.817  25.948  1.00 26.20           C  \n\
ATOM     46  CZ3 TRP A  14       4.227   6.691  25.618  1.00 24.08           C  \n\
ATOM     47  CH2 TRP A  14       3.086   7.497  25.565  1.00 26.45           C  \n\
ATOM     48  N   LEU A  15      10.013  10.648  25.859  1.00 17.77           N  \n\
ATOM     49  CA  LEU A  15      10.381  11.890  25.190  1.00 16.70           C  \n\
ATOM     50  C   LEU A  15      10.792  11.596  23.763  1.00 17.15           C  \n\
ATOM     51  O   LEU A  15      10.178  12.096  22.825  1.00 18.62           O  \n\
ATOM     52  CB  LEU A  15      11.528  12.595  25.932  1.00 15.06           C  \n\
ATOM     53  CG  LEU A  15      12.051  13.912  25.347  1.00 11.37           C  \n\
ATOM     54  CD1 LEU A  15      12.949  14.603  26.344  1.00 11.28           C  \n\
ATOM     55  CD2 LEU A  15      12.817  13.651  24.079  1.00 10.42           C  \n\
END                                                                             \n";


	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Write two pdb file and read them into atom vectors   |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	
	// write the input test PDB file
	ofstream pdb_fs;
	pdb_fs.open("/tmp/testPdb1.pdb");
	pdb_fs << pdb1text;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file testPdb1.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file testPdb1.pdb" << endl;
	}
	pdb_fs.close();

	pdb_fs.open("/tmp/testPdb2.pdb");
	pdb_fs << pdb2text;
	if (pdb_fs.fail()) {
		cerr << "Cannot write test input pdb file testPdb2.pdb" << endl;
		exit(1);
	} else {
		cout << "Written test input pdb file testPdb2.pdb" << endl;
	}
	pdb_fs.close();

	PDBReader rAv1;
	rAv1.open("/tmp/testPdb1.pdb");
	rAv1.read();
	AtomPointerVector av1 = rAv1.getAtomPointers();
	cout << "Read atom vector 1 with size " << av1.size() << endl;
	rAv1.close();
	for (AtomPointerVector::iterator k = av1.begin(); k != av1.end() ; k++){
		cout << *(*k) << endl;
	}

	PDBReader rAv2;
	rAv2.open("/tmp/testPdb2.pdb");
	rAv2.read();
	AtomPointerVector av2 = rAv2.getAtomPointers();
	cout << "Read atom vector 2 with size " << av2.size() << endl;
	rAv2.close();
	for (AtomPointerVector::iterator k = av2.begin(); k != av2.end() ; k++){
		cout << *(*k) << endl;
	}

	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Align with regular align using the atoms order       |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;
	Transforms tr;
	double RMSD1 = 0.0;
	if (tr.rmsdAlignment(av2, av1)) {
		RMSD1 = tr.getLastRMSD();
		cout << "Transform sucessful: RMSD = " << setprecision(15) << RMSD1 << endl;
		cout << "Rotation Matrix: " << tr.getLastRotationMatrix() << endl;
		cout << "Translation Vector: " << tr.getLastTranslation() <<endl;
	} else {
		cout << "Transform failed" << endl;
		result = false;
	}

	PDBWriter wAv2;
	wAv2.open("/tmp/testPdb2-aligned.pdb");
	wAv2.write(av2);
	wAv2.close();

//	cout << setprecision(15) << av2[3]->getX() << " " << av2[3]->getY() << " " << av2[3]->getZ() << endl;
//	cout << setprecision(15) << av2[25]->getX() << " " << av2[25]->getY() << " " << av2[25]->getZ() << endl;
//	cout << setprecision(15) << av2[48]->getX() << " " << av2[48]->getY() << " " << av2[48]->getZ() << endl;

	// test the coordinates against the expected result
	CartesianPoint c03(-40.8389357548435, 26.5315646815496, 13.1257933715629);
	CartesianPoint c25(-40.4135286156094, 26.5638110767258, 8.29151308015583);
	CartesianPoint c48(-39.0740041171713, 30.6736513643575, 9.0974323864413);


	cout << endl;
	cout << "=============================================================" << endl;
	cout << "|                                                           |" << endl;
	cout << "|      Align with smart align using shuffled atoms          |" << endl;
	cout << "|                                                           |" << endl;
	cout << "=============================================================" << endl;

	PDBReader rAv3;
	rAv3.open("/tmp/testPdb2.pdb");
	rAv3.read();
	AtomPointerVector av3 = rAv3.getAtomPointers();
	rAv3.close();

	// reshufle av3
	RandomNumberGenerator rgen;
	AtomPointerVector tmp;
	for (unsigned int i=0; i<av3.size(); i++) {
		tmp.push_back(av3[i]);
	}
	av3.clear();
	while (tmp.size()>0) {
		unsigned int rnum = rgen(tmp.size()-1);
		AtomPointerVector::iterator k=tmp.begin()+rnum;
		av3.push_back(*k);
		tmp.erase(k);
	}

	double RMSD2 = 0.0;
	if (tr.smartRmsdAlignment(av3, av1)) {
		RMSD2 = tr.getLastRMSD();
		cout << "Transform sucessful: RMSD = " << setprecision(15) << RMSD2 << endl;
		cout << "Rotation Matrix: " << tr.getLastRotationMatrix() << endl;
		cout << "Translation Vector: " << tr.getLastTranslation() <<endl;
	} else {
		cout << "Transform failed" << endl;
		result = false;
	}

	PDBWriter wAv3;
	wAv3.open("/tmp/testPdb3-aligned.pdb");
	wAv3.write(av3);
	wAv3.close();

	cout << endl;
	cout << "===================================================" << endl;
	cout << "Check the results against the expected values" << endl;
	cout << endl;

	if (av2[3]->getCoor().distance(c03) < epsilon) {
		cout << "Coordinate test 1 OK" << endl;
	} else {
		cout << "Coordinate test 1 NOT OK" << endl;
		result = false;
	}

	if (av2[25]->getCoor().distance(c25) < epsilon) {
		cout << "Coordinate test 2 OK" << endl;
	} else {
		cout << "Coordinate test 2 NOT OK" << endl;
		result = false;
	}

	if (av2[48]->getCoor().distance(c48) < epsilon) {
		cout << "Coordinate test 3 OK" << endl;
	} else {
		cout << "Coordinate test 3 NOT OK" << endl;
		result = false;
	}
	
	if (fabs(RMSD1 - 0.67384543920146) < epsilon) {
		cout << "RMSD test 1 OK" << endl;
	} else {
		cout << "RMSD test 1 NOT OK" << endl;
		result = false;
	}

	if (fabs(RMSD2 - 0.67384543920146) < epsilon) {
		cout << "RMSD test 2 OK" << endl;
	} else {
		cout << "RMSD test 2 NOT OK" << endl;
		result = false;
	}

	cout << endl;
	if (result) {
		cout << "GOLD" << endl;
	} else {
		cout << "LEAD" << endl;
	}

	return 0;



}
