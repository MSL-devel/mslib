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

#include "AtomPointerVector.h"
#include "AtomSelection.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "Frame.h"

using namespace MSL;
using namespace std;


int main() {
	
		string pdbtext = "\
ATOM      1  N   ASP A   2       9.293 -10.362  -4.795  1.00 54.85           N  \n\
ATOM      2  CA  ASP A   2      10.085 -10.223  -3.534  1.00 54.00           C  \n\
ATOM      3  C   ASP A   2       9.493  -9.176  -2.574  1.00 51.85           C  \n\
ATOM      4  O   ASP A   2      10.074  -8.895  -1.518  1.00 52.87           O  \n\
ATOM      5  CB  ASP A   2      11.550  -9.870  -3.869  1.00 56.19           C  \n\
ATOM      6  CG  ASP A   2      11.679  -8.661  -4.815  1.00 58.62           C  \n\
ATOM      7  OD1 ASP A   2      12.827  -8.280  -5.144  1.00 58.14           O  \n\
ATOM      8  OD2 ASP A   2      10.642  -8.091  -5.229  1.00 52.27           O  \n\
ATOM      9  N   TYR A   3       8.335  -8.617  -2.937  1.00 48.51           N  \n\
ATOM     10  CA  TYR A   3       7.679  -7.595  -2.120  1.00 45.01           C  \n\
ATOM     11  C   TYR A   3       7.326  -8.017  -0.702  1.00 40.31           C  \n\
ATOM     12  O   TYR A   3       7.682  -7.324   0.244  1.00 37.49           O  \n\
ATOM     13  CB  TYR A   3       6.413  -7.073  -2.805  1.00 48.98           C  \n\
ATOM     14  CG  TYR A   3       6.661  -5.949  -3.785  1.00 49.13           C  \n\
ATOM     15  CD1 TYR A   3       6.532  -6.156  -5.161  1.00 53.20           C  \n\
ATOM     16  CD2 TYR A   3       7.022  -4.673  -3.340  1.00 43.03           C  \n\
ATOM     17  CE1 TYR A   3       6.754  -5.125  -6.070  1.00 52.04           C  \n\
ATOM     18  CE2 TYR A   3       7.250  -3.628  -4.247  1.00 49.44           C  \n\
ATOM     19  CZ  TYR A   3       7.113  -3.864  -5.611  1.00 52.69           C  \n\
ATOM     20  OH  TYR A   3       7.339  -2.857  -6.525  1.00 60.68           O  \n\
ATOM     21  N   LEU A   4       6.617  -9.133  -0.554  1.00 34.74           N  \n\
ATOM     22  CA  LEU A   4       6.238  -9.615   0.773  1.00 31.53           C  \n\
ATOM     23  C   LEU A   4       7.456  -9.817   1.679  1.00 32.30           C  \n\
ATOM     24  O   LEU A   4       7.406  -9.534   2.877  1.00 32.30           O  \n\
ATOM     25  CB  LEU A   4       5.438 -10.918   0.658  1.00 31.12           C  \n\
ATOM     26  CG  LEU A   4       3.921 -10.791   0.472  1.00 35.49           C  \n\
ATOM     27  CD1 LEU A   4       3.599  -9.670  -0.481  1.00 28.81           C  \n\
ATOM     28  CD2 LEU A   4       3.356 -12.099  -0.033  1.00 32.47           C  \n\
ATOM     29  N   ARG A   5       8.545 -10.309   1.097  1.00 30.31           N  \n\
ATOM     30  CA  ARG A   5       9.788 -10.531   1.825  1.00 27.67           C  \n\
ATOM     31  C   ARG A   5      10.305  -9.205   2.405  1.00 28.79           C  \n\
ATOM     32  O   ARG A   5      10.636  -9.110   3.595  1.00 24.11           O  \n\
ATOM     33  CB  ARG A   5      10.830 -11.119   0.871  1.00 28.46           C  \n\
ATOM     34  CG  ARG A   5      12.247 -11.145   1.404  1.00 18.83           C  \n\
ATOM     35  CD  ARG A   5      12.591 -12.498   1.987  1.00 26.30           C  \n\
ATOM     36  NE  ARG A   5      13.887 -12.461   2.653  1.00 32.24           N  \n\
ATOM     37  CZ  ARG A   5      14.330 -13.404   3.477  1.00 34.40           C  \n\
ATOM     38  NH1 ARG A   5      13.581 -14.466   3.733  1.00 33.47           N  \n\
ATOM     39  NH2 ARG A   5      15.515 -13.270   4.060  1.00 35.78           N  \n\
ATOM     40  N   GLU A   6      10.374  -8.182   1.556  1.00 26.41           N  \n\
ATOM     41  CA  GLU A   6      10.862  -6.878   1.983  1.00 25.86           C  \n\
ATOM     42  C   GLU A   6       9.923  -6.160   2.943  1.00 21.37           C  \n\
ATOM     43  O   GLU A   6      10.371  -5.436   3.836  1.00 18.92           O  \n\
ATOM     44  CB  GLU A   6      11.153  -6.009   0.760  1.00 26.67           C  \n\
ATOM     45  CG  GLU A   6      12.607  -6.060   0.350  1.00 34.50           C  \n\
ATOM     46  CD  GLU A   6      13.167  -7.471   0.377  1.00 51.21           C  \n\
ATOM     47  OE1 GLU A   6      12.812  -8.274  -0.518  1.00 56.62           O  \n\
ATOM     48  OE2 GLU A   6      13.955  -7.777   1.303  1.00 57.21           O  \n\
ATOM     49  N   LEU A   7       8.620  -6.362   2.763  1.00 21.07           N  \n\
ATOM     50  CA  LEU A   7       7.634  -5.753   3.644  1.00 16.71           C  \n\
ATOM     51  C   LEU A   7       7.818  -6.328   5.031  1.00 18.64           C  \n\
ATOM     52  O   LEU A   7       7.813  -5.606   6.026  1.00 22.41           O  \n\
ATOM     53  CB  LEU A   7       6.218  -6.058   3.159  1.00 17.01           C  \n\
ATOM     54  CG  LEU A   7       5.581  -5.167   2.083  1.00 17.50           C  \n\
ATOM     55  CD1 LEU A   7       6.596  -4.222   1.451  1.00 23.59           C  \n\
ATOM     56  CD2 LEU A   7       4.948  -6.069   1.037  1.00 16.21           C  \n\
ATOM     57  N   TYR A   8       7.995  -7.645   5.073  1.00 20.14           N  \n\
ATOM     58  CA  TYR A   8       8.165  -8.386   6.312  1.00 17.57           C  \n\
ATOM     59  C   TYR A   8       9.367  -7.895   7.132  1.00 20.68           C  \n\
ATOM     60  O   TYR A   8       9.241  -7.639   8.336  1.00 22.34           O  \n\
ATOM     61  CB  TYR A   8       8.301  -9.875   5.995  1.00 16.35           C  \n\
ATOM     62  CG  TYR A   8       7.919 -10.778   7.141  1.00 14.97           C  \n\
ATOM     63  CD1 TYR A   8       6.607 -11.234   7.289  1.00 20.76           C  \n\
ATOM     64  CD2 TYR A   8       8.872 -11.171   8.092  1.00 19.60           C  \n\
ATOM     65  CE1 TYR A   8       6.250 -12.070   8.355  1.00 15.92           C  \n\
ATOM     66  CE2 TYR A   8       8.532 -11.991   9.158  1.00 19.88           C  \n\
ATOM     67  CZ  TYR A   8       7.219 -12.446   9.288  1.00 16.54           C  \n\
ATOM     68  OH  TYR A   8       6.892 -13.294  10.327  1.00 23.44           O  \n\
ATOM     69  N   LYS A   9      10.526  -7.755   6.489  1.00 22.43           N  \n\
ATOM     70  CA  LYS A   9      11.719  -7.275   7.191  1.00 24.39           C  \n\
ATOM     71  C   LYS A   9      11.465  -5.880   7.763  1.00 26.85           C  \n\
ATOM     72  O   LYS A   9      11.766  -5.597   8.924  1.00 28.31           O  \n\
ATOM     73  CB  LYS A   9      12.923  -7.208   6.248  1.00 25.83           C  \n\
ATOM     74  CG  LYS A   9      13.436  -8.548   5.762  1.00 27.05           C  \n\
ATOM     75  CD  LYS A   9      14.817  -8.380   5.139  1.00 24.90           C  \n\
ATOM     76  CE  LYS A   9      15.334  -9.662   4.508  1.00 31.82           C  \n\
ATOM     77  NZ  LYS A   9      16.625  -9.436   3.790  1.00 30.93           N  \n\
END \n";


	AtomPointerVector av;
	stringstream ss;
	ss.str(pdbtext);
	PDBReader rAv(ss);
	rAv.read();
	av = rAv.getAtomPointers();
	cout << "Read atom vector with size " << av.size() << endl;
	rAv.close();



	AtomSelection sel(av);

	AtomPointerVector n   = sel.select("n,  resi 2 and chain A and name N");
	AtomPointerVector ca  = sel.select("ca, resi 2 and chain A and name CA");
	AtomPointerVector cb  = sel.select("cb, resi 2 and chain A and name CB");

	Frame f;
	f.computeFrameFrom3Atoms(n(0),ca(0),cb(0));
	
	PDBWriter w("ats.pdb");
	w.open();
	w.write(av);
	w.close();


	ofstream fout("frames.py");
	fout << &f <<endl;
	fout.close();
	


}


