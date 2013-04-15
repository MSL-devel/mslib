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

#include "System.h"
#include "EZpotentialBuilder.h"
#include "SysEnv.h"

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

int main() {

	string inputPdb = "/tmp/orientedHelix.pdb";

	string pdbContent = "ATOM      1  N   LEU A   1      12.279   0.388  16.404  1.00  0.00           N  \n\
ATOM      2  CA  LEU A   1      12.266  -0.576  15.321  1.00  0.00           C  \n\
ATOM      3  CB  LEU A   1      13.167  -1.778  15.719  1.00  0.00           C  \n\
ATOM      4  CG  LEU A   1      14.680  -1.508  15.896  1.00  0.00           C  \n\
ATOM      5  CD1 LEU A   1      15.376  -2.729  16.515  1.00  0.00           C  \n\
ATOM      6  CD2 LEU A   1      15.366  -1.116  14.579  1.00  0.00           C  \n\
ATOM      7  C   LEU A   1      10.847  -1.009  14.995  1.00  0.00           C  \n\
ATOM      8  O   LEU A   1      10.484  -1.134  13.826  1.00  0.00           O  \n\
ATOM      9  N   GLY A   2      10.028  -1.243  16.037  1.00  0.00           N  \n\
ATOM     10  CA  GLY A   2       8.653  -1.638  15.825  1.00  0.00           C  \n\
ATOM     11  C   GLY A   2       7.895  -0.573  15.096  1.00  0.00           C  \n\
ATOM     12  O   GLY A   2       7.110  -0.862  14.196  1.00  0.00           O  \n\
ATOM     13  N   THR A   3       8.119   0.701  15.475  1.00  0.00           N  \n\
ATOM     14  CA  THR A   3       7.444   1.807  14.809  1.00  0.00           C  \n\
ATOM     15  CB  THR A   3       7.678   3.169  15.553  1.00  0.00           C  \n\
ATOM     16  OG1 THR A   3       7.183   3.143  16.889  1.00  0.00           O  \n\
ATOM     17  CG2 THR A   3       7.102   4.413  14.868  1.00  0.00           C  \n\
ATOM     18  C   THR A   3       7.820   1.868  13.342  1.00  0.00           C  \n\
ATOM     19  O   THR A   3       6.962   2.095  12.491  1.00  0.00           O  \n\
ATOM     20  N   LEU A   4       9.114   1.664  13.037  1.00  0.00           N  \n\
ATOM     21  CA  LEU A   4       9.565   1.683  11.657  1.00  0.00           C  \n\
ATOM     22  CB  LEU A   4      11.115   1.561  11.643  1.00  0.00           C  \n\
ATOM     23  CG  LEU A   4      11.927   2.721  12.266  1.00  0.00           C  \n\
ATOM     24  CD1 LEU A   4      13.409   2.337  12.384  1.00  0.00           C  \n\
ATOM     25  CD2 LEU A   4      11.778   4.033  11.481  1.00  0.00           C  \n\
ATOM     26  C   LEU A   4       8.898   0.581  10.852  1.00  0.00           C  \n\
ATOM     27  O   LEU A   4       8.490   0.799   9.712  1.00  0.00           O  \n\
ATOM     28  N   GLY A   5       8.783  -0.621  11.446  1.00  0.00           N  \n\
ATOM     29  CA  GLY A   5       8.142  -1.725  10.765  1.00  0.00           C  \n\
ATOM     30  C   GLY A   5       6.709  -1.412  10.466  1.00  0.00           C  \n\
ATOM     31  O   GLY A   5       6.212  -1.711   9.382  1.00  0.00           O  \n\
ATOM     32  N   TYR A   6       6.006  -0.795  11.436  1.00  0.00           N  \n\
ATOM     33  CA  TYR A   6       4.611  -0.431  11.229  1.00  0.00           C  \n\
ATOM     34  CB  TYR A   6       4.026   0.118  12.567  1.00  0.00           C  \n\
ATOM     35  CG  TYR A   6       4.098  -0.884  13.696  1.00  0.00           C  \n\
ATOM     36  CD1 TYR A   6       4.388  -0.463  15.007  1.00  0.00           C  \n\
ATOM     37  CD2 TYR A   6       3.949  -2.259  13.439  1.00  0.00           C  \n\
ATOM     38  CE1 TYR A   6       4.547  -1.404  16.035  1.00  0.00           C  \n\
ATOM     39  CE2 TYR A   6       4.112  -3.198  14.468  1.00  0.00           C  \n\
ATOM     40  CZ  TYR A   6       4.415  -2.769  15.763  1.00  0.00           C  \n\
ATOM     41  OH  TYR A   6       4.607  -3.706  16.794  1.00  0.00           O  \n\
ATOM     42  C   TYR A   6       4.471   0.565  10.084  1.00  0.00           C  \n\
ATOM     43  O   TYR A   6       3.564   0.453   9.264  1.00  0.00           O  \n\
ATOM     44  N   VAL A   7       5.383   1.557  10.025  1.00  0.00           N  \n\
ATOM     45  CA  VAL A   7       5.343   2.540   8.959  1.00  0.00           C  \n\
ATOM     46  CB  VAL A   7       6.418   3.650   9.210  1.00  0.00           C  \n\
ATOM     47  CG1 VAL A   7       6.467   4.736   8.113  1.00  0.00           C  \n\
ATOM     48  CG2 VAL A   7       6.236   4.338  10.577  1.00  0.00           C  \n\
ATOM     49  C   VAL A   7       5.538   1.882   7.605  1.00  0.00           C  \n\
ATOM     50  O   VAL A   7       4.857   2.221   6.638  1.00  0.00           O  \n\
ATOM     51  N   LEU A   8       6.482   0.924   7.524  1.00  0.00           N  \n\
ATOM     52  CA  LEU A   8       6.724   0.225   6.280  1.00  0.00           C  \n\
ATOM     53  CB  LEU A   8       7.954  -0.706   6.470  1.00  0.00           C  \n\
ATOM     54  CG  LEU A   8       9.322  -0.035   6.739  1.00  0.00           C  \n\
ATOM     55  CD1 LEU A   8      10.369  -1.087   7.135  1.00  0.00           C  \n\
ATOM     56  CD2 LEU A   8       9.821   0.787   5.541  1.00  0.00           C  \n\
ATOM     57  C   LEU A   8       5.491  -0.540   5.837  1.00  0.00           C  \n\
ATOM     58  O   LEU A   8       5.147  -0.543   4.656  1.00  0.00           O  \n\
ATOM     59  N   GLY A   9       4.810  -1.202   6.791  1.00  0.00           N  \n\
ATOM     60  CA  GLY A   9       3.609  -1.942   6.467  1.00  0.00           C  \n\
ATOM     61  C   GLY A   9       2.547  -1.033   5.932  1.00  0.00           C  \n\
ATOM     62  O   GLY A   9       1.858  -1.366   4.971  1.00  0.00           O  \n\
ATOM     63  N   ILE A  10       2.391   0.153   6.553  1.00  0.00           N  \n\
ATOM     64  CA  ILE A  10       1.398   1.113   6.090  1.00  0.00           C  \n\
ATOM     65  CB  ILE A  10       1.291   2.308   7.100  1.00  0.00           C  \n\
ATOM     66  CG2 ILE A  10       0.462   3.500   6.573  1.00  0.00           C  \n\
ATOM     67  CG1 ILE A  10       0.858   1.880   8.525  1.00  0.00           C  \n\
ATOM     68  CD1 ILE A  10       1.002   2.987   9.583  1.00  0.00           C  \n\
ATOM     69  C   ILE A  10       1.696   1.562   4.670  1.00  0.00           C  \n\
ATOM     70  O   ILE A  10       0.789   1.682   3.848  1.00  0.00           O  \n\
ATOM     71  N   THR A  11       2.984   1.816   4.370  1.00  0.00           N  \n\
ATOM     72  CA  THR A  11       3.367   2.229   3.038  1.00  0.00           C  \n\
ATOM     73  CB  THR A  11       4.868   2.683   2.974  1.00  0.00           C  \n\
ATOM     74  OG1 THR A  11       5.130   3.794   3.828  1.00  0.00           O  \n\
ATOM     75  CG2 THR A  11       5.403   3.033   1.581  1.00  0.00           C  \n\
ATOM     76  C   THR A  11       3.036   1.154   2.020  1.00  0.00           C  \n\
ATOM     77  O   THR A  11       2.548   1.458   0.933  1.00  0.00           O  \n\
ATOM     78  N   MET A  12       3.302  -0.117   2.371  1.00  0.00           N  \n\
ATOM     79  CA  MET A  12       3.000  -1.215   1.470  1.00  0.00           C  \n\
ATOM     80  CB  MET A  12       3.543  -2.533   2.091  1.00  0.00           C  \n\
ATOM     81  CG  MET A  12       5.073  -2.600   2.302  1.00  0.00           C  \n\
ATOM     82  SD  MET A  12       5.498  -3.981   3.411  1.00  0.00           S  \n\
ATOM     83  CE  MET A  12       5.215  -5.364   2.261  1.00  0.00           C  \n\
ATOM     84  C   MET A  12       1.508  -1.295   1.190  1.00  0.00           C  \n\
ATOM     85  O   MET A  12       1.096  -1.516   0.054  1.00  0.00           O  \n\
ATOM     86  N   MET A  13       0.682  -1.111   2.238  1.00  0.00           N  \n\
ATOM     87  CA  MET A  13      -0.754  -1.147   2.064  1.00  0.00           C  \n\
ATOM     88  CB  MET A  13      -1.428  -1.044   3.461  1.00  0.00           C  \n\
ATOM     89  CG  MET A  13      -1.117  -2.192   4.448  1.00  0.00           C  \n\
ATOM     90  SD  MET A  13      -1.605  -1.721   6.139  1.00  0.00           S  \n\
ATOM     91  CE  MET A  13      -3.412  -1.857   5.960  1.00  0.00           C  \n\
ATOM     92  C   MET A  13      -1.218  -0.043   1.129  1.00  0.00           C  \n\
ATOM     93  O   MET A  13      -2.074  -0.263   0.275  1.00  0.00           O  \n\
ATOM     94  N   VAL A  14      -0.646   1.167   1.286  1.00  0.00           N  \n\
ATOM     95  CA  VAL A  14      -1.008   2.274   0.428  1.00  0.00           C  \n\
ATOM     96  CB  VAL A  14      -0.306   3.586   0.917  1.00  0.00           C  \n\
ATOM     97  CG1 VAL A  14      -0.616   4.824   0.049  1.00  0.00           C  \n\
ATOM     98  CG2 VAL A  14      -0.645   3.914   2.384  1.00  0.00           C  \n\
ATOM     99  C   VAL A  14      -0.665   1.974  -1.021  1.00  0.00           C  \n\
ATOM    100  O   VAL A  14      -1.444   2.272  -1.925  1.00  0.00           O  \n\
ATOM    101  N   ILE A  15       0.519   1.375  -1.254  1.00  0.00           N  \n\
ATOM    102  CA  ILE A  15       0.922   1.031  -2.600  1.00  0.00           C  \n\
ATOM    103  CB  ILE A  15       2.409   0.533  -2.607  1.00  0.00           C  \n\
ATOM    104  CG2 ILE A  15       2.852  -0.103  -3.944  1.00  0.00           C  \n\
ATOM    105  CG1 ILE A  15       3.421   1.593  -2.102  1.00  0.00           C  \n\
ATOM    106  CD1 ILE A  15       4.836   1.042  -1.854  1.00  0.00           C  \n\
ATOM    107  C   ILE A  15      -0.039   0.030  -3.218  1.00  0.00           C  \n\
ATOM    108  O   ILE A  15      -0.401   0.150  -4.387  1.00  0.00           O  \n\
ATOM    109  N   ILE A  16      -0.464  -0.973  -2.427  1.00  0.00           N  \n\
ATOM    110  CA  ILE A  16      -1.397  -1.962  -2.920  1.00  0.00           C  \n\
ATOM    111  CB  ILE A  16      -1.567  -3.115  -1.872  1.00  0.00           C  \n\
ATOM    112  CG2 ILE A  16      -2.732  -4.081  -2.185  1.00  0.00           C  \n\
ATOM    113  CG1 ILE A  16      -0.250  -3.871  -1.562  1.00  0.00           C  \n\
ATOM    114  CD1 ILE A  16      -0.340  -4.820  -0.355  1.00  0.00           C  \n\
ATOM    115  C   ILE A  16      -2.717  -1.317  -3.310  1.00  0.00           C  \n\
ATOM    116  O   ILE A  16      -3.299  -1.656  -4.339  1.00  0.00           O  \n\
ATOM    117  N   ILE A  17      -3.201  -0.374  -2.480  1.00  0.00           N  \n\
ATOM    118  CA  ILE A  17      -4.440   0.311  -2.779  1.00  0.00           C  \n\
ATOM    119  CB  ILE A  17      -4.869   1.205  -1.565  1.00  0.00           C  \n\
ATOM    120  CG2 ILE A  17      -6.032   2.173  -1.879  1.00  0.00           C  \n\
ATOM    121  CG1 ILE A  17      -5.113   0.405  -0.260  1.00  0.00           C  \n\
ATOM    122  CD1 ILE A  17      -5.276   1.281   0.993  1.00  0.00           C  \n\
ATOM    123  C   ILE A  17      -4.331   1.090  -4.079  1.00  0.00           C  \n\
ATOM    124  O   ILE A  17      -5.257   1.088  -4.889  1.00  0.00           O  \n\
ATOM    125  N   ALA A  18      -3.187   1.769  -4.287  1.00  0.00           N  \n\
ATOM    126  CA  ALA A  18      -2.984   2.523  -5.505  1.00  0.00           C  \n\
ATOM    127  CB  ALA A  18      -1.660   3.317  -5.432  1.00  0.00           C  \n\
ATOM    128  C   ALA A  18      -3.009   1.599  -6.735  1.00  0.00           C  \n\
ATOM    129  O   ALA A  18      -3.584   1.903  -7.779  1.00  0.00           O  \n\
ATOM    130  N   ILE A  19      -2.361   0.413  -6.624  1.00  0.00           N  \n\
ATOM    131  CA  ILE A  19      -2.348  -0.539  -7.733  1.00  0.00           C  \n\
ATOM    132  CB  ILE A  19      -1.372  -1.724  -7.414  1.00  0.00           C  \n\
ATOM    133  CG2 ILE A  19      -1.482  -2.910  -8.399  1.00  0.00           C  \n\
ATOM    134  CG1 ILE A  19       0.097  -1.279  -7.200  1.00  0.00           C  \n\
ATOM    135  CD1 ILE A  19       1.011  -2.378  -6.633  1.00  0.00           C  \n\
ATOM    136  C   ILE A  19      -3.755  -1.005  -8.066  1.00  0.00           C  \n\
ATOM    137  O   ILE A  19      -4.115  -1.119  -9.237  1.00  0.00           O  \n\
ATOM    138  N   GLY A  20      -4.567  -1.279  -7.028  1.00  0.00           N  \n\
ATOM    139  CA  GLY A  20      -5.931  -1.708  -7.248  1.00  0.00           C  \n\
ATOM    140  C   GLY A  20      -6.718  -0.652  -7.959  1.00  0.00           C  \n\
ATOM    141  O   GLY A  20      -7.495  -0.947  -8.863  1.00  0.00           O  \n\
ATOM    142  N   ALA A  21      -6.529   0.621  -7.557  1.00  0.00           N  \n\
ATOM    143  CA  ALA A  21      -7.234   1.719  -8.203  1.00  0.00           C  \n\
ATOM    144  CB  ALA A  21      -6.950   3.049  -7.467  1.00  0.00           C  \n\
ATOM    145  C   ALA A  21      -6.855   1.818  -9.692  1.00  0.00           C  \n\
ATOM    146  O   ALA A  21      -7.685   2.034 -10.573  1.00  0.00           O  \n\
ATOM    147  N   GLY A  22      -5.545   1.654 -10.002  1.00  0.00           N  \n\
ATOM    148  CA  GLY A  22      -5.092   1.709 -11.390  1.00  0.00           C  \n\
ATOM    149  C   GLY A  22      -5.720   0.620 -12.203  1.00  0.00           C  \n\
ATOM    150  O   GLY A  22      -6.136   0.840 -13.338  1.00  0.00           O  \n\
ATOM    151  N   ILE A  23      -5.802  -0.597 -11.629  1.00  0.00           N  \n\
ATOM    152  CA  ILE A  23      -6.416  -1.714 -12.334  1.00  0.00           C  \n\
ATOM    153  CB  ILE A  23      -6.203  -3.042 -11.529  1.00  0.00           C  \n\
ATOM    154  CG2 ILE A  23      -7.035  -4.234 -12.054  1.00  0.00           C  \n\
ATOM    155  CG1 ILE A  23      -4.712  -3.415 -11.330  1.00  0.00           C  \n\
ATOM    156  CD1 ILE A  23      -4.479  -4.560 -10.329  1.00  0.00           C  \n\
ATOM    157  C   ILE A  23      -7.879  -1.431 -12.632  1.00  0.00           C  \n\
ATOM    158  O   ILE A  23      -8.363  -1.726 -13.723  1.00  0.00           O  \n\
ATOM    159  N   ILE A  24      -8.598  -0.850 -11.653  1.00  0.00           N  \n\
ATOM    160  CA  ILE A  24      -9.993  -0.523 -11.852  1.00  0.00           C  \n\
ATOM    161  CB  ILE A  24     -10.633  -0.046 -10.503  1.00  0.00           C  \n\
ATOM    162  CG2 ILE A  24     -12.041   0.572 -10.661  1.00  0.00           C  \n\
ATOM    163  CG1 ILE A  24     -10.589  -1.115  -9.382  1.00  0.00           C  \n\
ATOM    164  CD1 ILE A  24     -10.969  -0.582  -7.990  1.00  0.00           C  \n\
ATOM    165  C   ILE A  24     -10.160   0.485 -12.976  1.00  0.00           C  \n\
ATOM    166  O   ILE A  24     -11.066   0.361 -13.799  1.00  0.00           O  \n\
ATOM    167  N   LEU A  25      -9.277   1.501 -13.017  1.00  0.00           N  \n\
ATOM    168  CA  LEU A  25      -9.344   2.498 -14.063  1.00  0.00           C  \n\
ATOM    169  CB  LEU A  25      -8.295   3.604 -13.758  1.00  0.00           C  \n\
ATOM    170  CG  LEU A  25      -8.506   4.455 -12.483  1.00  0.00           C  \n\
ATOM    171  CD1 LEU A  25      -7.270   5.325 -12.208  1.00  0.00           C  \n\
ATOM    172  CD2 LEU A  25      -9.766   5.330 -12.560  1.00  0.00           C  \n\
ATOM    173  C   LEU A  25      -9.131   1.869 -15.429  1.00  0.00           C  \n\
ATOM    174  O   LEU A  25      -9.822   2.208 -16.388  1.00  0.00           O  \n\
ATOM    175  N   GLY A  26      -8.162   0.940 -15.526  1.00  0.00           N  \n\
ATOM    176  CA  GLY A  26      -7.901   0.270 -16.781  1.00  0.00           C  \n\
ATOM    177  C   GLY A  26      -9.096  -0.510 -17.232  1.00  0.00           C  \n\
ATOM    178  O   GLY A  26      -9.444  -0.505 -18.411  1.00  0.00           O  \n\
ATOM    179  N   TYR A  27      -9.759  -1.207 -16.288  1.00  0.00           N  \n\
ATOM    180  CA  TYR A  27     -10.946  -1.979 -16.628  1.00  0.00           C  \n\
ATOM    181  CB  TYR A  27     -11.395  -2.795 -15.376  1.00  0.00           C  \n\
ATOM    182  CG  TYR A  27     -10.331  -3.746 -14.879  1.00  0.00           C  \n\
ATOM    183  CD1 TYR A  27     -10.148  -3.957 -13.501  1.00  0.00           C  \n\
ATOM    184  CD2 TYR A  27      -9.459  -4.377 -15.786  1.00  0.00           C  \n\
ATOM    185  CE1 TYR A  27      -9.099  -4.766 -13.039  1.00  0.00           C  \n\
ATOM    186  CE2 TYR A  27      -8.409  -5.183 -15.323  1.00  0.00           C  \n\
ATOM    187  CZ  TYR A  27      -8.229  -5.373 -13.950  1.00  0.00           C  \n\
ATOM    188  OH  TYR A  27      -7.162  -6.156 -13.477  1.00  0.00           O  \n\
ATOM    189  C   TYR A  27     -12.052  -1.075 -17.156  1.00  0.00           C  \n\
ATOM    190  OXT TYR A  27     -11.916   0.142 -17.265  1.00  0.00           O  \n\
ATOM    191  O   TYR A  27     -13.149  -1.507 -17.507  1.00  0.00           O  \n\
TER     192                                                                     \n\
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
		cerr << "Cannot read PDB file " << inputPdb << endl;
		exit(1);
	}

	// test the N-terminal C-terminal function created in Atom in support of the EZ potential
	cout << "Test the isPositionNterminal and isPositionCterminal functions (used in EZ)" << endl;
	AtomPointerVector atoms = sys.getAtomPointers();
	for (unsigned int i=0; i<atoms.size(); i++) {
		if (atoms[i]->getName() == "CA") {
			cout << atoms[i]->getPositionId();
			if (atoms[i]->isPositionNterminal() && atoms[i]->isPositionCterminal()) {
				cout << " is BOTH N and C terminal" << endl;
			} else if (atoms[i]->isPositionNterminal()) {
				cout << " is N terminal" << endl;
			} else if (atoms[i]->isPositionCterminal()) {
				cout << " is C terminal" << endl;
			} else {
				cout << " is not terminal" << endl;
			}
		}
	}
	cout << endl;

	EZpotentialBuilder EZbuild(sys);
	EZbuild.setAddTermini(true);
	EZbuild.buildInteractions();
	
	double E = sys.calcEnergy();

	cout << "The EZ energy of " << inputPdb << " is " << E << endl;
	cout << endl;
	cout << "Summary:" << endl;
	cout << sys.getEnergySummary();

	cout << endl;
	cout << "Test result:" << endl;
	if (fabs(E - -8.646074) > 0.000001) {
		cout << "LEAD" << endl;
	} else {
		cout << "GOLD" << endl;
	}

}


