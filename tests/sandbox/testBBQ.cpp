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

#include <vector>
#include <map>

#include "BBQTableReader.h"
#include "BBQTableWriter.h"
#include "BBQTable.h"
#include "PDBReader.h"
#include "System.h"
#include "Chain.h"
#include "Position.h"
#include "Residue.h"
#include "AtomSelection.h"

#include "testData.h"

using namespace std;

using namespace MSL;


// Taken from 1jlu.pdb
string pdbTestData = "\
ATOM      1  N   GLN E  12      -6.184 -17.144 -12.468  1.00 45.23           N  \n\
ATOM      2  CA  GLN E  12      -4.907 -16.995 -13.174  1.00 99.99           C  \n\
ATOM      3  C   GLN E  12      -5.077 -16.275 -14.492  1.00 99.99           C  \n\
ATOM      4  O   GLN E  12      -4.186 -15.617 -15.026  1.00 31.45           O  \n\
ATOM      5  CB  GLN E  12      -4.161 -18.307 -13.398  1.00 46.50           C  \n\
ATOM      6  N   GLU E  13      -6.246 -16.378 -15.044  1.00 99.99           N  \n\
ATOM      7  CA  GLU E  13      -6.415 -15.674 -16.274  1.00 20.81           C  \n\
ATOM      8  C   GLU E  13      -6.996 -14.311 -15.943  1.00 16.85           C  \n\
ATOM      9  O   GLU E  13      -6.885 -13.359 -16.706  1.00 33.39           O  \n\
ATOM     10  CB  GLU E  13      -7.480 -16.396 -17.125  1.00 99.99           C  \n\
ATOM     11  CG  GLU E  13      -6.988 -17.618 -17.918  1.00 33.84           C  \n\
ATOM     12  CD  GLU E  13      -5.889 -17.306 -18.924  1.00 99.99           C  \n\
ATOM     13  OE1 GLU E  13      -5.066 -18.183 -19.217  1.00 95.55           O  \n\
ATOM     14  OE2 GLU E  13      -5.853 -16.073 -19.462  1.00 36.15           O  \n\
ATOM     15  N   SER E  14      -7.662 -14.290 -14.781  1.00 33.83           N  \n\
ATOM     16  CA  SER E  14      -8.321 -13.124 -14.253  1.00 48.46           C  \n\
ATOM     17  C   SER E  14      -7.315 -12.021 -14.019  1.00 57.74           C  \n\
ATOM     18  O   SER E  14      -7.662 -10.850 -14.202  1.00 25.39           O  \n\
ATOM     19  CB  SER E  14      -9.109 -13.445 -12.979  1.00 34.26           C  \n\
ATOM     20  N   VAL E  15      -6.096 -12.471 -13.617  1.00 22.70           N  \n\
ATOM     21  CA  VAL E  15      -4.906 -11.688 -13.312  1.00 53.28           C  \n\
ATOM     22  C   VAL E  15      -4.381 -11.040 -14.578  1.00 25.24           C  \n\
ATOM     23  O   VAL E  15      -4.260  -9.807 -14.651  1.00 22.94           O  \n\
ATOM     24  CB  VAL E  15      -3.889 -12.516 -12.510  1.00 15.85           C  \n\
ATOM     25  CG1 VAL E  15      -2.560 -11.811 -12.271  1.00 26.27           C  \n\
ATOM     26  CG2 VAL E  15      -4.489 -12.895 -11.170  1.00 62.65           C  \n\
ATOM     27  N   LYS E  16      -4.175 -11.892 -15.578  1.00 33.67           N  \n\
ATOM     28  CA  LYS E  16      -3.734 -11.472 -16.903  1.00 54.17           C  \n\
ATOM     29  C   LYS E  16      -4.690 -10.412 -17.468  1.00  9.04           C  \n\
ATOM     30  O   LYS E  16      -4.305  -9.377 -18.027  1.00 43.47           O  \n\
ATOM     31  CB  LYS E  16      -3.612 -12.642 -17.862  1.00 27.02           C  \n\
ATOM     32  N   GLU E  17      -5.974 -10.654 -17.300  1.00 19.04           N  \n\
ATOM     33  CA  GLU E  17      -6.946  -9.704 -17.792  1.00 17.73           C  \n\
ATOM     34  C   GLU E  17      -6.950  -8.424 -16.980  1.00 20.40           C  \n\
ATOM     35  O   GLU E  17      -7.244  -7.300 -17.475  1.00 19.09           O  \n\
ATOM     36  CB  GLU E  17      -8.365 -10.328 -17.869  1.00 48.83           C  \n\
ATOM     37  N   PHE E  18      -6.634  -8.582 -15.707  1.00 22.72           N  \n\
ATOM     38  CA  PHE E  18      -6.614  -7.426 -14.857  1.00 38.84           C  \n\
ATOM     39  C   PHE E  18      -5.529  -6.506 -15.324  1.00 26.49           C  \n\
ATOM     40  O   PHE E  18      -5.727  -5.330 -15.645  1.00 17.12           O  \n\
ATOM     41  CB  PHE E  18      -6.372  -7.800 -13.392  1.00 19.58           C  \n\
ATOM     42  CG  PHE E  18      -6.156  -6.586 -12.523  1.00 67.16           C  \n\
ATOM     43  CD1 PHE E  18      -7.235  -5.902 -11.966  1.00 22.30           C  \n\
ATOM     44  CD2 PHE E  18      -4.874  -6.128 -12.233  1.00 37.88           C  \n\
ATOM     45  CE1 PHE E  18      -7.067  -4.801 -11.127  1.00 38.15           C  \n\
ATOM     46  CE2 PHE E  18      -4.674  -5.025 -11.405  1.00 42.76           C  \n\
ATOM     47  CZ  PHE E  18      -5.772  -4.360 -10.862  1.00 39.10           C  \n\
ATOM     48  N   LEU E  19      -4.365  -7.098 -15.336  1.00 25.19           N  \n\
ATOM     49  CA  LEU E  19      -3.179  -6.425 -15.755  1.00 37.85           C  \n\
ATOM     50  C   LEU E  19      -3.408  -5.654 -17.074  1.00 37.15           C  \n\
ATOM     51  O   LEU E  19      -3.047  -4.484 -17.222  1.00 37.91           O  \n\
ATOM     52  CB  LEU E  19      -2.074  -7.495 -15.861  1.00 10.25           C  \n\
ATOM     53  CG  LEU E  19      -1.065  -7.536 -14.700  1.00 30.67           C  \n\
ATOM     54  CD1 LEU E  19      -1.583  -6.841 -13.448  1.00 18.11           C  \n\
ATOM     55  CD2 LEU E  19      -0.617  -8.957 -14.416  1.00 10.09           C  \n\
ATOM     56  N   ALA E  20      -4.032  -6.320 -18.031  1.00 34.57           N  \n\
ATOM     57  CA  ALA E  20      -4.349  -5.736 -19.329  1.00 12.02           C  \n\
ATOM     58  C   ALA E  20      -5.082  -4.401 -19.248  1.00 18.86           C  \n\
ATOM     59  O   ALA E  20      -4.649  -3.389 -19.867  1.00 16.24           O  \n\
ATOM     60  CB  ALA E  20      -5.180  -6.733 -20.111  1.00 29.39           C  \n\
ATOM     61  N   LYS E  21      -6.207  -4.389 -18.501  1.00 16.22           N  \n\
ATOM     62  CA  LYS E  21      -6.902  -3.134 -18.382  1.00 16.85           C  \n\
ATOM     63  C   LYS E  21      -6.004  -2.161 -17.643  1.00 45.19           C  \n\
ATOM     64  O   LYS E  21      -5.862  -1.016 -18.027  1.00 43.70           O  \n\
ATOM     65  CB  LYS E  21      -8.304  -3.187 -17.792  1.00 27.63           C  \n\
ATOM     66  N   ALA E  22      -5.363  -2.613 -16.585  1.00 25.81           N  \n\
ATOM     67  CA  ALA E  22      -4.493  -1.726 -15.860  1.00 17.17           C  \n\
ATOM     68  C   ALA E  22      -3.403  -1.084 -16.707  1.00 23.40           C  \n\
ATOM     69  O   ALA E  22      -3.072   0.082 -16.500  1.00 17.75           O  \n\
ATOM     70  CB  ALA E  22      -3.873  -2.414 -14.694  1.00 13.88           C  \n\
ATOM     71  N   LYS E  23      -2.829  -1.828 -17.657  1.00  9.78           N  \n\
ATOM     72  CA  LYS E  23      -1.796  -1.235 -18.505  1.00  9.22           C  \n\
ATOM     73  C   LYS E  23      -2.412  -0.127 -19.333  1.00 49.13           C  \n\
ATOM     74  O   LYS E  23      -1.804   0.889 -19.579  1.00 23.33           O  \n\
ATOM     75  CB  LYS E  23      -1.161  -2.235 -19.444  1.00  5.60           C  \n\
ATOM     76  CG  LYS E  23       0.125  -1.730 -20.074  1.00 16.68           C  \n\
ATOM     77  CD  LYS E  23       0.995  -2.839 -20.669  1.00 43.24           C  \n\
ATOM     78  CE  LYS E  23       2.286  -2.304 -21.290  1.00 41.28           C  \n\
ATOM     79  NZ  LYS E  23       3.509  -2.747 -20.613  1.00 22.33           N  \n\
ATOM     80  N   GLU E  24      -3.641  -0.332 -19.781  1.00 20.62           N  \n\
ATOM     81  CA  GLU E  24      -4.284   0.679 -20.589  1.00 29.83           C  \n\
ATOM     82  C   GLU E  24      -4.457   1.978 -19.826  1.00 36.66           C  \n\
ATOM     83  O   GLU E  24      -4.155   3.065 -20.318  1.00 30.43           O  \n\
ATOM     84  CB  GLU E  24      -5.577   0.173 -21.315  1.00 44.07           C  \n\
ATOM     85  N   ASP E  25      -4.936   1.865 -18.601  1.00 21.17           N  \n\
ATOM     86  CA  ASP E  25      -5.100   3.044 -17.800  1.00 20.70           C  \n\
ATOM     87  C   ASP E  25      -3.747   3.769 -17.647  1.00 13.66           C  \n\
ATOM     88  O   ASP E  25      -3.593   4.952 -17.947  1.00 27.83           O  \n\
ATOM     89  CB  ASP E  25      -5.843   2.752 -16.491  1.00 30.22           C  \n\
ATOM     90  N   PHE E  26      -2.735   3.047 -17.211  1.00 25.02           N  \n\
ATOM     91  CA  PHE E  26      -1.433   3.632 -17.041  1.00 22.43           C  \n\
ATOM     92  C   PHE E  26      -0.969   4.360 -18.272  1.00  5.18           C  \n\
ATOM     93  O   PHE E  26      -0.517   5.491 -18.185  1.00 22.06           O  \n\
ATOM     94  CB  PHE E  26      -0.424   2.537 -16.764  1.00 20.40           C  \n\
ATOM     95  CG  PHE E  26       1.015   2.967 -16.853  1.00 24.95           C  \n\
ATOM     96  CD1 PHE E  26       1.587   3.700 -15.808  1.00 26.88           C  \n\
ATOM     97  CD2 PHE E  26       1.817   2.577 -17.933  1.00 16.81           C  \n\
ATOM     98  CE1 PHE E  26       2.923   4.100 -15.850  1.00 14.40           C  \n\
ATOM     99  CE2 PHE E  26       3.164   2.951 -17.988  1.00 24.00           C  \n\
ATOM    100  CZ  PHE E  26       3.704   3.710 -16.942  1.00 19.17           C  \n\
ATOM    101  N   LEU E  27      -1.067   3.680 -19.418  1.00 14.45           N  \n\
ATOM    102  CA  LEU E  27      -0.623   4.306 -20.647  1.00 21.66           C  \n\
ATOM    103  C   LEU E  27      -1.328   5.598 -20.958  1.00 21.45           C  \n\
ATOM    104  O   LEU E  27      -0.700   6.525 -21.440  1.00 32.95           O  \n\
ATOM    105  CB  LEU E  27      -0.540   3.419 -21.892  1.00 15.60           C  \n\
ATOM    106  CG  LEU E  27       0.282   2.174 -21.697  1.00 58.43           C  \n\
ATOM    107  CD1 LEU E  27       0.010   1.224 -22.872  1.00 19.84           C  \n\
ATOM    108  CD2 LEU E  27       1.743   2.565 -21.632  1.00 29.78           C  \n\
ATOM    109  N   LYS E  28      -2.626   5.703 -20.686  1.00 20.88           N  \n\
ATOM    110  CA  LYS E  28      -3.273   6.976 -20.992  1.00 30.58           C  \n\
ATOM    111  C   LYS E  28      -2.579   8.100 -20.247  1.00 23.42           C  \n\
ATOM    112  O   LYS E  28      -2.256   9.136 -20.822  1.00 24.59           O  \n\
ATOM    113  CB  LYS E  28      -4.771   6.943 -20.744  1.00 36.86           C  \n\
ATOM    114  N   LYS E  29      -2.325   7.866 -18.960  1.00 17.86           N  \n\
ATOM    115  CA  LYS E  29      -1.644   8.828 -18.096  1.00 16.12           C  \n\
ATOM    116  C   LYS E  29      -0.157   8.952 -18.389  1.00 30.84           C  \n\
ATOM    117  O   LYS E  29       0.379  10.025 -18.213  1.00 19.71           O  \n\
ATOM    118  CB  LYS E  29      -1.841   8.460 -16.650  1.00 16.49           C  \n\
ATOM    119  CG  LYS E  29      -3.339   8.348 -16.365  1.00 13.84           C  \n\
ATOM    120  CD  LYS E  29      -3.634   8.284 -14.903  1.00 25.76           C  \n\
ATOM    121  CE  LYS E  29      -4.966   7.600 -14.601  1.00 54.63           C  \n\
ATOM    122  NZ  LYS E  29      -4.698   6.369 -13.856  1.00 99.99           N  \n\
ATOM    123  N   TRP E  30       0.504   7.877 -18.823  1.00 13.80           N  \n\
ATOM    124  CA  TRP E  30       1.939   7.927 -19.125  1.00  9.86           C  \n\
ATOM    125  C   TRP E  30       2.185   8.817 -20.304  1.00  7.56           C  \n\
ATOM    126  O   TRP E  30       3.135   9.585 -20.315  1.00 15.06           O  \n\
ATOM    127  CB  TRP E  30       2.575   6.544 -19.337  1.00 51.85           C  \n\
ATOM    128  CG  TRP E  30       4.001   6.548 -19.831  1.00 43.91           C  \n\
ATOM    129  CD1 TRP E  30       4.384   6.354 -21.119  1.00 20.48           C  \n\
ATOM    130  CD2 TRP E  30       5.234   6.699 -19.080  1.00 26.59           C  \n\
ATOM    131  NE1 TRP E  30       5.744   6.346 -21.237  1.00 30.73           N  \n\
ATOM    132  CE2 TRP E  30       6.296   6.556 -20.008  1.00 26.31           C  \n\
ATOM    133  CE3 TRP E  30       5.550   6.929 -17.739  1.00 26.34           C  \n\
ATOM    134  CZ2 TRP E  30       7.631   6.661 -19.628  1.00 20.43           C  \n\
ATOM    135  CZ3 TRP E  30       6.869   7.060 -17.359  1.00 39.32           C  \n\
ATOM    136  CH2 TRP E  30       7.897   6.905 -18.291  1.00 25.52           C  \n\
ATOM    137  N   GLU E  31       1.321   8.737 -21.292  1.00 20.83           N  \n\
ATOM    138  CA  GLU E  31       1.492   9.569 -22.475  1.00 47.17           C  \n\
ATOM    139  C   GLU E  31       0.952  10.989 -22.327  1.00 17.65           C  \n\
ATOM    140  O   GLU E  31       1.322  11.894 -23.056  1.00 38.89           O  \n\
ATOM    141  CB  GLU E  31       1.103   8.828 -23.752  1.00 29.52           C  \n\
ATOM    142  N   THR E  32       0.103  11.233 -21.337  1.00 58.21           N  \n\
ATOM    143  CA  THR E  32      -0.411  12.594 -21.163  1.00 16.99           C  \n\
ATOM    144  C   THR E  32      -0.544  13.028 -19.668  1.00 14.72           C  \n\
ATOM    145  O   THR E  32      -1.646  13.035 -19.092  1.00 35.18           O  \n\
ATOM    146  CB  THR E  32      -1.636  12.854 -22.082  1.00 34.37           C  \n\
ATOM    147  OG1 THR E  32      -2.449  13.883 -21.556  1.00 88.48           O  \n\
ATOM    148  CG2 THR E  32      -2.470  11.570 -22.263  1.00 29.57           C  \n\
ATOM    149  N   PRO E  33       0.617  13.383 -19.039  1.00 24.28           N  \n\
ATOM    150  CA  PRO E  33       0.674  13.752 -17.631  1.00 39.43           C  \n\
ATOM    151  C   PRO E  33       0.160  15.086 -17.160  1.00 99.99           C  \n\
ATOM    152  O   PRO E  33       0.087  16.069 -17.893  1.00 14.38           O  \n\
ATOM    153  CB  PRO E  33       2.105  13.528 -17.127  1.00 27.57           C  \n\
ATOM    154  CG  PRO E  33       2.925  13.125 -18.329  1.00 19.63           C  \n\
ATOM    155  CD  PRO E  33       2.001  13.086 -19.535  1.00 21.95           C  \n\
ATOM    156  N   SER E  34      -0.145  15.039 -15.862  1.00 33.63           N  \n\
ATOM    157  CA  SER E  34      -0.655  16.126 -15.061  1.00 18.42           C  \n\
ATOM    158  C   SER E  34       0.415  17.124 -14.739  1.00 96.78           C  \n\
ATOM    159  O   SER E  34       1.525  16.690 -14.465  1.00 54.82           O  \n\
ATOM    160  CB  SER E  34      -1.006  15.531 -13.743  1.00 35.97           C  \n\
ATOM    161  OG  SER E  34      -2.385  15.408 -13.717  1.00 29.56           O  \n\
ATOM    162  N   GLN E  35       0.069  18.417 -14.719  1.00 14.61           N  \n\
ATOM    163  CA  GLN E  35       1.018  19.479 -14.405  1.00 50.74           C  \n\
ATOM    164  C   GLN E  35       0.379  20.624 -13.622  1.00 41.51           C  \n\
ATOM    165  O   GLN E  35      -0.734  21.034 -13.929  1.00 32.21           O  \n\
ATOM    166  CB  GLN E  35       1.545  20.084 -15.698  1.00 27.51           C  \n\
ATOM    167  CG  GLN E  35       2.695  19.327 -16.342  1.00 38.12           C  \n\
ATOM    168  CD  GLN E  35       3.321  20.221 -17.401  1.00 36.36           C  \n\
ATOM    169  OE1 GLN E  35       2.611  20.715 -18.295  1.00 99.99           O  \n\
ATOM    170  NE2 GLN E  35       4.636  20.457 -17.307  1.00 58.37           N  \n\
ATOM    171  N   ASN E  36       1.098  21.134 -12.630  1.00 40.13           N  \n\
ATOM    172  CA  ASN E  36       0.700  22.242 -11.745  1.00 27.60           C  \n\
ATOM    173  C   ASN E  36      -0.753  22.234 -11.288  1.00 73.69           C  \n\
ATOM    174  O   ASN E  36      -1.610  22.991 -11.745  1.00 30.71           O  \n\
ATOM    175  CB  ASN E  36       1.373  23.607 -12.062  1.00 64.56           C  \n\
ATOM    176  CG  ASN E  36       2.910  23.574 -11.821  1.00 99.99           C  \n\
ATOM    177  OD1 ASN E  36       3.483  24.446 -11.137  1.00 39.67           O  \n\
ATOM    178  ND2 ASN E  36       3.599  22.571 -12.401  1.00 53.02           N  \n\
ATOM    179  N   THR E  37      -1.037  21.337 -10.344  1.00 22.66           N  \n\
ATOM    180  CA  THR E  37      -2.390  21.168  -9.823  1.00 53.58           C  \n\
ATOM    181  C   THR E  37      -2.655  21.769  -8.437  1.00 21.98           C  \n\
ATOM    182  O   THR E  37      -3.726  21.619  -7.841  1.00 44.07           O  \n\
ATOM    183  CB  THR E  37      -2.814  19.698  -9.941  1.00 21.15           C  \n\
ATOM    184  OG1 THR E  37      -1.839  18.865  -9.339  1.00 45.15           O  \n\
ATOM    185  CG2 THR E  37      -2.872  19.344 -11.414  1.00 15.72           C  \n\
ATOM    186  N   ALA E  38      -1.658  22.449  -7.926  1.00 28.45           N  \n\
ATOM    187  CA  ALA E  38      -1.776  23.054  -6.643  1.00 30.56           C  \n\
ATOM    188  C   ALA E  38      -0.764  24.172  -6.527  1.00 21.88           C  \n\
ATOM    189  O   ALA E  38      -0.075  24.505  -7.505  1.00 20.21           O  \n\
ATOM    190  CB  ALA E  38      -1.701  22.033  -5.520  1.00 21.83           C  \n\
ATOM    191  N   GLN E  39      -0.736  24.761  -5.320  1.00 18.68           N  \n\
ATOM    192  CA  GLN E  39       0.135  25.871  -4.993  1.00  7.72           C  \n\
ATOM    193  C   GLN E  39       0.562  25.661  -3.591  1.00 16.20           C  \n\
ATOM    194  O   GLN E  39      -0.124  25.045  -2.783  1.00 30.98           O  \n\
ATOM    195  CB  GLN E  39      -0.664  27.165  -4.934  1.00 52.29           C  \n\
ATOM    196  CG  GLN E  39      -0.467  28.125  -6.107  1.00 63.37           C  \n\
ATOM    197  CD  GLN E  39      -1.732  28.928  -6.326  1.00 22.65           C  \n\
ATOM    198  OE1 GLN E  39      -2.435  29.285  -5.372  1.00 57.31           O  \n\
ATOM    199  NE2 GLN E  39      -2.054  29.194  -7.585  1.00 99.99           N  \n\
ATOM    200  N   LEU E  40       1.704  26.212  -3.319  1.00 17.38           N  \n\
ATOM    201  CA  LEU E  40       2.275  26.104  -2.026  1.00 24.67           C  \n\
ATOM    202  C   LEU E  40       1.377  26.589  -0.921  1.00 37.40           C  \n\
ATOM    203  O   LEU E  40       1.253  25.965   0.120  1.00 33.99           O  \n\
ATOM    204  CB  LEU E  40       3.559  26.924  -1.989  1.00 23.74           C  \n\
ATOM    205  CG  LEU E  40       4.280  26.683  -0.682  1.00 53.59           C  \n\
ATOM    206  CD1 LEU E  40       4.256  25.178  -0.315  1.00 27.52           C  \n\
ATOM    207  CD2 LEU E  40       5.711  27.183  -0.766  1.00 11.00           C  \n\
ATOM    208  N   ASP E  41       0.784  27.745  -1.170  1.00 25.67           N  \n\
ATOM    209  CA  ASP E  41      -0.103  28.435  -0.247  1.00 21.65           C  \n\
ATOM    210  C   ASP E  41      -1.395  27.726   0.114  1.00 23.82           C  \n\
ATOM    211  O   ASP E  41      -2.001  27.965   1.152  1.00 42.42           O  \n\
ATOM    212  CB  ASP E  41      -0.414  29.820  -0.786  1.00 23.24           C  \n\
ATOM    213  CG  ASP E  41      -1.231  29.741  -2.032  1.00 35.84           C  \n\
ATOM    214  OD1 ASP E  41      -2.313  29.173  -2.097  1.00 99.99           O  \n\
ATOM    215  OD2 ASP E  41      -0.656  30.361  -3.030  1.00 99.99           O  \n\
ATOM    216  N   GLN E  42      -1.819  26.861  -0.749  1.00 25.17           N  \n\
ATOM    217  CA  GLN E  42      -3.027  26.120  -0.505  1.00 19.26           C  \n\
ATOM    218  C   GLN E  42      -2.853  25.132   0.673  1.00 18.41           C  \n\
ATOM    219  O   GLN E  42      -3.768  24.387   1.027  1.00 20.38           O  \n\
ATOM    220  CB  GLN E  42      -3.222  25.288  -1.785  1.00 13.77           C  \n\
ATOM    221  CG  GLN E  42      -3.567  26.091  -3.026  1.00 17.91           C  \n\
ATOM    222  CD  GLN E  42      -4.152  25.215  -4.128  1.00 25.76           C  \n\
ATOM    223  OE1 GLN E  42      -4.847  24.204  -3.875  1.00 99.99           O  \n\
ATOM    224  NE2 GLN E  42      -3.866  25.594  -5.365  1.00 51.20           N  \n\
ATOM    225  N   PHE E  43      -1.662  25.123   1.262  1.00 20.24           N  \n\
ATOM    226  CA  PHE E  43      -1.348  24.196   2.359  1.00 40.92           C  \n\
ATOM    227  C   PHE E  43      -0.598  24.846   3.532  1.00 16.36           C  \n\
ATOM    228  O   PHE E  43       0.040  25.899   3.396  1.00 44.92           O  \n\
ATOM    229  CB  PHE E  43      -0.431  23.008   1.889  1.00 31.78           C  \n\
ATOM    230  CG  PHE E  43      -0.708  22.380   0.552  1.00 30.47           C  \n\
ATOM    231  CD1 PHE E  43      -1.533  21.258   0.460  1.00 13.65           C  \n\
ATOM    232  CD2 PHE E  43      -0.136  22.899  -0.613  1.00 35.27           C  \n\
ATOM    233  CE1 PHE E  43      -1.824  20.691  -0.782  1.00 14.54           C  \n\
ATOM    234  CE2 PHE E  43      -0.401  22.333  -1.861  1.00 36.87           C  \n\
ATOM    235  CZ  PHE E  43      -1.251  21.229  -1.937  1.00 11.36           C  \n\
ATOM    236  N   ASP E  44      -0.627  24.152   4.678  1.00 42.67           N  \n\
ATOM    237  CA  ASP E  44       0.082  24.578   5.875  1.00 30.35           C  \n\
ATOM    238  C   ASP E  44       1.165  23.560   6.194  1.00 30.79           C  \n\
ATOM    239  O   ASP E  44       0.943  22.356   6.045  1.00 41.30           O  \n\
ATOM    240  CB  ASP E  44      -0.841  24.788   7.068  1.00 51.67           C  \n\
ATOM    241  CG  ASP E  44      -1.658  26.004   6.835  1.00 64.52           C  \n\
ATOM    242  OD1 ASP E  44      -1.184  27.016   6.331  1.00 64.64           O  \n\
ATOM    243  OD2 ASP E  44      -2.903  25.830   7.200  1.00 41.43           O  \n\
ATOM    244  N   ARG E  45       2.330  24.035   6.600  1.00 36.47           N  \n\
ATOM    245  CA  ARG E  45       3.387  23.111   6.878  1.00 52.69           C  \n\
ATOM    246  C   ARG E  45       3.347  22.709   8.311  1.00 25.36           C  \n\
ATOM    247  O   ARG E  45       3.725  23.469   9.182  1.00 70.02           O  \n\
ATOM    248  CB  ARG E  45       4.753  23.646   6.472  1.00 99.99           C  \n\
ATOM    249  N   ILE E  46       2.901  21.494   8.558  1.00 49.22           N  \n\
ATOM    250  CA  ILE E  46       2.807  20.982   9.895  1.00  8.57           C  \n\
ATOM    251  C   ILE E  46       4.129  20.545  10.537  1.00 35.28           C  \n\
ATOM    252  O   ILE E  46       4.589  21.167  11.491  1.00 52.39           O  \n\
ATOM    253  CB  ILE E  46       1.777  19.876   9.951  1.00 49.42           C  \n\
ATOM    254  CG1 ILE E  46       0.442  20.417   9.474  1.00 17.67           C  \n\
ATOM    255  CG2 ILE E  46       1.681  19.317  11.363  1.00 27.63           C  \n\
ATOM    256  CD1 ILE E  46      -0.062  21.639  10.198  1.00 38.28           C  \n\
ATOM    257  N   LYS E  47       4.725  19.457  10.022  1.00 33.87           N  \n\
ATOM    258  CA  LYS E  47       5.966  18.909  10.549  1.00 36.34           C  \n\
ATOM    259  C   LYS E  47       6.772  18.194   9.497  1.00 26.47           C  \n\
ATOM    260  O   LYS E  47       6.189  17.681   8.585  1.00 21.28           O  \n\
ATOM    261  CB  LYS E  47       5.645  17.847  11.586  1.00 26.74           C  \n\
ATOM    262  CG  LYS E  47       6.868  17.413  12.389  1.00 99.99           C  \n\
ATOM    263  CD  LYS E  47       6.490  16.527  13.574  1.00 62.42           C  \n\
ATOM    264  CE  LYS E  47       7.587  16.416  14.630  1.00 38.45           C  \n\
ATOM    265  NZ  LYS E  47       7.230  15.534  15.759  1.00 99.99           N  \n\
ATOM    266  N   THR E  48       8.107  18.137   9.637  1.00 13.37           N  \n\
ATOM    267  CA  THR E  48       8.923  17.407   8.685  1.00 16.62           C  \n\
ATOM    268  C   THR E  48       8.787  15.936   9.059  1.00 35.22           C  \n\
ATOM    269  O   THR E  48       8.743  15.654  10.239  1.00 27.02           O  \n\
ATOM    270  CB  THR E  48      10.367  17.882   8.753  1.00  7.52           C  \n\
ATOM    271  OG1 THR E  48      10.441  19.132   8.143  1.00 29.85           O  \n\
ATOM    272  CG2 THR E  48      11.298  16.891   8.058  1.00 40.63           C  \n\
ATOM    273  N   LEU E  49       8.683  14.988   8.085  1.00 21.43           N  \n\
ATOM    274  CA  LEU E  49       8.512  13.544   8.419  1.00 18.31           C  \n\
ATOM    275  C   LEU E  49       9.713  12.721   8.015  1.00 20.67           C  \n\
ATOM    276  O   LEU E  49       9.844  11.547   8.344  1.00 28.06           O  \n\
ATOM    277  CB  LEU E  49       7.291  12.890   7.708  1.00 14.70           C  \n\
ATOM    278  CG  LEU E  49       5.950  13.506   8.067  1.00  3.82           C  \n\
ATOM    279  CD1 LEU E  49       4.908  13.307   6.987  1.00  5.22           C  \n\
ATOM    280  CD2 LEU E  49       5.398  12.968   9.363  1.00 20.42           C  \n\
ATOM    281  N   GLY E  50      10.595  13.340   7.267  1.00 31.72           N  \n\
ATOM    282  CA  GLY E  50      11.729  12.591   6.842  1.00 15.92           C  \n\
ATOM    283  C   GLY E  50      12.602  13.478   6.023  1.00 21.85           C  \n\
ATOM    284  O   GLY E  50      12.132  14.477   5.498  1.00 24.77           O  \n\
TER     285                                                                     \n\
END                                                                             \n";

void getAtomsOfInterest(map<string, bool> &atomsOfInterest);
void getChainWithJustCA(vector<Residue *> &_resVec, vector<Position *> &_posVec);
void testWriter(System &sys, string fileName);
void testReader(System &sys, string fileName);
void writeResidueVecToPDB(PDBWriter &_pdbWriter, vector<Residue *> &_resVec);

#include "SysEnv.h"
static SysEnv SYSENV;

int main() {
    PDBReader pdb;
    System sys;

    pdb.read(pdbTestData);
    sys.addAtoms(pdb.getAtomPointers());

    testWriter(sys, "bbqTableTest.txt");
    testReader(sys, "bbqTableTest.txt");


    // Test using a system/chain approach
    System justCAsys;
    AtomSelection sel(sys.getAtomPointers());
    AtomPointerVector chA = sel.select("name CA");	

    justCAsys.addAtoms(chA);

    justCAsys.writePdb("/tmp/ChainBBQ.test.pre.pdb");
    cout << "BBQ TABLE: "<<SYSENV.getEnv("MSL_BBQ_TABLE")<<endl;
    BBQTable bbqTable(SYSENV.getEnv("MSL_BBQ_TABLE"));
    cout << "CHAIN E SIZE1: "<<justCAsys("E").atomSize()<<endl;


    bbqTable.fillInMissingBBAtoms(justCAsys("E"));

    cout << "CHAIN E SIZE2: "<<justCAsys("E").getAtomPointers().size()<<endl;

    justCAsys.writePdb("/tmp/ChainBBQ.test.post.pdb");
};

void testWriter(System &sys, string fileName) {
    vector<Chain *> allChains;
    BBQTable bbqTable;
    BBQTableWriter bbqTableWriter;
    
    map<string, bool> atomsOfInterest;

    getAtomsOfInterest(atomsOfInterest);    
    
    allChains = sys.getChains();

    bbqTable.setBinSizes(0.02f, 0.02f, 0.02f);

    for (unsigned int chainIndex = 0; chainIndex < allChains.size(); ++chainIndex) {
        Chain *currChain = allChains[chainIndex];
        vector<Position *> posVec = currChain->getPositions();
        vector<Residue *> resVec;

        for (unsigned int positionIndex = 0; positionIndex < posVec.size(); ++positionIndex) {
            Residue &currRes = posVec[positionIndex]->getCurrentIdentity();
            resVec.push_back(&currRes);
        }

        //bbqTable.addQuadrilateralInfoFromResidues(resVec, atomsOfInterest);
	bbqTable.addQuadrilateralInfoFromResidues(resVec);
    }

    bbqTable.normalize();

    bbqTableWriter.open(fileName);
    bbqTableWriter.write(bbqTable);
    bbqTableWriter.close();
}

void testReader(System &sys, string fileName) {
    BBQTableReader bbqTableReader(fileName);
    BBQTable bbqTable;
    PDBWriter pdbWriter(fileName + ".pdb");
    vector<Chain *> allChains = sys.getChains();
    vector<Residue *> resVec;

    bbqTableReader.open();
    bbqTableReader.read(bbqTable);
    bbqTableReader.close();
    pdbWriter.open();

    for (unsigned int chainIndex = 0; chainIndex < allChains.size(); ++chainIndex) {
        Chain *currChain = allChains[chainIndex];
        Chain newChain;
        vector<Position *> posVec = currChain->getPositions();
        vector<Residue *> resVec;
        
        getChainWithJustCA(resVec, posVec);
        bbqTable.fillInMissingBBAtoms(resVec);
        
        writeResidueVecToPDB(pdbWriter, resVec);
    }
    
    pdbWriter.close();

    cout << "Done.\n";
}

void writeResidueVecToPDB(PDBWriter &_pdbWriter, vector<Residue *> &_resVec) {
    for(unsigned int index = 0; index < _resVec.size(); ++index) {
        _pdbWriter.write( _resVec[index]->getAtomPointers(), (index == (_resVec.size()))?true:false );
        
        delete _resVec[index];
    }
}

void getChainWithJustCA(vector<Residue *> &_resVec, vector<Position *> &_posVec) {
    for (unsigned int positionIndex = 0; positionIndex < _posVec.size(); ++positionIndex) {
        Residue &currRes = _posVec[positionIndex]->getCurrentIdentity();
        Residue *newRes = new Residue();
        newRes->addAtom( currRes.getAtom("CA") );
        _resVec.push_back(newRes);
    }
}

void getAtomsOfInterest(map<string, bool> &atomsOfInterest) {
    atomsOfInterest["N"] = true;
    atomsOfInterest["C"] = true;
    atomsOfInterest["O"] = true;
}
