/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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

#include "CoiledCoils.h"
#include "PDBWriter.h"
#include "PDBReader.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;


int main() {

	// OFFERS FORMULATION
	CoiledCoils cc;
	cc.useBBQTable("/home/dwkulp/software/msl/tables/PiscesBBQTable.txt");
	cc.offersCoiledCoil(6.5, 1.5115, 194.0, 2.25, 30,1.79);
	AtomPointerVector offers = cc.getAtomPointers();


	PDBWriter pout("offers.pdb");
	pout.open();
	pout.write(offers);
	pout.close();

	System &offersSys  = cc.getSystem();
	offersSys.writePdb("offersSys.pdb");
	

	// CRICKS FORMULATION
	cc.cricksCoiledCoil(6.5,1.5115, 194.0, 2.25, 30,1.79);
	AtomPointerVector cricks = cc.getAtomPointers();

	pout.open("cricks.pdb");
	pout.write(cricks);
	pout.close();

	System &cricksSys  = cc.getSystem();
	offersSys.writePdb("cricksSys.pdb");

	
	// BEN NORTHS FORMULATION
	cc.northCoiledCoils(6.5,1.5115,194,2.25, 30,102, 0  );
	AtomPointerVector norths = cc.getAtomPointers();
	
	pout.open("norths.pdb");
	pout.write(norths);
	pout.close();


	// Coiling an ideal, Z-aligned helix
	string idealString = "ATOM     56  N   ALA A   1      -0.570   1.463 -14.265  1.00  0.00              \n\
ATOM     57  CA  ALA A   1      -1.716   1.560 -13.361  1.00  0.00              \n\
ATOM     58  C   ALA A   1      -1.670   0.475 -12.312  1.00  0.00              \n\
ATOM     59  O   ALA A   1      -1.991   0.687 -11.137  1.00  0.00              \n\
ATOM     60  CB  ALA A   1      -2.992   1.508 -14.219  1.00  0.00              \n\
ATOM     61  N   ALA A   2      -1.282  -0.718 -12.730  1.00  0.00              \n\
ATOM     62  CA  ALA A   2      -1.190  -1.860 -11.825  1.00  0.00              \n\
ATOM     63  C   ALA A   2      -0.126  -1.638 -10.776  1.00  0.00              \n\
ATOM     64  O   ALA A   2      -0.283  -1.987  -9.600  1.00  0.00              \n\
ATOM     65  CB  ALA A   2      -0.930  -3.112 -12.682  1.00  0.00              \n\
ATOM     66  N   ALA A   3       0.984  -1.060 -11.194  1.00  0.00              \n\
ATOM     67  CA  ALA A   3       2.096  -0.781 -10.289  1.00  0.00              \n\
ATOM     68  C   ALA A   3       1.704   0.229  -9.239  1.00  0.00              \n\
ATOM     69  O   ALA A   3       2.073   0.132  -8.065  1.00  0.00              \n\
ATOM     70  CB  ALA A   3       3.289  -0.320 -11.145  1.00  0.00              \n\
ATOM     71  N   ALA A   4       0.952   1.232  -9.655  1.00  0.00              \n\
ATOM     72  CA  ALA A   4       0.493   2.285  -8.752  1.00  0.00              \n\
ATOM     73  C   ALA A   4      -0.440   1.730  -7.703  1.00  0.00              \n\
ATOM     74  O   ALA A   4      -0.405   2.111  -6.527  1.00  0.00              \n\
ATOM     75  CB  ALA A   4      -0.157   3.384  -9.609  1.00  0.00              \n\
ATOM     76  N   ALA A   5      -1.305   0.822  -8.120  1.00  0.00              \n\
ATOM     77  CA  ALA A   5      -2.268   0.199  -7.217  1.00  0.00              \n\
ATOM     78  C   ALA A   5      -1.568  -0.630  -6.167  1.00  0.00              \n\
ATOM     79  O   ALA A   5      -1.950  -0.659  -4.991  1.00  0.00              \n\
ATOM     80  CB  ALA A   5      -3.248  -0.622  -8.073  1.00  0.00              \n\
ATOM     81  N   ALA A   6      -0.532  -1.334  -6.583  1.00  0.00              \n\
ATOM     82  CA  ALA A   6       0.241  -2.183  -5.680  1.00  0.00              \n\
ATOM     83  C   ALA A   6       0.943  -1.356  -4.631  1.00  0.00              \n\
ATOM     84  O   ALA A   6       1.034  -1.729  -3.455  1.00  0.00              \n\
ATOM     85  CB  ALA A   6       1.213  -3.012  -6.536  1.00  0.00              \n\
ATOM     86  N   ALA A   7       1.468  -0.217  -5.046  1.00  0.00              \n\
ATOM     87  CA  ALA A   7       2.179   0.684  -4.144  1.00  0.00              \n\
ATOM     88  C   ALA A   7       1.247   1.243  -3.095  1.00  0.00              \n\
ATOM     89  O   ALA A   7       1.600   1.393  -1.918  1.00  0.00              \n\
ATOM     90  CB  ALA A   7       2.837   1.779  -5.000  1.00  0.00              \n\
ATOM     91  N   ALA A   8       0.038   1.573  -3.511  1.00  0.00              \n\
ATOM     92  CA  ALA A   8      -0.966   2.125  -2.608  1.00  0.00              \n\
ATOM     93  C   ALA A   8      -1.365   1.113  -1.559  1.00  0.00              \n\
ATOM     94  O   ALA A   8      -1.573   1.438  -0.383  1.00  0.00              \n\
ATOM     95  CB  ALA A   8      -2.154   2.597  -3.464  1.00  0.00              \n\
ATOM     96  N   ALA A   9      -1.493  -0.132  -1.974  1.00  0.00              \n\
ATOM     97  CA  ALA A   9      -1.872  -1.214  -1.072  1.00  0.00              \n\
ATOM     98  C   ALA A   9      -0.812  -1.442  -0.023  1.00  0.00              \n\
ATOM     99  O   ALA A   9      -1.096  -1.698   1.154  1.00  0.00              \n\
ATOM    100  CB  ALA A   9      -2.143  -2.464  -1.927  1.00  0.00              \n\
ATOM    101  N   ALA A  10       0.438  -1.363  -0.437  1.00  0.00              \n\
ATOM    102  CA  ALA A  10       1.568  -1.561   0.466  1.00  0.00              \n\
ATOM    103  C   ALA A  10       1.617  -0.475   1.514  1.00  0.00              \n\
ATOM    104  O   ALA A  10       1.918  -0.715   2.690  1.00  0.00              \n\
ATOM    105  CB  ALA A  10       2.846  -1.621  -0.391  1.00  0.00              \n\
ATOM    106  N   ALA A  11       1.337   0.746   1.098  1.00  0.00              \n\
ATOM    107  CA  ALA A  11       1.345   1.893   2.002  1.00  0.00              \n\
ATOM    108  C   ALA A  11       0.265   1.764   3.051  1.00  0.00              \n\
ATOM    109  O   ALA A  11       0.453   2.097   4.228  1.00  0.00              \n\
ATOM    110  CB  ALA A  11       1.196   3.163   1.146  1.00  0.00              \n\
ATOM    111  N   ALA A  12      -0.891   1.283   2.635  1.00  0.00              \n\
ATOM    112  CA  ALA A  12      -2.023   1.105   3.538  1.00  0.00              \n\
ATOM    113  C   ALA A  12      -1.720   0.062   4.586  1.00  0.00              \n\
ATOM    114  O   ALA A  12      -2.082   0.192   5.763  1.00  0.00              \n\
ATOM    115  CB  ALA A  12      -3.250   0.751   2.680  1.00  0.00              \n\
ATOM    116  N   ALA A  13      -1.057  -1.002   4.171  1.00  0.00              \n\
ATOM    117  CA  ALA A  13      -0.696  -2.090   5.074  1.00  0.00              \n\
ATOM    118  C   ALA A  13       0.285  -1.618   6.123  1.00  0.00              \n\
ATOM    119  O   ALA A  13       0.215  -1.994   7.300  1.00  0.00              \n\
ATOM    120  CB  ALA A  13      -0.145  -3.243   4.218  1.00  0.00              \n\
ATOM    121  N   ALA A  14       1.223  -0.790   5.707  1.00  0.00              \n\
ATOM    122  CA  ALA A  14       2.238  -0.255   6.611  1.00  0.00              \n\
ATOM    123  C   ALA A  14       1.612   0.635   7.659  1.00  0.00              \n\
ATOM    124  O   ALA A  14       1.994   0.628   8.836  1.00  0.00              \n\
ATOM    125  CB  ALA A  14       3.285   0.480   5.754  1.00  0.00              \n\
ATOM    126  N   ALA A  15       0.642   1.428   7.243  1.00  0.00              \n\
ATOM    127  CA  ALA A  15      -0.054   2.337   8.147  1.00  0.00              \n\
ATOM    128  C   ALA A  15      -0.828   1.577   9.195  1.00  0.00              \n\
ATOM    129  O   ALA A  15      -0.885   1.955  10.371  1.00  0.00              \n\
ATOM    130  CB  ALA A  15      -0.950   3.250   7.290  1.00  0.00              \n\
ATOM    131  N   ALA A  16      -1.450   0.488   8.780  1.00  0.00              \n\
ATOM    132  CA  ALA A  16      -2.235  -0.346   9.683  1.00  0.00              \n\
ATOM    133  C   ALA A  16      -1.356  -0.986  10.732  1.00  0.00              \n\
ATOM    134  O   ALA A  16      -1.722  -1.104  11.909  1.00  0.00              \n\
ATOM    135  CB  ALA A  16      -2.988  -1.379   8.827  1.00  0.00              \n\
ATOM    136  N   ALA A  17      -0.182  -1.420  10.315  1.00  0.00              \n\
ATOM    137  CA  ALA A  17       0.772  -2.060  11.220  1.00  0.00              \n\
ATOM    138  C   ALA A  17       1.257  -1.085  12.268  1.00  0.00              \n\
ATOM    139  O   ALA A  17       1.432  -1.425  13.445  1.00  0.00              \n\
ATOM    140  CB  ALA A  17       1.912  -2.632  10.363  1.00  0.00              \n\
ATOM    141  N   ALA A  18       1.494   0.141  11.853  1.00  0.00              \n\
ATOM    142  CA  ALA A  18       1.966   1.188  12.757  1.00  0.00              \n\
ATOM    143  C   ALA A  18       0.928   1.507  13.804  1.00  0.00              \n\
ATOM    144  O   ALA A  18       1.233   1.738  14.982  1.00  0.00              \n\
ATOM    145  CB  ALA A  18       2.343   2.408  11.901  1.00  0.00              \n\
ATOM    146  N   ALA A  19      -0.323   1.540  13.390  1.00  0.00              \n\
ATOM    147  CA  ALA A  19      -1.432   1.834  14.293  1.00  0.00              \n\
ATOM    148  C   ALA A  19      -1.579   0.757  15.342  1.00  0.00              \n\
ATOM    149  O   ALA A  19      -1.857   1.023  16.517  1.00  0.00              \n\
ATOM    150  CB  ALA A  19      -2.699   2.007  13.435  1.00  0.00              \n\
ATOM    151  N   ALA A  20      -1.403  -0.484  14.925  1.00  0.00              \n\
ATOM    152  CA  ALA A  20      -1.514  -1.623  15.829  1.00  0.00              \n\
ATOM    153  C   ALA A  20      -0.428  -1.590  16.878  1.00  0.00              \n\
ATOM    154  O   ALA A  20      -0.642  -1.908  18.053  1.00  0.00              \n\
ATOM    155  CB  ALA A  20      -1.477  -2.903  14.973  1.00  0.00              \n\
ATOM    156  N   ALA A  21       0.765  -1.215  16.461  1.00  0.00              \n\
ATOM    157  CA  ALA A  21       1.911  -1.137  17.366  1.00  0.00              \n\
ATOM    158  C   ALA A  21       1.700  -0.070  18.414  1.00  0.00              \n\
ATOM    159  O   ALA A  21       2.048  -0.233  19.591  1.00  0.00              \n\
ATOM    160  CB  ALA A  21       3.166  -0.890  16.509  1.00  0.00              \n\
ATOM    161  N   ALA A  22       1.134   1.047  17.999  1.00  0.00              \n\
ATOM    162  CA  ALA A  22       0.868   2.163  18.902  1.00  0.00              \n\
ATOM    163  C   ALA A  22      -0.150   1.779  19.950  1.00  0.00              \n\
ATOM    164  O   ALA A  22      -0.049   2.149  21.126  1.00  0.00              \n\
ATOM    165  CB  ALA A  22       0.419   3.360  18.045  1.00  0.00              \n\
ATOM    166  N   ALA A  23      -1.158   1.038  19.534  1.00  0.00              \n\
ATOM    167  CA  ALA A  23      -2.217   0.592  20.438  1.00  0.00              \n\
ATOM    168  C   ALA A  23      -1.670  -0.348  21.486  1.00  0.00              \n\
ATOM    169  O   ALA A  23      -2.052  -0.309  22.662  1.00  0.00              \n\
ATOM    170  CB  ALA A  23      -3.321  -0.048  19.581  1.00  0.00              \n\
ATOM    171  N   ALA A  24      -0.774  -1.223  21.071  1.00  0.00              \n\
ATOM    172  CA  ALA A  24      -0.159  -2.191  21.974  1.00  0.00              \n\
ATOM    173  C   ALA A  24       0.676  -1.500  23.023  1.00  0.00              \n\
ATOM    174  O   ALA A  24       0.701  -1.884  24.200  1.00  0.00              \n\
ATOM    175  CB  ALA A  24       0.652  -3.179  21.118  1.00  0.00              \n\
ATOM    176  N   ALA A  25       1.391  -0.473  22.608  1.00  0.00              \n\
ATOM    177  CA  ALA A  25       2.246   0.291  23.511  1.00  0.00              \n\
ATOM    178  C   ALA A  25       1.429   1.005  24.559  1.00  0.00              \n\
ATOM    179  O   ALA A  25       1.801   1.092  25.735  1.00  0.00              \n\
ATOM    180  CB  ALA A  25       3.088   1.255  22.654  1.00  0.00              \n\
ATOM    181  N   ALA A  26       0.296   1.541  24.143  1.00  0.00              \n\
ATOM    182  CA  ALA A  26      -0.596   2.261  25.046  1.00  0.00              \n\
ATOM    183  C   ALA A  26      -1.167   1.335  26.095  1.00  0.00              \n\
ATOM    184  O   ALA A  26      -1.313   1.689  27.271  1.00  0.00              \n\
ATOM    185  CB  ALA A  26      -1.686   2.931  24.190  1.00  0.00              \n\
ATOM    186  N   ALA A  27      -1.510   0.129  25.679  1.00  0.00              \n\
ATOM    187  CA  ALA A  27      -2.073  -0.869  26.583  1.00  0.00              \n\
ATOM    188  C   ALA A  27      -1.067  -1.279  27.632  1.00  0.00              \n\
ATOM    189  O   ALA A  27      -1.390  -1.481  28.808  1.00  0.00              \n\
ATOM    190  CB  ALA A  27      -2.555  -2.052  25.727  1.00  0.00              \n\
ATOM    191  N   ALA A  28       0.179  -1.418  27.216  1.00  0.00              \n\
ATOM    192  CA  ALA A  28       1.256  -1.811  28.119  1.00  0.00              \n\
ATOM    193  C   ALA A  28       1.494  -0.751  29.168  1.00  0.00              \n\
ATOM    194  O   ALA A  28       1.748  -1.040  30.345  1.00  0.00              \n\
ATOM    195  CB  ALA A  28       2.503  -2.094  27.263  1.00  0.00              \n\
ATOM    196  N   ALA A  29       1.430   0.499  28.754  1.00  0.00              \n\
ATOM    197  CA  ALA A  29       1.638   1.626  29.656  1.00  0.00              \n\
ATOM    198  C   ALA A  29       0.554   1.688  30.705  1.00  0.00              \n\
ATOM    199  O   ALA A  29       0.797   1.986  31.881  1.00  0.00              \n\
ATOM    200  CB  ALA A  29       1.714   2.904  28.798  1.00  0.00              \n\
END";
	AtomPointerVector ideal;
	stringstream ss;
	ss.str(idealString);
	PDBReader pin(ss);
	pin.read();
	ideal = pin.getAtomPointers();
	pin.close();
	
	Transforms tr;
	CartesianPoint x(6.5,0,0);
	//ideal.translate(x);
	tr.translate(ideal, x);



	cc.applyCoiledCoil(ideal,194);


	pout.open("idealCoiled.pdb");
	pout.write(cc.getAtomPointers());
	pout.close();

	System &idealSys  = cc.getSystem();
	offersSys.writePdb("idealSys.pdb");

}
