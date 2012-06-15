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

#include "CharmmEEF1ParameterReader.h"

using namespace std;
using namespace MSL;

int main(int argc,char *argv[]) {

/*
! Parameters for solvation free energy calculation in EEF1
! Volume in A^3, energies in Kcal/mol, cp in cal/molK
! Sigw (A) is the distance within which 84% of Gsolv arises
!      Volume     Gref    Gfree      Href   CPref	Sigw
CHEX
C       14.7     0.000     0.00     0.000   0.00        3.5
CT1     23.7    -0.645    -0.86    -1.038   0.00        3.5
CT2     22.4    -0.720    -1.01    -1.160   0.00        3.5
CT3     30.0    -0.665    -0.92    -1.071   0.00        3.5
CC      14.7     0.000     0.00     0.000   0.00        3.5
NH1      4.4    -1.145    -1.72    -1.843   0.00        3.5
NR1      4.4    -1.145    -1.72    -1.843   0.00        3.5
NY       4.4    -1.145    -1.72    -1.843   0.00        3.5
NR2      4.4    -1.630    -1.71    -2.624   0.00        3.5
NH2     11.2    -1.145    -1.64    -1.843   0.00        3.5
NH3     11.2    -1.145    -1.15    -1.843   0.00        6.0
NC2     11.2    -0.200    -0.20    -0.322   0.00        6.0
N        0.0    -1.145    -1.72    -1.843   0.00        3.5
CPH1    18.4    -0.410    -0.57    -0.660   0.00        3.5
CPH2    18.4    -0.410    -0.57    -0.660   0.00        3.5
CA      18.4    -0.410    -0.57    -0.660   0.00        3.5
CY      18.4    -0.410    -0.57    -0.660   0.00        3.5
CPT     18.4    -0.410    -0.57    -0.660   0.00        3.5
CP1     23.7    -0.645    -0.86    -1.038   0.00        3.5
CP2     22.4    -0.720    -1.01    -1.160   0.00        3.5
CP3     22.4    -0.720    -1.01    -1.160   0.00        3.5
S       14.7    -1.780    -2.25    -2.866   0.00        3.5
SM      14.7    -1.780    -2.25    -2.866   0.00        3.5
OH1     10.8    -0.960    -1.09    -1.546   0.00        3.5
O       10.8    -1.270    -1.39    -2.045   0.00        3.5
OC      10.8    -0.900    -0.90    -1.449   0.00        6.0
END
WATER
C       14.7     0.000     0.00     0.000   0.00        3.5
CT1     23.7    -0.187    -0.25	    0.876   0.00 	3.5
CT2     22.4     0.372     0.52	   -0.610  18.60 	3.5
CT3     30.0     1.089     1.50	   -1.779  35.60 	3.5
CC      14.7     0.000     0.00	    0.000   0.00 	3.5
NH1      4.4    -5.950    -8.90    -9.059  -8.80        3.5
NR1      4.4    -5.950    -8.90    -9.059  -8.80 	3.5
NY       4.4    -5.950    -8.90    -9.059  -8.80 	3.5
NR2      4.4    -3.820    -4.00    -4.654  -8.80 	3.5
NH2     11.2    -5.450    -7.80    -9.028  -7.00        3.5
NH3     11.2   -20.000   -20.00   -25.000 -18.00        6.0
NC2     11.2   -10.000   -10.00   -12.000  -7.00        6.0
N        0.0    -1.000    -1.55    -1.250   8.80        3.5
CPH1    18.4     0.057     0.08    -0.973   6.90 	3.5
CPH2    18.4     0.057     0.08    -0.973   6.90 	3.5
CA      18.4     0.057     0.08    -0.973   6.90 	3.5
CY      18.4     0.057     0.08    -0.973   6.90 	3.5
CPT     18.4     0.057     0.08    -0.973   6.90 	3.5
CP1     23.7    -0.187    -0.25	    0.876   0.00 	3.5
CP2     22.4     0.372     0.52	   -0.610  18.60 	3.5
CP3     22.4     0.372     0.52	   -0.610  18.60 	3.5
S       14.7    -3.240    -4.10    -4.475 -39.90        3.5
SM      14.7    -3.240    -4.10    -4.475 -39.90 	3.5
OH1     10.8    -5.920    -6.70    -9.264 -11.20        3.5
O       10.8    -5.330    -5.85    -5.787  -8.80        3.5
OC      10.8   -10.000   -10.00   -12.000  -9.40        6.0
END
*/

	if (argc < 2) {
		cout << "Usage: testCharmmEEF1ParameterReader <parfile>" << endl;
		exit(0);
	}

	CharmmEEF1ParameterReader parRead(argv[1]);
	if (!parRead.read()) {
		cerr << "Cannot read " << argv[1] << endl;
		exit(1);
	}
	vector<string> solvents;
	solvents.push_back("WATER");
	solvents.push_back("CHEX");

	vector<string> atomTypes;
	atomTypes.push_back("C");
	atomTypes.push_back("CT1");
	atomTypes.push_back("CT2");
	atomTypes.push_back("CT3");
	atomTypes.push_back("CC");
	atomTypes.push_back("NH1");
	atomTypes.push_back("NR1");
	atomTypes.push_back("NY");
	atomTypes.push_back("NR2");
	atomTypes.push_back("NH2");
	atomTypes.push_back("NH3");
	atomTypes.push_back("NC2");
	atomTypes.push_back("N");
	atomTypes.push_back("CPH1");
	atomTypes.push_back("CPH2");
	atomTypes.push_back("CA");
	atomTypes.push_back("CY");
	atomTypes.push_back("CPT");
	atomTypes.push_back("CP1");
	atomTypes.push_back("CP2");
	atomTypes.push_back("CP3");
	atomTypes.push_back("S");
	atomTypes.push_back("SM");
	atomTypes.push_back("OH1");
	atomTypes.push_back("O");
	atomTypes.push_back("OC");

	for (unsigned int i=0; i<solvents.size(); i++) {
		cout << "Solvent: " << solvents[i] << endl;
		if (!parRead.solventExists(solvents[i])) {
			cerr << "Solvent " << solvents[i] << " not found" << endl;
			continue;
		}
		for (unsigned int j=0; j<atomTypes.size(); j++) {
			vector<double> params;
			if (parRead.EEF1Param(params, atomTypes[j], solvents[i])) {
				cout << "    " << solvents[i] << ", " << atomTypes[j] << " OK" << endl;
			} else {
				cout << "    " << solvents[i] << ", " << atomTypes[j] << " NOT OK" << endl;
				continue;
			}
			for (unsigned int k=0; k<params.size(); k++) {
				cout << params[k] << " ";
			}
			cout << endl;
		}
	}



	return 0;
}
