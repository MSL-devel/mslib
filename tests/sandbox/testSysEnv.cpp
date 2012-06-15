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

#include "SysEnv.h"
#include <iostream>

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

/*
  Run this program twice to test. First just run with no variables set. Then set MSL_DIR to a random value and see that the output changes.
 */

int main(){

	cout << "MSL_DIR is set to: "<<SYSENV.getEnv("MSL_DIR")<<endl;

	cout << "MSL_CHARMM_PAR is set to: "<<SYSENV.getEnv("MSL_CHARMM_PAR")  <<endl;
	cout << "MSL_CHARMM_TOP is set to: "<<SYSENV.getEnv("MSL_CHARMM_TOP")  <<endl;
                  
	cout << "MSL_HBOND_PAR is set to: "<<SYSENV.getEnv("MSL_HBOND_PAR")   	<<endl;
	cout << "MSL_HBOND_CA_PAR is set to: "<<SYSENV.getEnv("MSL_HBOND_CA_PAR")   	<<endl;
                  
	cout << "MSL_PDB_TOP is set to: "<<SYSENV.getEnv("MSL_PDB_TOP")   	<<endl;
                  
	cout << "MSL_ROTLIB is set to: "<<SYSENV.getEnv("MSL_ROTLIB")     <<endl;

	cout << "MSL_EXAMPLE_FILE_DIR is set to: "<<SYSENV.getEnv("MSL_EXAMPLE_FILE_DIR")     <<endl;

	cout << "FOOBAR is set to: "<<SYSENV.getEnv("FOOBAR") <<endl;
}
