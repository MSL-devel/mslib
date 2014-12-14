/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2014 The MSL Developer Group (see README.TXT)
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

#include "testData.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "SysEnv.h"

using namespace std;
using namespace MSL;

static SysEnv SYSENV;

int main(int argc, char* argv[]) {
	if(argc < 2) {
		cerr << "Usage: " << argv[0] << " <pdbFileName>" << endl;
		exit(0);
	}

	string pdbFile = string(argv[1]);
	System sys;
	CharmmSystemBuilder csb(sys, SYSENV.getEnv("MSL_CHARMM_TOP"), SYSENV.getEnv("MSL_CHARMM_PAR"));
	csb.setBuildNonBondedInteractions(false);
	if(!csb.buildSystemFromPDB(pdbFile)) {
		cerr << "Unable to build system from " << pdbFile << endl;
		exit(0);
	}

	SystemRotamerLoader sysRotLoad(sys,SYSENV.getEnv("MSL_ROTLIB"), "library/BEBL_11-2014.txt");

	vector<Position*> positions = sys.getPositions();
	for(int i = 0; i < positions.size(); i++) {
		string resName = positions[i]->getResidueName();
		double phi = positions[i]->getPhi();
		double psi = positions[i]->getPsi();
		cout << positions[i]->getPositionId() << "," << positions[i]->getResidueName();
		if(phi != MslTools::doubleMax) {
			cout << " Phi: " << phi;
		} else {
			cout << " Phi: *" ;
		}
		if(psi != MslTools::doubleMax) {
			cout << " Psi: " << psi;
		} else {
			cout << " Psi: *" ;
		}
		cout << endl;
		if(!sysRotLoad.loadRotamers(positions[i],resName, "SL60.00")) {
			cerr << "Unable to load rotamers for " << positions[i]->getPositionId() << "," << resName << endl;
		} else {
			cout << "Loaded "<< positions[i]->getTotalNumberOfRotamers() << " rotamers " << endl;
		}
	}	

	if(!sys.writePdb("tmp.pdb")) {
		cerr << "Unable to write tmp.pdb" << endl;
		exit(0);
	}

	
	
}
