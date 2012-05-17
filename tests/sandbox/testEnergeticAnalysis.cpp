/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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
#include "EnergeticAnalysis.h"
#include "CharmmSystemBuilder.h"
#include "PolymerSequence.h"
#include "testData.h"

using namespace MSL;
using namespace std;

#include "SysEnv.h"
static SysEnv SYSENV;

int main() {
	writePdbFile();

	System sys;
	sys.readPdb("/tmp/xtalLattice.pdb");

	PolymerSequence pseq(sys);

	System outSys;
	string topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
	CharmmSystemBuilder CSB(outSys,topfile,parfile);

	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	int numAssignedAtoms = outSys.assignCoordinates(sys.getAtomPointers(),false);
	fprintf(stdout,"Number of assigned atoms: %d",numAssignedAtoms);

	// Build the all atoms without coordinates (not in initial PDB)
	outSys.buildAllAtoms();

	outSys.writePdb("/tmp/preEA.pdb");
	EnergeticAnalysis ea;

	cout << "Analyze "<<outSys.getResidue(15).toString()<<endl;
	ea.analyzePosition(outSys, 15);
}
