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
#include <cstdlib>

#include "AtomContainer.h"
#include "System.h"
#include "MslTools.h"
#include "CoiledCoils.h"
#include "Symmetry.h"


using namespace std;
using namespace MSL;

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 1) {
		cerr << "USAGE:\nexample_coiled_coil_and_symmetric_bundles" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to generate coiled-coils in MSL (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	// CoiledCoils object used to create a coiled helix (by a number of algorithms)
	CoiledCoils cc;

	// A super-helical radius
	double sr = 6.5;

	// An alpha-helical phase angle
	double aph = 0.0;

	// A super-helical pitch angle, and pitch distance
	double shpa = 190;
	double shPitch = (2*M_PI*sr)/tan(M_PI*shpa/180);

	// Hard code values of h (rise/residue) = 1.51, r1 (alpha-helical radius), and theta (alpha helical frequency)
	double risePerResidue = 1.51;

        // Use observed medians by Gevorg Grigoryan (Probing Deisgnability via a Generalized Model of Helical Bundle Geometry, JMB 2010)
	double alphaHelicalRadius = 2.26;
	double alphaHelicalFrequency = 102.8;

	// Number of residues in coil ( lets do 4 heptads = 28 residues )
	double numberOfResidues = 28;

	// Generate a coiled coil, using specified parameters
    	//cc.northCoiledCoils(sr, risePerResidue, shPitch, alphaHelicalRadius, numberOfResidues, alphaHelicalFrequency, aph);

	// dZ
	double dZ = 0.0;

	// Get the atoms from the CoiledCoils object back (this is a single coiled-coil helix)
	AtomPointerVector coil = cc.getCoiledCoil(sr, risePerResidue, shPitch, alphaHelicalRadius, alphaHelicalFrequency,dZ, aph,numberOfResidues);

	cout << "Writing /tmp/singleHelixCoil.pdb"<<endl;
	PDBWriter pout;
	pout.open("/tmp/singleHelixCoil.pdb");
	pout.write(coil);
	pout.close();

	// Create a symmtery object to generate a coiled bundle ( C4 symmetric )
	Symmetry sym;

	// Apply C4 to "coil"
	sym.applyCN(coil,4);
	cout << "Writing /tmp/C4HelixCoil.pdb"<<endl;
	pout.open("/tmp/C4HelixCoil.pdb");
	pout.write(sym.getAtomPointers());
	pout.close();

	
}


