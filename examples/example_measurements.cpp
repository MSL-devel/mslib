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
#include "PhiPsiStatistics.h"
#include "MslTools.h"

using namespace std;
using namespace MSL;

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_measurements <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to make measurements in MSL (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string refFile = "example0005.pdb";
	refFile = (string)argv[1] + "/" + refFile;

	AtomContainer atoms;
	if (!atoms.readPdb(refFile)) {
		cerr << endl;
		cerr << "File " << refFile << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	}


	// Distance and Distance Squared bewteen two atoms atoms[0] and atoms[size()-1]   (first and last atom)
	double distance        = atoms[0].distance(atoms[atoms.size()-1]);
	double distanceSquared = atoms[0].distance2(atoms[atoms.size()-1]);

	// Angle between three atoms ( atoms[0] - atoms[1] - atoms[2] angle)
	double angleDegrees    = atoms[0].angle(atoms[1],atoms[2]);
	double angleRadians    = atoms[0].angleRadians(atoms[1],atoms[2]);
	

	// Dihedral between four atoms ( atoms[0] - atoms[1] - atoms[2] - atoms[3] )
	double dihedralDegrees    = atoms[0].dihedral(atoms[1],atoms[2],atoms[3]);
	double dihedralRadians    = atoms[0].dihedralRadians(atoms[1],atoms[2],atoms[3]);
	
	
	// Print out measurements
	fprintf(stdout, "%-15s = %8.3f\n","Distance",distance);
	fprintf(stdout, "%-15s = %8.3f\n","DistanceSquared",distanceSquared);
	fprintf(stdout, "%-15s = %8.3f\n","AngleDegrees",angleDegrees);
	fprintf(stdout, "%-15s = %8.3f\n","AngleRadians",angleRadians);
	fprintf(stdout, "%-15s = %8.3f\n","DihedralDegrees",dihedralDegrees);
	fprintf(stdout, "%-15s = %8.3f\n","DihedralRadians",dihedralRadians);


	// Measurement Namespace : CartesianGeometry


	// The above Atom-based functions are wrappers to the utility functions in the CartesianGeometry namespace.

	CartesianPoint pt1(0.0 , 0.0 , 0.0);
	CartesianPoint pt2(1.0 , 0.0 , 0.0);
	CartesianPoint pt3(1.0 , 1.0 , 0.0);
	CartesianPoint pt4(1.0 , 1.0 , 1.0);

	distance        = CartesianGeometry::distance(pt1,pt2);
	angleDegrees    = CartesianGeometry::angle(pt1,pt2,pt3);
	dihedralDegrees = CartesianGeometry::dihedral(pt1,pt2,pt3,pt4);

	cout << endl;

	fprintf(stdout, "%-15s = %8.3f\n","Point Distance",distance);
	fprintf(stdout, "%-15s = %8.3f\n","Point Angle",angleDegrees);
	fprintf(stdout, "%-15s = %8.3f\n","Point Dihedral",dihedralDegrees);


	// Get specific dihedrals... Phi/Psi of alpha amino acid proteins

	// First create a System object
	System sys;
	sys.addAtoms(atoms.getAtomPointers());

	cout << endl << endl << "Phi-Psi of residue "<<sys.getResidue(1).toString()<<endl;

	double phi = PhiPsiStatistics::getPhi(sys.getResidue(0), sys.getResidue(1));
	double psi = PhiPsiStatistics::getPsi(sys.getResidue(1), sys.getResidue(2));

	fprintf(stdout,"\tPHI: %8.3f\n\tPSI: %8.3f\n",phi,psi);
	
}


