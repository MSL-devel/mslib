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
#include "MslTools.h"
#include "Transforms.h"

using namespace std;
using namespace MSL;
int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cerr << "USAGE:\nexample_molecular_alignment <path_of_exampleFiles_directory>" << endl;
		exit(0);
	}

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     How to align two molecules (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;


	string refFile = "example0005.pdb";
	refFile = (string)argv[1] + "/" + refFile;
	cout << "Create an AtomContainer and read the atoms from " << refFile << endl;


	string moveFile = "example0006.pdb";
	moveFile = (string)argv[1] + "/" + moveFile;
	cout << "Create an AtomContainer and read the atoms from " << moveFile << endl;


 

	// read the PDB into a new AtomContainer "refAtoms"
	AtomContainer refAtoms;
	if (!refAtoms.readPdb(refFile)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << refFile << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
	  cout << "Reference PDB has "<<refAtoms.size()<<" atoms."<<endl;
	}

	// read the PDB into a new AtomContainer "moveAtoms"
	AtomContainer moveAtoms;
	if (!moveAtoms.readPdb(moveFile)) {
		// error checking, the PDB could not be read
		cerr << endl;
		cerr << "File " << moveFile << " cannot be found, please speficy the path of the \"exampleFiles\" directory as an argument" << endl;
		exit(1);
	} else {
	  cout << "Move PDB has "<<moveAtoms.size()<<" atoms."<<endl;
	}

	// Get RMSD before alignment
	double preRMSD = refAtoms.getAtomPointers().rmsd(moveAtoms.getAtomPointers());

	// Create a Transforms object (includes functions to do molecular alignments)
	Transforms alignAtoms;

	// Do molecular alignment first atoms are moving, second atoms are the reference
	if (!alignAtoms.rmsdAlignment(moveAtoms.getAtomPointers(),refAtoms.getAtomPointers())){
	  cerr << "ERROR alignment of molecules using all the atoms of reference: "<<refFile<<" and move: "<<moveFile<<" files."<<endl;
	  exit(0);
	}

	// Get RMSD after alignment
	double postRMSD = refAtoms.getAtomPointers().rmsd(moveAtoms.getAtomPointers());

	// Query Transform object for the rotation and translations required for alignment
	Matrix rotMatrix           = alignAtoms.getLastRotationMatrix();
	CartesianPoint translation = alignAtoms.getLastTranslation();

	// Output RMSD, Rotation Matrix and Translation Vector
	cout << endl << endl;
	fprintf(stdout,"Before alignment RMSD: %8.3f\n",preRMSD);
	fprintf(stdout,"After  alignment RMSD: %8.3f\n",postRMSD);
	
	cout << "Rotation Matix: " <<endl<<rotMatrix << endl;
	cout << "Translation Vector: "<<endl<<translation << endl;
}


