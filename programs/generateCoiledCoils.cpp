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

#include <string>
#include <vector>
#include <ostream>
#include <fstream>
#include <cmath>

#include "generateCoiledCoils.h"
#include "CoiledCoils.h"
#include "Symmetry.h"
#include "OptionParser.h"
#include "PDBWriter.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;


int main(int argc, char *argv[]){
	// Option Parser
	Options opt = setupOptions(argc,argv);
	Transforms tr;

	// Super-helical Radius Loop
	for (double sr = opt.superHelicalRadius[0]; sr <= opt.superHelicalRadius[1]; sr += opt.superHelicalRadius[2]){

		// Alpha-helical Phase Angle Loop
		for (double aph = opt.alphaHelicalPhaseAngle[0]; aph < opt.alphaHelicalPhaseAngle[1];aph+=opt.alphaHelicalPhaseAngle[2]){

		  // Super-helical Pitch Angle loop added by David Slochower
		  for(double shpa = opt.superHelicalPitchAngle[0]; shpa < opt.superHelicalPitchAngle[1]; shpa+=opt.superHelicalPitchAngle[2]) {
		    double shPitch = (2*M_PI*sr)/tan(M_PI*shpa/180);

			// Generate a coiled helix
			CoiledCoils cc;
			cc.northCoiledCoils(sr, 1.5232, shPitch, 2.25, opt.numberOfResidues, 103.195, aph);

			AtomPointerVector coil = cc.getAtomPointers();

			// Apply symmetry operations to create a bundle
			if (opt.symmetry == "C2"){

				Symmetry sym;
				sym.applyCN(coil,2);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}

			if (opt.symmetry == "C3"){

				Symmetry sym;
				sym.applyCN(coil,3);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}
			
			if (opt.symmetry == "C4"){
				Symmetry sym;
				sym.applyCN(coil,4);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}

			if (opt.symmetry == "C5"){
				Symmetry sym;
				sym.applyCN(coil,5);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}

			if (opt.symmetry == "C6"){
				Symmetry sym;
				sym.applyCN(coil,6);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}


			if (opt.symmetry == "C7"){
				Symmetry sym;
				sym.applyCN(coil,7);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}
			
			if (opt.symmetry == "C8"){
				Symmetry sym;
				sym.applyCN(coil,8);

				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%02d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.numberOfResidues, sr, aph, shpa);
				
				cout << "Writing "<<filename<<"."<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}

			if (opt.symmetry == "D2"){
				// Z Rotate 
				for (double spa = opt.superHelicalPhaseAngle[0]; spa < opt.superHelicalPhaseAngle[1]; spa += opt.superHelicalPhaseAngle[2]){
					coil.clearSavedCoor();
					coil.saveCoor("preSPA");

					Matrix zRot = CartesianGeometry::instance()->getZRotationMatrix(spa);
					//coil.rotate(zRot);
					tr.rotate(coil, zRot);
						
					// Z Trans
					for (double ztrans = opt.d2zTranslation[0];ztrans < opt.d2zTranslation[1]; ztrans += opt.d2zTranslation[2]){
						coil.saveCoor("preZtrans");

						CartesianPoint z(0,0,ztrans);
						//coil.translate(z);
						tr.translate(coil, z);

						Symmetry sym;
						sym.applyD2(coil);
							
						// Write out bundle
						char filename[80];
						sprintf(filename, "%s_%05.2f_%05.2f_%05.2f_%05.2f.pdb", opt.name.c_str(),sr, aph, spa, ztrans);
				
						cout << "Writing "<<filename<<"."<<endl;
						PDBWriter pout;
						pout.open(filename);
						pout.write(sym.getAtomPointers());
						pout.close();

						coil.applySavedCoor("preZtrans");
					} // Ztrans

					coil.applySavedCoor("preSPA");
				} // SHPA
			} // SPA

			if (opt.symmetry == "D3"){
					
				// Z Rotate 
				for (double spa = opt.superHelicalPhaseAngle[0]; spa < opt.superHelicalPhaseAngle[1]; spa += opt.superHelicalPhaseAngle[2]){
					coil.clearSavedCoor();
					coil.saveCoor("preSPA");

					Matrix zRot = CartesianGeometry::instance()->getZRotationMatrix(spa);
					//coil.rotate(zRot);
					tr.rotate(coil, zRot);
						
					// Z Trans
					for (double ztrans = opt.d2zTranslation[0];ztrans < opt.d2zTranslation[1]; ztrans += opt.d2zTranslation[2]){
						coil.saveCoor("preZtrans");

						CartesianPoint z(0,0,ztrans);
						//coil.translate(z);
						tr.translate(coil, z);

						Symmetry sym;
						sym.applyDN(coil,3);
							
						// Write out bundle
						char filename[80];
						sprintf(filename, "%s_%05.2f_%05.2f_shp%05.2f_%05.2f_%05.2f.pdb", opt.name.c_str(),sr, aph, shpa, spa, ztrans);
				
						cout << "Writing "<<filename<<"."<<endl;
						PDBWriter pout;
						pout.open(filename);
						pout.write(sym.getAtomPointers());
						pout.close();

						coil.applySavedCoor("preZtrans");
					} // Ztrans
					coil.applySavedCoor("preSPA");
				} // SHPA
			} // SPA
                } // D2
           }
	}
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "generateCoiledBundles --symmetry SYM --superHelicalRadius LOW HIGH STEP --alphaHelicalPhaseAngle LOW HIGH STEP --superHelicalPitchAngle LOW HIGH STEP --numberOfResidues NUM [ --d2zTranslation LOW HIGH STEP --superHelicalPhaseAngle LOW HIGH STEP --name OUTFILE]\n";
		exit(0);
	}

	opt.symmetry = OP.getString("symmetry");
	if (OP.fail()){
		cerr << "ERROR 1111 symmetry not specified.\n";
		exit(1111);
	}
	opt.superHelicalRadius = OP.getDoubleVector("superHelicalRadius");
	if (OP.fail()){		   
		cerr << "ERROR 1111 superHelicalRadius not specified.\n";
		exit(1111);
	}

	opt.alphaHelicalPhaseAngle = OP.getDoubleVector("alphaHelicalPhaseAngle");
	if (OP.fail()){		   
		cerr << "ERROR 1111 alphaHelicalPhaseAngle not specified.\n";
		exit(1111);
	}

	opt.superHelicalPitchAngle = OP.getDoubleVector("superHelicalPitchAngle");
	if (OP.fail()){		   
		cerr << "ERROR 1111 superHelicalPitchAngle not specified.\n";
		exit(1111);
	}

	opt.numberOfResidues = OP.getInt("numberOfResidues");
	if (OP.fail()){		   
		cerr << "ERROR 1111 numberOfResidues not specified.\n";
		exit(1111);
	}

	opt.d2zTranslation         = OP.getDoubleVector("d2zTranslation");
	opt.superHelicalPhaseAngle = OP.getDoubleVector("superHelicalPhaseAngle");
	if (opt.d2zTranslation.size() == 0) { vector<double> tmp; tmp.push_back(0.0); tmp.push_back(1.0); tmp.push_back(100.0); opt.d2zTranslation = tmp; cerr << "Warning,  d2zTranslation not defined (required for D2, D3)" << endl; }
	if (opt.superHelicalPhaseAngle.size() == 0) { vector<double> tmp; tmp.push_back(0.0); tmp.push_back(1.0); tmp.push_back(100.0); opt.superHelicalPhaseAngle = tmp; cerr << "Warning,  superHelicalPhaseAngle not defined (required for D2, D3)" << endl; }

	opt.name = OP.getString("name");
	if (OP.fail()){
		opt.name = "CoiledBundle";
	}

	return opt;
}

