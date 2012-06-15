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
			// Values used for previous work: cc.northCoiledCoils(sr, 1.5232, shPitch, 2.25, opt.numberOfResidues, 103.195, aph);
			// March 31, 2010: Jason Donald
			// Hard code values of h (rise/residue) = 1.51, r1 (alpha-helical radius), and theta (alpha helical frequency)
                        // based on median values observed by Gevorg Grigoryan
			//cc.northCoiledCoils(sr, 1.51, shPitch, 2.26, opt.numberOfResidues, 102.8, aph);
			
			AtomPointerVector coil = cc.getCoiledCoil(sr, 1.51, shPitch, 2.26, 102.8, aph, 0.0,opt.numberOfResidues); 

			// Apply symmetry operations to create a bundle
			int C_axis = atoi(opt.symmetry.substr(1,(opt.symmetry.length()-1)).c_str());
			if (opt.symmetry.substr(0,1) == "C"){
				Symmetry sym;
				sym.applyCN(coil,C_axis);
	
				// Write out bundle
				char filename[80];
				sprintf(filename, "%s_%s_%03d_%05.2f_%05.2f_shp%05.2f.pdb", opt.name.c_str(),opt.symmetry.c_str(),opt.numberOfResidues, sr, aph, shpa);
			
				cout << "Writing "<<filename<<endl;
				PDBWriter pout;
				pout.open(filename);
				pout.write(sym.getAtomPointers());
				pout.close();
			}	
			else if (opt.symmetry.substr(0,1) == "D"){
				// Z Rotate 
				for (double spa = opt.superHelicalPhaseAngle[0]; spa < opt.superHelicalPhaseAngle[1]; spa += opt.superHelicalPhaseAngle[2]){
					coil.clearSavedCoor();
					coil.saveCoor("preSPA");

					Matrix zRot = CartesianGeometry::getZRotationMatrix(spa);
					//coil.rotate(zRot);
					tr.rotate(coil, zRot);
					
					// Z Trans
					for (double ztrans = opt.d2zTranslation[0];ztrans < opt.d2zTranslation[1]; ztrans += opt.d2zTranslation[2]){
						coil.saveCoor("preZtrans");

						CartesianPoint z(0,0,ztrans);
						//coil.translate(z);
						tr.translate(coil, z);

						Symmetry sym;
						sym.applyDN(coil,C_axis);
								
						// Write out bundle
						char filename[80];
						sprintf(filename, "%s_%s_%03d_%05.2f_%05.2f_shp%05.2f_%05.2f_%05.2f.pdb", opt.name.c_str(),opt.symmetry.c_str(),opt.numberOfResidues,sr, aph, shpa, spa, ztrans);
			
						cout << "Writing "<<filename<<endl;
						PDBWriter pout;
						pout.open(filename);
						pout.write(sym.getAtomPointers());
						pout.close();

						coil.applySavedCoor("preZtrans");
					} // Ztrans
					coil.applySavedCoor("preSPA");
				} // SHA
			} 
                }
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
		cout << "Recommended parameters:" << endl;
		cout << "--symmetry: C2, C3, D2, D3, etc." << endl;
		cout << "--superHelicalRadius: radius from center in Angstroms" << endl;
		cout << "--alphaHelicalPhaseAngle: phase angle of the helices (0-360 degrees)" << endl;
		cout << "--superHelicalPitchAngle: often 5-20 degrees, but depends on the coil (average around 12)" << endl;
		cout << "--d2zTranslation: z offset between D_N symmetric coils, usually small" << endl;
		cout << "--superHelicalPhaseAngle: angle in degrees.  Value of 90/N for D_N typically works well, but larger and smaller values are possible" << endl;
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
	if (opt.superHelicalPhaseAngle.size() == 0) {
		int _N = atoi(opt.symmetry.substr(1,(opt.symmetry.length()-1)).c_str());
		double recommendedValue = 90./double(_N);
		vector<double> tmp; tmp.push_back(recommendedValue); tmp.push_back((recommendedValue + 1.)); tmp.push_back(100.0); opt.superHelicalPhaseAngle = tmp; cerr << "Warning,  superHelicalPhaseAngle not defined (required for D2, D3)" << endl;
	}

	opt.name = OP.getString("name");
	if (OP.fail()){
		opt.name = "CoiledBundle";
	}

	return opt;
}

