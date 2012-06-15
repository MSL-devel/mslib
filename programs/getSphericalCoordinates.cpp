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

//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes
#include "getSphericalCoordinates.h"
#include "System.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Frame.h"



using namespace std;

using namespace MSL;



int main(int argc, char *argv[]){
    
	// Option Parser
        Options opt = setupOptions(argc,argv);

	// Read PDB
	System sys;
	sys.readPdb(opt.pdb);

	int centerResidueIndex = -1;
	for (uint i = 0; i < sys.positionSize();i++){
		if ((sys.getPosition(i).getResidueNumber() == opt.resnum) && (sys.getPosition(i).getChainId() == opt.chain)) {
			centerResidueIndex = i;
			break;
		}
	}
	if (centerResidueIndex == -1){
		cerr << "ERROR Couldn't find central residue"<<endl;
		cerr << opt.pdb << endl;
		cerr << opt.resnum << "\t" << opt.chain << endl;
		exit(3);
	}
	cout << "Compute Frame.h"<<endl;
 	Frame f;
	if (!f.computeFrameFromFunctionalGroup(sys.getPosition(centerResidueIndex).getCurrentIdentity())){
			cerr << "Problem creating frame from central residue"<<endl;
			cerr << opt.pdb << endl;
			cerr << opt.resnum << "\t" << opt.chain << endl;
			exit(324);
	}

	if (opt.printFrames){
		cout << "Write out basic frame"<<endl;
		ofstream fout;
		fout.open("basicFrame.py");
		fout << f.toString()<<endl;
		fout.close();
	}


	// Align frame and atoms of sys to origin.
	AtomPointerVector &av = sys.getAtomPointers();
	f.transformToGlobalBasis(av);


	Transforms t;
	for (uint i = 0; i < sys.positionSize();i++){
		if (i == centerResidueIndex) continue;
		Residue &r  = sys.getResidue(i);
		CartesianPoint cp;

		double angleBetweenFrames = MslTools::doubleMax;
		if (r.getResidueName() == "ASP" &&
			    r.atomExists("CG") &&
			    r.atomExists("OD1") &&
			    r.atomExists("OD2")){
				
				cp.setCoor( (r("CG").getX()+r("OD1").getX()+r("OD2").getX()) / 3, 
					    (r("CG").getY()+r("OD1").getY()+r("OD2").getY()) / 3, 
					    (r("CG").getZ()+r("OD1").getZ()+r("OD2").getZ()) / 3);


				if (i != centerResidueIndex) {

					f.transformFromGlobalBasis(r.getAtomPointers());

					Frame floatingFrame;
					floatingFrame.computeFrameFrom3Atoms(r("OD1"),r("CG"),r("OD2"));

					if (opt.printFrames){

						
						char name[80];
						sprintf(name,"aspFrame%1s%03d.py",r.getChainId().c_str(),r.getResidueNumber());

						ofstream fout;
						fout.open(name);
						fout << floatingFrame.toString()<<endl;
						fout.close();
					}



					//cout << floatingFrame.toString()<<endl;

					Matrix m = f.anglesBetweenFrame(floatingFrame);
					angleBetweenFrames = m[2][2];// z vs z

					f.transformToGlobalBasis(r.getAtomPointers());

				}
				
		}

		if (r.getResidueName() == "ASN" &&
			    r.atomExists("CG") &&
			    r.atomExists("OD1") &&
			    r.atomExists("ND2")){
				
				cp.setCoor( (r("CG").getX()+r("OD1").getX()+r("ND2").getX()) / 3, 
					    (r("CG").getY()+r("OD1").getY()+r("ND2").getY()) / 3, 
					    (r("CG").getZ()+r("OD1").getZ()+r("ND2").getZ()) / 3);



				if (i != centerResidueIndex) {
					f.transformFromGlobalBasis(r.getAtomPointers());
					Frame floatingFrame;
					floatingFrame.computeFrameFrom3Atoms(r("OD1"),r("CG"),r("ND2"));

					if (opt.printFrames){
						char name[80];
						sprintf(name,"asnFrame%1s%03d.py",r.getChainId().c_str(),r.getResidueNumber());

						ofstream fout;
						fout.open(name);
						fout << floatingFrame.toString()<<endl;
						fout.close();
					}

					Matrix m = f.anglesBetweenFrame(floatingFrame);
					angleBetweenFrames = m[2][2];// z vs z
					f.transformToGlobalBasis(r.getAtomPointers());
				}
				
		}


		if (r.getResidueName() == "GLU" &&
			    r.atomExists("CD") &&
			    r.atomExists("OE1") &&
			    r.atomExists("OE2")){
				
				cp.setCoor( (r("CD").getX()+r("OE1").getX()+r("OE2").getX()) / 3, 
					    (r("CD").getY()+r("OE1").getY()+r("OE2").getY()) / 3, 
					    (r("CD").getZ()+r("OE1").getZ()+r("OE2").getZ()) / 3);



				if (i != centerResidueIndex) {
					f.transformFromGlobalBasis(r.getAtomPointers());
					Frame floatingFrame;
					floatingFrame.computeFrameFrom3Atoms(r("OE1"),r("CD"),r("OE2"));

					if (opt.printFrames){
						char name[80];
						sprintf(name,"gluFrame%1s%03d.py",r.getChainId().c_str(),r.getResidueNumber());

						ofstream fout;
						fout.open(name);
						fout << floatingFrame.toString()<<endl;
						fout.close();
					}

					Matrix m = f.anglesBetweenFrame(floatingFrame);
					angleBetweenFrames = m[2][2];// z vs z
					f.transformToGlobalBasis(r.getAtomPointers());
				}
				
		}

		if (r.getResidueName() == "GLN" &&
			    r.atomExists("CD") &&
			    r.atomExists("OE1") &&
			    r.atomExists("NE2")){
				
				cp.setCoor( (r("CD").getX()+r("OE1").getX()+r("NE2").getX()) / 3, 
					    (r("CD").getY()+r("OE1").getY()+r("NE2").getY()) / 3, 
					    (r("CD").getZ()+r("OE1").getZ()+r("NE2").getZ()) / 3);


				if (i != centerResidueIndex) {
					f.transformFromGlobalBasis(r.getAtomPointers());
					Frame floatingFrame;
					floatingFrame.computeFrameFrom3Atoms(r("OE1"),r("CD"),r("NE2"));
					if (opt.printFrames){
						char name[80];
						sprintf(name,"glnFrame%1s%03d.py",r.getChainId().c_str(),r.getResidueNumber());

						ofstream fout;
						fout.open(name);
						fout << floatingFrame.toString()<<endl;
						fout.close();
					}

					Matrix m = f.anglesBetweenFrame(floatingFrame);
					angleBetweenFrames = m[2][2];// z vs z
					f.transformToGlobalBasis(r.getAtomPointers());
				}

				
		}

		
		SphericalPoint spRes = t.transform(cp);
		if (angleBetweenFrames == MslTools::doubleMax){
			angleBetweenFrames = 0;
		} else {
			if (angleBetweenFrames > 90){
				angleBetweenFrames = 180 - angleBetweenFrames;
			}

		}

 		fprintf(stdout, "RES %10s %1s %04d %3s %8.3f %8.3f %8.3f %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),r.getChainId().c_str(),r.getResidueNumber(),r.getResidueName().c_str(),spRes.getRadius(), spRes.getSigma(),spRes.getTheta(),angleBetweenFrames*M_PI/180);


		AtomPointerVector &ats = r.getAtomPointers();
		for (uint j = 0; j < ats.size();j++){
			SphericalPoint sp = t.transform(ats(j).getCoor());
			
			//fprintf(stdout, "ATM %10s %1s %04d %3s %4s %8.3f %8.3f %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),ats(j).getChainId().c_str(),ats(j).getResidueNumber(),r.getResidueName().c_str(),ats(j).getName().c_str(),sp.getRadius(), sp.getSigma(),sp.getTheta());
			Residue &cent = sys.getPosition(centerResidueIndex).getCurrentIdentity();
 			if (cent.getResidueName() == "LYS"){
	 			double lysDihedral = 0;
 				double lysAngle    = 0;
 				lysDihedral = cent("CD").getCoor().dihedral(cent("CE").getCoor(),cent("NZ").getCoor(),ats(j).getCoor());
 				lysAngle    = cent("CE").getCoor().angle(cent("NZ").getCoor(),ats(j).getCoor());
				fprintf(stdout, "ATM %10s %04d %1s %04d %3s %4s %8.3f %8.3f %8.3f %8.3f %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),opt.resnum,ats(j).getChainId().c_str(),ats(j).getResidueNumber(),r.getResidueName().c_str(),ats(j).getName().c_str(),sp.getRadius(), sp.getSigma(),sp.getTheta(),lysAngle,lysDihedral);
 			}
			else {
				// JEDONALD WAY..
				fprintf(stdout, "ATM %10s %04d %1s %04d %3s %4s %8.3f %8.3f %8.3f %8.3f\n",MslTools::getFileName(opt.pdb).c_str(),opt.resnum,ats(j).getChainId().c_str(),ats(j).getResidueNumber(),r.getResidueName().c_str(),ats(j).getName().c_str(),sp.getRadius(), sp.getSigma(),sp.getTheta(), angleBetweenFrames*M_PI/180);
			}
		}

		
	}
	
		
	cout << "Done."<<endl;
}

Options setupOptions(int theArgc, char * theArgv[]){

	// Create the options
	Options opt;
	
	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: getSphericalCoordinates conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.configFile = OP.getString("configfile");
	
	if (opt.configFile != "") {
		OP.readFile(opt.configFile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.resnum = OP.getInt("resnum");
	if (OP.fail()){
		cerr << "ERRROR 1111 no resnum\n";
		exit(1111);
	}

	opt.chain = OP.getString("chain");
	if (OP.fail()){
		cerr << "ERRROR 1111 no chain\n";
		exit(1111);
	}

	opt.negativeRes = OP.getBool("neg");
        if (OP.fail()){
                opt.negativeRes = false;
        }

        if (opt.negativeRes) {
                opt.resnum = -1*opt.resnum;
        }

	opt.printFrames = OP.getBool("printFrames");
	return opt;
}
