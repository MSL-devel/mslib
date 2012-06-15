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
#include "PDBWriter.h"
#include "PDBReader.h"
#include "testData.h"
#include "AtomPointerVector.h"
#include "CrystalLattice.h"
#include "Transforms.h"

using namespace MSL;
using namespace std;

void copyAtoms(const AtomPointerVector & _atoms, AtomPointerVector *newAts) ;
void invert(Matrix &a,Matrix &in);

int main(){

	// Write a pdb with proper REMARK 290
	writePdbFile();

	CrystalLattice cl("/tmp/xtalLattice.pdb");

	cl.generateCrystal();

	cl.writeCrystalUnits("/tmp/E",true,true,"A",false);
	
	exit(1);	

	PDBReader rin;
	rin.open("/tmp/xtalLattice.pdb");
	rin.read();
	rin.close();
	cout << "Done reading"<<endl;
	AtomPointerVector &ats = rin.getAtomPointers();

	PDBWriter wout;
	wout.open("/tmp/monomer.pdb");
	wout.write(ats);
	wout.close();


	Transforms tr;



	// Do something..
	vector<Matrix  *> &symMats          = rin.getSymmetryRotations();
	vector<CartesianPoint  *> &symTrans = rin.getSymmetryTranslations();
	map<string,double> &bounds          = rin.getBoundingCoordinates();
	Matrix &scaleMat                    = rin.getScaleRotation();
	CartesianPoint &scaleTrans          = rin.getScaleTranslation();
	vector<double>& unitCellParams      = rin.getUnitCellParameters();

	double a      = unitCellParams[0];
	double b      = unitCellParams[1];
	double c      = unitCellParams[2];
	double alpha  = unitCellParams[3] * M_PI / 180;
	double beta   = unitCellParams[4] * M_PI / 180;
	double gamma  = unitCellParams[5] * M_PI / 180;
	
	Matrix fractToOrtho (3,3,0.0);

	double cosAlpha  = ( ( cos(beta) * cos(gamma) - cos(alpha) ) / (sin(beta)*sin(gamma)));
	double sinAlpha  = sqrt(1 - cosAlpha*cosAlpha) ;
	fractToOrtho[0][0] = a;
	fractToOrtho[0][1] = b * cos(gamma);
	fractToOrtho[0][2] = c * cos(beta);
	fractToOrtho[1][0] = 0.0;
	fractToOrtho[1][1] = b *sin(gamma);
	fractToOrtho[1][2] = -c*sin(beta)* cosAlpha;
	fractToOrtho[2][0] = 0.0;
	fractToOrtho[2][1] = 0.0;
	fractToOrtho[2][2] = c*sin(beta) * sinAlpha;

        Matrix identityIHope = fractToOrtho*scaleMat;
        cout << identityIHope.toString() << endl;
	Matrix orthoToFract(3,3,0.0);
	orthoToFract[0][0] = 1.0 / a;
	orthoToFract[0][1] = -cos(gamma)/(sin(gamma) * a);
	orthoToFract[0][2] = -(cos(gamma)*sin(beta)*cosAlpha+cos(beta)*sin(gamma)) / ( sin(beta) * sinAlpha *sin(gamma)*a);
	orthoToFract[1][0] = 0.0;
	orthoToFract[1][1] = 1.0 / (sin(gamma) * b);
	orthoToFract[1][2] = cosAlpha / (sinAlpha * sin(gamma) * b);
	orthoToFract[2][0] = 0.0;
	orthoToFract[2][1] = 0.0;
	orthoToFract[2][2] =  1.0 / (sin(beta)*sinAlpha*c);

	CartesianPoint p1(1,0,0);
	CartesianPoint p2(0,1,0);
	CartesianPoint p3(0,0,1);

	CartesianPoint transP1 = CartesianGeometry::matrixTimesCartesianPoint(p1,fractToOrtho);
	CartesianPoint transP2 = CartesianGeometry::matrixTimesCartesianPoint(p2,fractToOrtho);
	CartesianPoint transP3 = CartesianGeometry::matrixTimesCartesianPoint(p3,fractToOrtho);

	cout << "P1: "<<transP1<<endl;
	cout << "P2: "<<transP2<<endl;
	cout << "P3: "<<transP3<<endl;

	ats.saveCoor("pre");
//	ats.translate(transP1);
	tr.translate(ats, transP1);

	wout.open("/tmp/p1p.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");

	transP1 *= -1;
//	ats.translate(transP1);
	tr.translate(ats, transP1);

	wout.open("/tmp/p1m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


//	ats.translate(transP2);
	tr.translate(ats, transP2);

	wout.open("/tmp/p2p.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


	transP2 *= -1;
//	ats.translate(transP2);
	tr.translate(ats, transP2);

	wout.open("/tmp/p2m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


//	ats.translate(transP3);
	tr.translate(ats, transP3);

	wout.open("/tmp/p3.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");

	transP3 *= -1;
//	ats.translate(transP3);
	tr.translate(ats, transP3);

	wout.open("/tmp/p3m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");
	


	Matrix scaleMatInv(3,3,0.0);
	invert(scaleMat, scaleMatInv);
	

	CartesianPoint s1(1,0,0);
	CartesianPoint s2(0,1,0);
	CartesianPoint s3(0,0,1);

	CartesianPoint transS1 = CartesianGeometry::matrixTransposeTimesCartesianPoint(s1,scaleMatInv);
	CartesianPoint transS2 = CartesianGeometry::matrixTransposeTimesCartesianPoint(s2,scaleMatInv);
	CartesianPoint transS3 = CartesianGeometry::matrixTransposeTimesCartesianPoint(s3,scaleMatInv);

	cout << "S1: "<<transS1<<endl;
	cout << "S2: "<<transS2<<endl;
	cout << "S3: "<<transS3<<endl;


	ats.saveCoor("pre");
//	ats.translate(transS1);
	tr.translate(ats, transS1);

	wout.open("/tmp/s1p.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");

	transS1 *= -1;
//	ats.translate(transS1);
	tr.translate(ats, transS1);

	wout.open("/tmp/s1m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


//	ats.translate(transS2);
	tr.translate(ats, transS2);

	wout.open("/tmp/s2p.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


	transS2 *= -1;
//	ats.translate(transS2);
	tr.translate(ats, transS2);

	wout.open("/tmp/s2m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");


//	ats.translate(transS3);
	tr.translate(ats, transS3);

	wout.open("/tmp/s3.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");

	transS3 *= -1;
//	ats.translate(transS3);
	tr.translate(ats, transS3);

	wout.open("/tmp/s3m.pdb");
	wout.write(ats);
	wout.close();

	ats.applySavedCoor("pre");




	/*

	// Loop 50 times

	ats.updateGeometricCenter();

	vector<AtomPointerVector *> symMates;	
	map<double,int> symMateDist;
	symMates.push_back(&ats);

	int symMatesStart = 0;
	int symMatesEnd   = symMates.size();

	int closeSymmetryMates = 0;
	int iterations = 0;
	cout << "Iterate now."<<endl;
	while (iterations++ < 2){


		// For each current symmetry mate, generate next level symmetry mates, test them,store them, print out.
		cout << "Create transforms "<<symMatesStart<<","<<symMatesEnd<<endl;
		for (uint i = symMatesStart; i < symMatesEnd;i++){
			

			for (uint j = 0; j < symMats.size();j++){

				AtomPointerVector *tmp = new AtomPointerVector();
				copyAtoms(*symMates[i],tmp);

				tmp->translate( (*symTrans[j] * -1) );
				tmp->rotate(*symMats[j]);

				tmp->updateGeometricCenter();
				
				double dist = ats.getGeometricCenter().distance(tmp->getGeometricCenter());
				cout << "Dist: "<<dist<<" "<<bounds["maxDelta"]<<endl;
				if (1){//(symMateDist[dist] != 1){ //&& dist != 0  && dist < 2*bounds["maxDelta"]){

					symMateDist[dist] = 1;
					symMates.push_back(tmp);


					char name[80];
					sprintf(name,"/tmp/lattice%06d-%06d-%06d-%06d.pdb", closeSymmetryMates++,iterations-1,i,j);
					fprintf(stdout,"Generated %s\n",name);
					wout.open(name);
					wout.write(*tmp);
					wout.close();
					
					
				} else {
					tmp->deletePointers();
				}

			}

			
			
		}
		symMatesStart = symMatesEnd;
		symMatesEnd   = symMates.size();
		
		
	}
	/*
	fprintf(stdout,"Number of symmetry mates for original monomer: %-4d\n",symMats.size());
	for (uint i = 0; i < symMats.size();i++){

		ats.saveCoor("pre");
		ats.translate( (*symTrans[i] * -1) );
		ats.rotate(*symMats[i]);

		for (uint j = 0; j < symMats.size();j++){

			ats.saveCoor("pre2");
			ats.translate( (*symTrans[j] * -1) );
			ats.rotate(*symMats[j]);

			char name[80];
			sprintf(name,"/tmp/lattice%04d-%04d.pdb",i,j);
			fprintf(stdout,"Generated %s\n",name);
			wout.open(name);
			wout.write(ats);
			wout.close();

			ats.applySavedCoor("pre2");
		}
		ats.applySavedCoor("pre");

	}
	*/
	
}
void copyAtoms(const AtomPointerVector & _atoms, AtomPointerVector *newAts) {
	for (uint i = 0; i < _atoms.size();i++){
		newAts->push_back(new Atom(*_atoms[i]));
	}
}

void invert(Matrix &a,Matrix &in){

    	double det=a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);//adjoin
    	in[0][0]=(a[1][1]*a[2][2]-a[2][1]*a[1][2])/det;
    	in[0][1]=-(a[1][0]*a[2][2]-a[1][2]*a[2][0])/det;
    	in[0][2]=(a[1][0]*a[2][1]-a[2][0]*a[1][1])/det;
    	in[1][0]=-(a[0][1]*a[2][2]-a[0][2]*a[2][1])/det;
    	in[1][1]=(a[0][0]*a[2][2]-a[0][2]*a[2][0])/det;
    	in[1][2]=-(a[0][0]*a[2][1]-a[2][0]*a[0][1])/det;
    	in[2][0]=(a[0][1]*a[1][2]-a[0][2]*a[1][1])/det;
    	in[2][1]=-(a[0][0]*a[1][2]-a[1][0]*a[0][2])/det;
    	in[2][2]=(a[0][0]*a[1][1]-a[1][0]*a[0][1])/det;

}
