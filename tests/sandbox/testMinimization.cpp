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

#include "SysEnv.h"
#include "System.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "CharmmBondInteraction.h"
#include "PyMolVisualization.h"
#include "EnergySet.h"
#include "testData.h"
#include "AtomContainer.h"

#include "GSLMinimizer.h"
using namespace MSL;
using namespace std;


// These test Minimization with/without gradient.
void twoAtomMin();
void threeAtomMin();
void fourAtomMin();

// These test the gradients of the energy for different size(complex) systems.
void triPepTest();
void twoAtomTest();
void threeAtomTest();
void fourAtomTest();


static SysEnv SYSENV;

int main(){
/*
	// Just compute the gradient
	cout << "\n\n==== TWO ATOM CASE ====\n\n";
	twoAtomTest();


	cout << "\n\n==== THREE ATOM CASE ====\n\n";
	threeAtomTest();

	cout << "\n\n==== FOUR ATOM CASE ====\n\n";
	fourAtomTest();

	cout << "\n\n==== TRI PEP CASE ====\n\n";
	triPepTest();


	// Now try minimization..
	cout << "\n\n==== TWO ATOM MIN ====\n\n";
	twoAtomMin();

	cout << "\n\n==== THREE ATOM MIN ====\n\n";
	threeAtomMin();
*/
	cout << "\n\n==== FOUR ATOM MIN ====\n\n";
	fourAtomMin();

}



void twoAtomMin(){
	
	//Atom N ("N" ,  -0.557,   1.474, -12.560);
	//Atom CA("CA",  -1.686,   1.490, -11.633);
	//Atom CA("CA",  -2.686,   1.490, -11.633);

	Atom N ("N" , 0,0,0);
	Atom CA("CA", 5,0,0);

	AtomContainer atContainer;
	atContainer.addAtom(N);
	atContainer.addAtom(CA);
	cout << "Added to container"<<endl;
	AtomPointerVector &ats = atContainer.getAtomPointers();
	PDBWriter pout;
	pout.open("twoAtomsA.pdb");
	pout.write(ats);
	pout.close();

	CharmmBondInteraction *pCBI = new CharmmBondInteraction(ats(0), ats(1), 1, 2);
	double e = pCBI->getEnergy();
	cout << "Bond energy: "<<e<<endl;

	EnergySet es;
	es.addInteraction(pCBI);
        es.setTermActive("CHARMM_BOND",true);


	ats.saveCoor("preMin");
	double energy = es.calcEnergy();	
	fprintf(stdout,"Prior to minimization energy: %8.3f\n",energy);

	GSLMinimizer min(&es,&ats);
	min.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	min.minimize();
	

	pout.open("twoAtomsA.min.noGradient.pdb");
	pout.write(ats);
	pout.close();

	
	ats.applySavedCoor("preMin");
	GSLMinimizer minG(&es,&ats);
	minG.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	minG.minimize();

	pout.open("twoAtomsA.min.Gradient.pdb");
	pout.write(ats);
	pout.close();

	
}

void twoAtomTest(){
	
	//Atom N ("N" ,  -0.557,   1.474, -12.560);
	//Atom CA("CA",  -1.686,   1.490, -11.633);
	//Atom CA("CA",  -2.686,   1.490, -11.633);

	Atom N ("N" , 0,0,0);
	Atom CA("CA", 5,0,0);

	AtomPointerVector ats;
	ats.push_back(&N);
	ats.push_back(&CA);

	PDBWriter pout;
	pout.open("twoAtomsA.pdb");
	pout.write(ats);
	pout.close();

	CharmmBondInteraction *pCBI = new CharmmBondInteraction(N, CA, 1, 2);

	EnergySet es;
	es.addInteraction(pCBI);
	
	double energy = es.calcEnergy();	

	vector<double> gradient(2*3,0.0);
	es.calcEnergyGradient(gradient);


	es.printSummary();

	
	fprintf(stdout, "ENERGY: %8.3f\n",energy);

	// PyMolViz...
	PyMolVisualization py;
	int index = 0;
	for (uint i = 0; i < gradient.size()-2;i+=3){
		
		CartesianPoint grad(gradient[i],gradient[i+1],gradient[i+2]);

		if (abs(grad.length()) < pow(10.0,-6.0)){
			fprintf(stdout, "%s %12.3f %12.3f %12.3f ****ZERO****\n", ats(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2]);
			index++;
			continue;
			
		}
		grad = grad.getUnit();

		fprintf(stdout, "%s %12.3f %12.3f %12.3f - %12.3f %12.3f %12.3f\n", ats(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2],grad[0],grad[1],grad[2]);

		char name[80];
		sprintf(name,"grad%04d",index+1);
		py.createArrow(ats(index).getCoor(),grad,(string)name);
		index++;
	}

	ofstream fout("gradTwoAtomA.py");
	fout << py;
	fout <<endl;
	fout.close();
	
	// TOO CLOSE , move CA closer
	CA.setCoor(-1.4,   1.490, -11.933);

	pout.open("twoAtomsB.pdb");
	pout.write(ats);
	pout.close();

	energy = es.calcEnergy();
	for (uint i = 0; i < gradient.size();i++){
		gradient[i] = 0.0;
	}

	es.calcEnergyGradient(gradient);
	index = 0;
	fprintf(stdout, "ENERGY: %8.3f\n",energy);
	for (uint i = 0; i < gradient.size()-2;i+=3){
		
		CartesianPoint grad(gradient[i],gradient[i+1],gradient[i+2]);

		if (abs(grad.length()) < pow(10.0,-6.0)){
			fprintf(stdout, "%s %12.3f %12.3f %12.3f ****ZERO****\n", ats(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2]);
			index++;
			continue;
			
		}
		grad = grad.getUnit();

		fprintf(stdout, "%s %12.3f %12.3f %12.3f - %12.3f %12.3f %12.3f\n", ats(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2],grad[0],grad[1],grad[2]);

		char name[80];
		sprintf(name,"grad%04d",index+1);
		py.createArrow(ats(index).getCoor(),grad,(string)name);
		index++;
	}
	
	fout.open("gradTwoAtomB.py");
	fout << py;
	fout <<endl;
	fout.close();
	
}

void threeAtomMin() {

	Atom N ("N" ,  -0.557,   1.474, -12.560);
	Atom CA("CA",  -1.686,   1.490, -11.633);
	Atom C("C" ,  -1.553,   0.374, -10.582);

	AtomPointerVector ats;
	ats.push_back(&N);
	ats.push_back(&CA);
	ats.push_back(&C);

	PDBWriter pout;
	pout.open("threeAtomsA.pdb");
	pout.write(ats);
	pout.close();

	CharmmBondInteraction *pCBI1 = new CharmmBondInteraction(N, CA, 300, 2);
	CharmmBondInteraction *pCBI2 = new CharmmBondInteraction(CA, C, 300, 2);
	CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(N, CA, C, 65, 109.5*(M_PI/180));

	EnergySet esBond;
	EnergySet esAngle;
	EnergySet esBoth;

	esBond.addInteraction(pCBI1);
	esBond.addInteraction(pCBI2);

	esAngle.addInteraction(pCAI);

	esBoth.addInteraction(pCBI1);
	esBoth.addInteraction(pCBI2);
	esBoth.addInteraction(pCAI);
	
	double energyBond = esBond.calcEnergy();	
	vector<double> gradientBond(3*3,0.0);
	esBond.calcEnergyGradient(gradientBond);

	double energyAngle = esAngle.calcEnergy();
	vector<double> gradientAngle(3*3,0.0);
	esAngle.calcEnergyGradient(gradientAngle);


	double energyBoth = esBoth.calcEnergy();	
	vector<double> gradientBoth(3*3,0.0);
	esBoth.calcEnergyGradient(gradientBoth);


	fprintf(stdout, "ENERGIES (BOND,ANGLE,BOTH): %8.3f,%8.3f,%8.3f\n",energyBond,energyAngle,energyBoth);

	// PyMolViz...
	PyMolVisualization py;
	int index = 0;
	for (uint i = 0; i < gradientBoth.size()-2;i+=3){
		
		CartesianPoint gradBond(gradientBond[i],gradientBond[i+1],gradientBond[i+2]);
		CartesianPoint gradAngle(gradientAngle[i],gradientAngle[i+1],gradientAngle[i+2]);
		CartesianPoint gradBoth(gradientBoth[i],gradientBoth[i+1],gradientBoth[i+2]);

		if (abs(gradBond.length()) > pow(10.0,-6.0)){
			gradBond = gradBond.getUnit();
			gradBond *= -1;

			char name[80];
			sprintf(name,"gradBond%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradBond,(string)name);

		}

		if (abs(gradAngle.length()) > pow(10.0,-6.0)){
			gradAngle = gradAngle.getUnit();
			gradAngle *= -1;

			char name[80];
			sprintf(name,"gradAngle%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradAngle,(string)name);

		}

		if (abs(gradBoth.length()) > pow(10.0,-6.0)){
			gradBoth = gradBoth.getUnit();
			gradBoth *= -1;

			char name[80];
			sprintf(name,"gradBoth%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradBoth,(string)name);

		}

		fprintf(stdout, "%s %12.3f %12.3f %12.3f : %12.3f %12.3f %12.3f : %12.3f %12.3f %12.3f\n", ats(index).toString().c_str(),gradientBond[i],gradientBond[i+1],gradientBond[i+2],gradientAngle[i],gradientAngle[i+1],gradientAngle[i+2],gradientBoth[i],gradientBoth[i+1],gradientBoth[i+2]);

		index++;
	}
	cout << "DONE"<<endl;
	ofstream fout("gradThreeAtomA.py");
	fout << py;
	fout <<endl;
	fout.close();



	ats.saveCoor("preMin");

	cout << "BOND ENERGY: 	"<<esBond.calcEnergy()<<endl;
	GSLMinimizer minBonds(&esBond,&ats);
	minBonds.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD1);
	minBonds.minimize();

	pout.open("threeAtomsA.minBonds.noGradient.pdb");
	pout.write(ats);
	pout.close();	


	ats.applySavedCoor("preMin");
	cout << "ANGLE ENERGY: 	"<<esAngle.calcEnergy()<<endl;
	GSLMinimizer minAngle(&esAngle,&ats);
	minAngle.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD1);
	minAngle.minimize();
	

	pout.open("threeAtomsA.minAngle.noGradient.pdb");
	pout.write(ats);
	pout.close();	

	ats.applySavedCoor("preMin");
	cout << "BOTH ENERGY: 	"<<esBoth.calcEnergy()<<endl;
	GSLMinimizer minBoth(&esBoth,&ats);
	minBoth.setMinimizeAlgorithm(GSLMinimizer::NELDERMEAD1);
	minBoth.minimize();
	

	pout.open("threeAtomsA.minBoth.noGradient.pdb");
	pout.write(ats);
	pout.close();	



	// With Gradient...
	ats.applySavedCoor("preMin");
	cout << "BOND ENERGY: 	"<<esBond.calcEnergy()<<endl;
	GSLMinimizer minBondsG(&esBond,&ats);
	minBondsG.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	minBondsG.minimize();

	pout.open("threeAtomsA.minBonds.Gradient.pdb");
	pout.write(ats);
	pout.close();	


	ats.applySavedCoor("preMin");
	cout << "ANGLE ENERGY: 	"<<esAngle.calcEnergy()<<endl;
	GSLMinimizer minAngleG(&esAngle,&ats);
	minAngleG.setMinimizeAlgorithm(GSLMinimizer::STEEPEST_DESCENT);
	minAngleG.minimize();
	

	pout.open("threeAtomsA.minAngle.Gradient.pdb");
	pout.write(ats);
	pout.close();	

	ats.applySavedCoor("preMin");
	cout << "BOTH ENERGY: 	"<<esBoth.calcEnergy()<<endl;
	GSLMinimizer minBothG(&esBoth,&ats);
	minBothG.setMinimizeAlgorithm(GSLMinimizer::BFGS);
	minBothG.minimize();
	

	pout.open("threeAtomsA.minBoth.Gradient.pdb");
	pout.write(ats);
	pout.close();	


	/* EnergySets delete the memory for all interactions they contain. But here the same interaction is part of multiple
	*  EnergySets, so let's clear them in all but one EnergySets.
	*/
	esBond.clearAllInteractions();
	esAngle.clearAllInteractions();

	
}

void threeAtomTest() {

	Atom N ("N" ,  -0.557,   1.474, -12.560);
	Atom CA("CA",  -1.686,   1.490, -11.633);
	Atom C("C" ,  -1.553,   0.374, -10.582);


	AtomPointerVector ats;
	ats.push_back(&N);
	ats.push_back(&CA);
	ats.push_back(&C);

	PDBWriter pout;
	pout.open("threeAtomsA.pdb");
	pout.write(ats);
	pout.close();

	CharmmBondInteraction *pCBI1 = new CharmmBondInteraction(N, CA, 1, 2);
	CharmmBondInteraction *pCBI2 = new CharmmBondInteraction(CA, C, 1, 2);
	CharmmAngleInteraction *pCAI = new CharmmAngleInteraction(N, CA, C, 65, 109.5*(M_PI/180));

	EnergySet esBond;
	EnergySet esAngle;
	EnergySet esBoth;

	esBond.addInteraction(pCBI1);
	esBond.addInteraction(pCBI2);

	esAngle.addInteraction(pCAI);

	esBoth.addInteraction(pCBI1);
	esBoth.addInteraction(pCBI2);
	esBoth.addInteraction(pCAI);
	
	double energyBond = esBond.calcEnergy();	
	vector<double> gradientBond(3*3,0.0);
	esBond.calcEnergyGradient(gradientBond);

	double energyAngle = esAngle.calcEnergy();
	vector<double> gradientAngle(3*3,0.0);
	esAngle.calcEnergyGradient(gradientAngle);


	double energyBoth = esBoth.calcEnergy();	
	vector<double> gradientBoth(3*3,0.0);
	esBoth.calcEnergyGradient(gradientBoth);


	fprintf(stdout, "ENERGIES (BOND,ANGLE,BOTH): %8.3f,%8.3f,%8.3f\n",energyBond,energyAngle,energyBoth);

	// PyMolViz...
	PyMolVisualization py;
	int index = 0;
	for (uint i = 0; i < gradientBoth.size()-2;i+=3){
		
		CartesianPoint gradBond(gradientBond[i],gradientBond[i+1],gradientBond[i+2]);
		CartesianPoint gradAngle(gradientAngle[i],gradientAngle[i+1],gradientAngle[i+2]);
		CartesianPoint gradBoth(gradientBoth[i],gradientBoth[i+1],gradientBoth[i+2]);

		if (abs(gradBond.length()) > pow(10.0,-6.0)){
			gradBond = gradBond.getUnit();

			char name[80];
			sprintf(name,"gradBond%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradBond,(string)name);

		}

		if (abs(gradAngle.length()) > pow(10.0,-6.0)){
			gradAngle = gradAngle.getUnit();

			char name[80];
			sprintf(name,"gradAngle%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradAngle,(string)name);

		}

		if (abs(gradBoth.length()) > pow(10.0,-6.0)){
			gradBoth = gradBoth.getUnit();

			char name[80];
			sprintf(name,"gradBoth%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradBoth,(string)name);

		}

		fprintf(stdout, "%s %12.3f %12.3f %12.3f : %12.3f %12.3f %12.3f : %12.3f %12.3f %12.3f\n", ats(index).toString().c_str(),gradientBond[i],gradientBond[i+1],gradientBond[i+2],gradientAngle[i],gradientAngle[i+1],gradientAngle[i+2],gradientBoth[i],gradientBoth[i+1],gradientBoth[i+2]);

		index++;
	}
	cout << "DONE"<<endl;
	ofstream fout("gradThreeAtomA.py");
	fout << py;
	fout <<endl;
	fout.close();

	/* EnergySets delete the memory for all interactions they contain. But here the same interaction is part of multiple
	*  EnergySets, so let's clear them in all but one EnergySets.
	*/
	esBond.clearAllInteractions();
	esAngle.clearAllInteractions();


}


void fourAtomTest(){
}
void fourAtomMin() {

	Atom N1("N1" ,  -0.557,   1.474, -12.560);
	Atom CA1("CA1",  -1.686,   1.490, -11.633);
	Atom C1("C1" ,  -1.553,   0.374, -10.582);
	Atom N2("N2" ,  -1.553,   0.374, -9.582);


	AtomPointerVector ats;
	ats.push_back(&N1);
	ats.push_back(&CA1);
	ats.push_back(&C1);
	ats.push_back(&N2);

	PDBWriter pout;
	pout.open("fourAtomsA.pdb");
	pout.write(ats);
	pout.close();

	CharmmBondInteraction *pCBI1    = new CharmmBondInteraction(N1, CA1, 320, 1.43);
	CharmmBondInteraction *pCBI2    = new CharmmBondInteraction(CA1, C1, 250, 1.49);
	CharmmBondInteraction *pCBI3    = new CharmmBondInteraction(C1, N2, 370, 1.345);
	CharmmAngleInteraction *pCAI1   = new CharmmAngleInteraction(N1, CA1, C1, 50, 107*(M_PI/180));
	CharmmAngleInteraction *pCAI2   = new CharmmAngleInteraction(N2, C1, CA1, 80, 116.5*(M_PI/180));
	CharmmDihedralInteraction *pCDI = new CharmmDihedralInteraction(N1, CA1, C1, N2,1000, 3, 0*(M_PI/180));

	EnergySet esDihedral;
	EnergySet esAll;

	esAll.addInteraction(pCBI1);
	esAll.addInteraction(pCBI2);
	esAll.addInteraction(pCBI3);
	esAll.addInteraction(pCAI1);
	esAll.addInteraction(pCAI2);
	esAll.addInteraction(pCDI);

	esDihedral.addInteraction(pCDI);
	
	double energyAll = esAll.calcEnergy();	
	vector<double> gradientAll(3*ats.size(),0.0);
	esAll.calcEnergyGradient(gradientAll);

	double energyDihedral = esDihedral.calcEnergy();
	vector<double> gradientDihedral(3*ats.size(),0.0);
	esDihedral.calcEnergyGradient(gradientDihedral);

	fprintf(stdout, "PREMINIMIZATION ENERGIES (ALL,DIHEDRAL): %8.3f,%8.3f\n",energyAll,energyDihedral);

	// PyMolViz...
	PyMolVisualization py;
	int index = 0;
	for (uint i = 0; i < gradientAll.size()-2;i+=3){
		
		CartesianPoint gradAll(gradientAll[i],gradientAll[i+1],gradientAll[i+2]);
		CartesianPoint gradDihedral(gradientDihedral[i],gradientDihedral[i+1],gradientDihedral[i+2]);


		if (abs(gradAll.length()) > pow(10.0,-6.0)){
			gradAll = gradAll.getUnit();
			gradAll *= -1;

			char name[80];
			sprintf(name,"gradAll%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradAll,(string)name);

		}

		if (abs(gradDihedral.length()) > pow(10.0,-6.0)){
			gradDihedral = gradDihedral.getUnit();
			gradDihedral *= -1;

			char name[80];
			sprintf(name,"gradDihedral%04d",index+1);
			py.createArrow(ats(index).getCoor(),gradDihedral,(string)name);

		}

		fprintf(stdout, "%s %12.3f %12.3f %12.3f : %12.3f %12.3f %12.3f\n", ats(index).toString().c_str(),gradientAll[i],gradientAll[i+1],gradientAll[i+2],gradientDihedral[i],gradientDihedral[i+1],gradientDihedral[i+2]);

		index++;
	}
	cout << "DONE"<<endl;
	ofstream fout("gradFourAtomA.py");
	fout << py;
	fout <<endl;
	fout.close();


	ats.saveCoor("preMin");

	GSLMinimizer minAll(&esAll,&ats);
	minAll.setMinimizeAlgorithm(GSLMinimizer::BFGS);
	minAll.minimize();
	


	pout.open("fourAtomsA.minAll.Gradient.pdb");
	pout.write(ats);
	pout.close();	

	ats.applySavedCoor("preMin");
	GSLMinimizer minDihedral(&esDihedral,&ats);
	minDihedral.setMinimizeAlgorithm(GSLMinimizer::BFGS);
	minDihedral.minimize();
	
	energyAll = esAll.calcEnergy();	
	energyDihedral = esDihedral.calcEnergy();
	fprintf(stdout, "POSTMINIMIZATION ENERGIES (ALL,DIHEDRAL): %8.3f,%8.3f\n",energyAll,energyDihedral);

	pout.open("fourAtomsA.minDihedral.Gradient.pdb");
	pout.write(ats);
	pout.close();	


	esDihedral.clearAllInteractions();
	
}




void triPepTest(){

	// This stream contains a triPeptide pdb from testData.h
	stringstream ss;
	ss << triPeptide;
	
	PDBReader pRead(ss);

	System initSys(pRead.getAtomPointers());

	PolymerSequence pseq(initSys);
	string topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	string parfile = SYSENV.getEnv("MSL_CHARMM_PAR");

	System sys;
	CharmmSystemBuilder CSB(sys,topfile,parfile);

	//CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).

	// Build the all atoms without coordinates (not in initial PDB)
	sys.buildAllAtoms();

	cout << "SET MIN" <<endl;
	for (uint i = 0; i < sys.getAtomPointers().size();i++){
		sys.getAtom(i).setMinimizationIndex(i+1);
	}
	
	EnergySet *es = sys.getEnergySet();
	// Compute Energy
	cout << "ENERGY CALC"<<endl;
	double energy = es->calcEnergy();

	
	// Compute Gradient
	cout << "ENERGY GRADIENT CALC"<<endl;
	vector<double> gradient(sys.getAtomPointers().size()*3,0.0);
	es->calcEnergyGradient(gradient);

	cout << "PYMOL"<<endl;

	fprintf(stdout, "ENERGY: %8.3f\n",energy);

	// PyMolViz...
	PyMolVisualization py;
	int index = 0;
	for (uint i = 0; i < gradient.size();i+=3){
		
		CartesianPoint grad(gradient[i],gradient[i+1],gradient[i+2]);

		if (abs(grad.length()) < pow(10.0,-6.0)){
			fprintf(stdout, "%s %12.3f %12.3f %12.3f ****ZERO****\n", sys.getAtom(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2]);
			index++;
			continue;
			
		}

		grad = grad.getUnit() * -1;

		fprintf(stdout, "%s %12.3f %12.3f %12.3f - %12.3f %12.3f %12.3f\n", sys.getAtom(index).toString().c_str(),gradient[i],gradient[i+1],gradient[i+2],grad[0],grad[1],grad[2]);

		char name[80];
		sprintf(name,"grad%04d",index+1);
		py.createArrow(sys.getAtom(index).getCoor(),grad,(string)name);
		index++;
	}
	cout << "DONE" << endl;

	ofstream fout("grad.py");
	fout << py;
	fout <<endl;

}
