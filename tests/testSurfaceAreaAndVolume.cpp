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


#include "SurfaceAreaAndVolume.h"
#include "CharmmParameterReader.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "System.h"
#include "testData.h"
#include <fstream>

using namespace MSL;
using namespace std;



int main(){
	

	/*****************************************************************************
	  TEST 1a:

	  2 Intersecting Spheres. One at 0,0,0 and one at 0,0,2 both a radius of 2.

	  V(r) = (9/4)*pi*r^3
	  A(r) = 6*pi*r*r

	  V(2) = 56.54866776461628
	  A(2) = 75.39822368615504
	  
	 ******************************************************************************/

	AtomPointerVector av;
	Atom a;
	Atom b;
	Atom c;
	vector<double> radii;


	a.setCoor(0,0,0); a.setChainId("A"); a.setResidueNumber(1);a.setResidueName("ALA");
	b.setCoor(0,0,2); b.setChainId("B"); b.setResidueNumber(1);b.setResidueName("ALA");


	av.push_back(&a);
	av.push_back(&b);

	radii.push_back(2);
	radii.push_back(2);
	

	SurfaceAreaAndVolume sav;
	//sav.setDebug(true);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);

	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1a",sav.getSurfaceArea());
	fprintf(stdout, " %-20s: Volume             %8.3f , expected 56.54866776461628\n","TEST1a",sav.getVolume());

	//sav.setDebug("true");
	sav.computeTest(av,radii);
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1a-atomwise1",sav.getSurfaceArea());
	double sumSA = 0.0;
	for (uint i = 0; i < av.size();i++){
		sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
	}
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1a-atomwise2",sumSA);


	// TRY READING USING addAtomsAndCharmmRadii + computeSurfaceAreaAndVolumeStereographicProjectIntegration()
	//                   get each atoms SASA and sum here to get the same number as above.

	System sys;

	// Add atoms to a system object
	sys.addAtoms(av);


	// Getting radii from parameter file 
	CharmmParameterReader par;
	par.reset();
	par.open("/library/charmmTopPar/par_all22_prot.inp");
	par.read();
	par.close();

	// atom type "C" has a radii of 2.0, so probeRadius will = 0
	sys.getAtomPointers()[0]->setType("C");
	sys.getAtomPointers()[1]->setType("C");
	sav.setProbeRadius(0.0); 


	// Add atoms to SAV object
	sav.addAtomsAndCharmmRadii(sys.getAtomPointers(),par);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration();

	// Get Atomic SASA and sum.
	double aSA = (sav.getRadiiSurfaceAreaAndVolume(sys.getAtomPointers()[0]))[1];
	double bSA = (sav.getRadiiSurfaceAreaAndVolume(sys.getAtomPointers()[1]))[1];
	fprintf(stdout," %-20s: Surface Area        %8.3f, expecting 75.39822368615504\n","TEST1-atomwise3",aSA+bSA);
	

	// OCCLUDING POINTS...
	//sav.computeSurfaceAreaOccludingPoints(av,radii);
	//fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1aNumeric",sav.getSurfaceArea());
	

	/*****************************************************************************
	  TEST 1b:

	  2 Intersecting Spheres. One at 0,0,0 and one at 0,1.73205080756888,1 both a radius of 2.

	  V(r) = (9/4)*pi*r^3
	  A(r) = 6*pi*r*r

	  V(2) = 56.54866776461628
	  A(2) = 75.39822368615504
	  
	 ******************************************************************************/
	av.clear();
	a.setCoor(0,0,0);
	b.setCoor(0,1.73205080756888,1);
	av.push_back(&a);
	av.push_back(&b);
	//sav.setDebug(true);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);

	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1b",sav.getSurfaceArea());
	
	sav.computeTest(av,radii);
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1b-atomwise1",sav.getSurfaceArea());
	sumSA = 0.0;
	for (uint i = 0; i < av.size();i++){
		sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
	}
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n","TEST1b-atomwise2",sumSA);



	/*****************************************************************************
	  TEST 2:

	  3 Intersecting Spheres. 
              1. 0,0,0 radius 1
              2. 2,0,0 radius 2
              3. -0.75,2.90473750965556,0 radius 3

         V(Sa union Sb union Sc)  = 0.5736544730318
         A(Sa union Sb union Sc)  = 4.2139413434876
	  
	 ******************************************************************************/


	a.setCoor(0,0,0);
	b.setCoor(2,0,0);
	c.setCoor(-0.75,2.90473750965556,0);


	// Get A-B
	av.clear();
	av.push_back(&a);
	av.push_back(&b);

	radii.clear();
	radii.push_back(1);
	radii.push_back(2);

	//cout << "GET SURFACE AREA AC\n";
	//sav.setDebug(true);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);
	double saAB = sav.getSurfaceArea();
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 54.9778714378214\n","TEST2-AB",saAB);


	// Get B-C
	av.clear();
	av.push_back(&b);
	av.push_back(&c);

	radii.clear();
	radii.push_back(2);
	radii.push_back(3);

	//cout << "GET SURFACE AREA BC\n";
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);
	double saBC = sav.getSurfaceArea();
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.4402528821177\n","TEST2-BC",saBC);

	// Get A-C
	av.clear();
	av.push_back(&a);
	av.push_back(&c);

	radii.clear();
	radii.push_back(1);
	radii.push_back(3);

	//cout << "GET SURFACE AREA AC\n";
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);
	double saAC = sav.getSurfaceArea();
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 117.2861257340189\n","TEST2-AC",saAC);

	// Get A-B-C
	av.clear();
	av.push_back(&a);
	av.push_back(&b);
	av.push_back(&c);

	radii.clear();
	radii.push_back(1);
	radii.push_back(2);
	radii.push_back(3);


	//cout << "GET SURFACE AREA ABC\n";
	//sav.setDebug(true);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);
	double saABC = sav.getSurfaceArea();
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.9890027964171\n","TEST2-ABC",saABC);
	sav.computeTest(av,radii);
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.9890027964171\n","TEST2-ABC-atom1",sav.getSurfaceArea());
	sumSA = 0.0;
	for (uint i = 0; i < av.size();i++){
		sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
	}
	fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.9890027964171\n","TEST2-ABC-atom2",sumSA);



	double Sa = 4*M_PI*radii[0]*radii[0];
	double Sb = 4*M_PI*radii[1]*radii[1];
	double Sc = 4*M_PI*radii[2]*radii[2];
	fprintf(stdout," %-20s: Surface Area       %8.3f , expected 175.9291886010284\n","TEST2-A+B+C",(Sa+Sb+Sc));

	fprintf(stdout," %-20s: Surface Area       %8.3f , expected 4.2139413434876\n","TEST2-3spheres",saABC-(saAB+saAC+saBC)+(Sa+Sb+Sc));


	/*****************************************************************************
	  TEST 3:
              8 Spheres

              Volume: 2329.934829835795 
              SA    : 1011.872531375812
	  


           NOTES: IS Atom C big and encapsulated atom f???
	 ******************************************************************************/

	av.clear();
	radii.clear();

	// a,b,c
	Atom d;
	Atom e;
	Atom f;
	Atom g; 
	Atom h;

	a.setCoor(0,0,0); av.push_back(&a);
	radii.push_back(2);

	b.setCoor(-5,0,0); av.push_back(&b);
	radii.push_back(6);

 	c.setCoor(5,0,0);av.push_back(&c);
	radii.push_back(6);

	d.setCoor(0,1,0);av.push_back(&d);
	radii.push_back(4);

	e.setCoor(-1,2,3);av.push_back(&e);
	radii.push_back(2);

	f.setCoor(1,1,1);av.push_back(&f);
	radii.push_back(1);

	g.setCoor(0,0,1);av.push_back(&g);
	radii.push_back(2);

	h.setCoor(10,0,0);av.push_back(&h);
	radii.push_back(6);

	PDBWriter out;
	out.open("/tmp/foo.pdb");
	out.write(av);
	out.close();

	//xosav.setDebug(true);
	sav.computeSurfaceAreaAndVolumeStereographicProjectIntegration(av,radii);



	fprintf(stdout," %-20s: Surface Area       %8.3f , expected 1011.872531375812\n","TEST3-8spheres",sav.getSurfaceArea());


	/*
	string abc="abcdefgh";
	for (uint i = 0; i < av.size();i++){
		cout << av(i).getCoor()<<endl;
	}

	*/
	





	writePdbFile();
	System sys2;
	sys2.readPdb("/tmp/triPep.pdb");
	//sys2.readPdb("/home/dwkulp/software/msl/sele.pdb");
	System sys3;
	PolymerSequence seq(sys2);
	string topFile = "/library/charmmTopPar/top_all22_prot.inp";
	string parFile = "/library/charmmTopPar/par_all22_prot.inp";
	CharmmSystemBuilder CSB(topFile,parFile);
	CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.
	CSB.buildSystem(sys3,seq);
	int numAssignedAtoms = sys3.assignCoordinates(sys2.getAtomPointers(),false);
	sys3.buildAllAtoms(); 

	radii.clear();
	for (uint i = 0; i < sys3.getAtomPointers().size();i++){
		vector<double> vPar;
		if (!par.vdwParam(vPar, sys3.getAtom(i).getType())) {
			cerr << "ERROR vdw parameters not found for atom " << sys3.getAtom(i) << " (" << sys3.getAtom(i).getType() << ")" << endl;
		} else {
			
			cout << sys3.getAtom(i) << " VDW: "<<vPar[1]<<endl;
			radii.push_back(vPar[1]/2+1.4);
		}
	}

	//sav.setDebug(true);
	sav.computeTest(sys2.getAtomPointers(),radii);

	fprintf(stdout, "SA of testPdb.pdb: %8.3f\n",sav.getSurfaceArea());



}
