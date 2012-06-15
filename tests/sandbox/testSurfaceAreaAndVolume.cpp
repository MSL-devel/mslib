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


#include "SurfaceAreaAndVolume.h"
#include "CharmmParameterReader.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "System.h"
#include "testData.h"
#include <fstream>

using namespace MSL;
using namespace std;

#include "SysEnv.h"
static SysEnv SYSENV;

void run_1j4m_Test();

int main() {


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

    a.setCoor(0, 0, 0);
    a.setChainId("A");
    a.setResidueNumber(1);
    a.setResidueName("ALA");
    b.setCoor(0, 0, 2);
    b.setChainId("B");
    b.setResidueNumber(1);
    b.setResidueName("ALA");


    av.push_back(&a);
    av.push_back(&b);

    radii.push_back(2);
    radii.push_back(2);


    SurfaceAreaAndVolume sav;
    sav.computeSurfaceAreaAndVolume(av, radii);

    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n", "TEST1a", sav.getSurfaceArea());
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 56.54866776461628\n", "TEST1a", sav.getVolume());

    sav.computeSurfaceAreaAndVolume(av, radii);
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n", "TEST1a-atomwise1", sav.getSurfaceArea());
    double sumSA = 0.0;
    for (uint i = 0; i < av.size(); i++) {
        sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
    }
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n", "TEST1a-atomwise2", sumSA);


    // TRY READING USING addAtomsAndCharmmRadii + computeSurfaceAreaAndVolumeStereographicProjectIntegration()
    //                   get each atoms SASA and sum here to get the same number as above.

    System sys;

    // Add atoms to a system object
    sys.addAtoms(av);


    // Getting radii from parameter file 
    CharmmParameterReader par;
    par.reset();
    par.open(SYSENV.getEnv("MSL_CHARMM_PAR"));
    par.read();
    par.close();

    // atom type "C" has a radii of 2.0, so probeRadius will = 0
    sys.getAtomPointers()[0]->setType("C");
    sys.getAtomPointers()[1]->setType("C");
    sav.setProbeRadius(0.0);


    // Add atoms to SAV object
    sav.addAtomsAndCharmmRadii(sys.getAtomPointers(), par);
    sav.computeSurfaceAreaAndVolume();

    // Get Atomic SASA and sum.
    double aSA = (sav.getRadiiSurfaceAreaAndVolume(sys.getAtomPointers()[0]))[1];
    double bSA = (sav.getRadiiSurfaceAreaAndVolume(sys.getAtomPointers()[1]))[1];
    fprintf(stdout, " %-20s: Surface Area       %8.3f, expecting 75.39822368615504\n", "TEST1-atomwise3", aSA + bSA);


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
    a.setCoor(0, 0, 0);
    b.setCoor(0, 1.73205080756888, 1);
    av.push_back(&a);
    av.push_back(&b);


    sav.computeSurfaceAreaAndVolume(av, radii);
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n", "TEST1b-atomwise1", sav.getSurfaceArea());

    sumSA = 0.0;
    for (uint i = 0; i < av.size(); i++) {
        sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
    }
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 75.39822368615504\n", "TEST1b-atomwise2", sumSA);



    /*****************************************************************************
      TEST 2:

      3 Intersecting Spheres. 
          1. 0,0,0 radius 1
          2. 2,0,0 radius 2
          3. -0.75,2.90473750965556,0 radius 3

     V(Sa union Sb union Sc)  = 0.5736544730318
     A(Sa union Sb union Sc)  = 4.2139413434876
	  
     ******************************************************************************/


    a.setCoor(0, 0, 0);
    b.setCoor(2, 0, 0);
    c.setCoor(-0.75, 2.90473750965556, 0);


    // Get A-B
    av.clear();
    av.push_back(&a);
    av.push_back(&b);

    radii.clear();
    radii.push_back(1);
    radii.push_back(2);

    //cout << "GET SURFACE AREA AC\n";
    sav.computeSurfaceAreaAndVolume(av, radii);
    double saAB = sav.getSurfaceArea();
    double volAB = sav.getVolume();
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 54.9778714378214\n", "TEST2-AB", saAB);
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 35.9974158223830\n", "TEST2-AB", volAB);


    // Get B-C
    av.clear();
    av.push_back(&b);
    av.push_back(&c);

    radii.clear();
    radii.push_back(2);
    radii.push_back(3);

    //cout << "GET SURFACE AREA BC\n";
    sav.computeSurfaceAreaAndVolume(av, radii);
    double saBC = sav.getSurfaceArea();
    double volBC = sav.getVolume();
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.4402528821177\n", "TEST2-BC", saBC);
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 143.1388152791849\n", "TEST2-BC", volBC);

    // Get A-C
    av.clear();
    av.push_back(&a);
    av.push_back(&c);

    radii.clear();
    radii.push_back(1);
    radii.push_back(3);

    //cout << "GET SURFACE AREA AC\n";
    sav.computeSurfaceAreaAndVolume(av, radii);
    double saAC = sav.getSurfaceArea();
    double volAC = sav.getVolume();
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 117.2861257340189\n", "TEST2-AC", saAC);
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 115.4535300194249\n", "TEST2-AC", volAC);

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
    sav.computeSurfaceAreaAndVolume(av, radii);
    double saABC = sav.getSurfaceArea();
    double volABC = sav.getVolume();
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.9890027964171\n", "TEST2-ABC", saABC);
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 144.3669682217146\n", "TEST2-ABC", volABC);
    //sav.setDebug(false);



    sumSA = 0.0;
    for (uint i = 0; i < av.size(); i++) {
        sumSA += sav.getRadiiSurfaceAreaAndVolume(av[i])[1];
    }
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 148.9890027964171\n", "TEST2-ABC-atom2", sumSA);



    double Sa = 4 * M_PI * radii[0] * radii[0];
    double Sb = 4 * M_PI * radii[1] * radii[1];
    double Sc = 4 * M_PI * radii[2] * radii[2];

    double Va = (4.0 / 3.0) * M_PI * radii[0] * radii[0] * radii[0];
    double Vb = (4.0 / 3.0) * M_PI * radii[1] * radii[1] * radii[1];
    double Vc = (4.0 / 3.0) * M_PI * radii[2] * radii[2] * radii[2];
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 175.9291886010284\n", "TEST2-A+B+C", (Sa + Sb + Sc));
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 150.7964473723101\n", "TEST2-A+B+C", (Va + Vb + Vc));

    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 4.2139413434876\n", "TEST2-3spheres", saABC - (saAB + saAC + saBC)+(Sa + Sb + Sc));
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 0.5736544730318\n", "TEST2-3spheres", volABC - (volAB + volAC + volBC)+(Va + Vb + Vc));

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

    a.setCoor(0, 0, 0);
    av.push_back(&a);
    radii.push_back(2);

    b.setCoor(-5, 0, 0);
    av.push_back(&b);
    radii.push_back(6);

    c.setCoor(5, 0, 0);
    av.push_back(&c);
    radii.push_back(6);

    d.setCoor(0, 1, 0);
    av.push_back(&d);
    radii.push_back(4);

    e.setCoor(-1, 2, 3);
    av.push_back(&e);
    radii.push_back(2);

    f.setCoor(1, 1, 1);
    av.push_back(&f);
    radii.push_back(1);

    g.setCoor(0, 0, 1);
    av.push_back(&g);
    radii.push_back(2);

    h.setCoor(10, 0, 0);
    av.push_back(&h);
    radii.push_back(6);

    PDBWriter out;
    out.open("/tmp/foo.pdb");
    out.write(av);
    out.close();

    sav.computeSurfaceAreaAndVolume(av, radii);


    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 1011.872531375812\n", "TEST3-8spheres", sav.getSurfaceArea());
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 2329.934829835795\n", "TEST3-8spheres", sav.getVolume());

    // Test	4
    av.clear();
    radii.clear();
    a.setCoor(0.0, 0.0, 0.0);
    av.push_back(&a);
    radii.push_back(2.0);
    sav.computeSurfaceAreaAndVolume(av, radii);
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 50.26548245743664\n", "TEST4-1sphere", sav.getSurfaceArea());
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 33.510321638291093333333333333333\n", "TEST4-1sphere", sav.getVolume());

    // Test 5
    av.clear();
    radii.clear();
    double delta = 0.025;

    a.setCoor(0.0, 0.0, 1.0);
    av.push_back(&a);
    radii.push_back(2.0);

    b.setCoor(delta, 2.0 * delta, 0.0);
    av.push_back(&b);
    radii.push_back(1.0);

    c.setCoor(-delta, 2.0 * delta, 0.0);
    av.push_back(&c);
    radii.push_back(1.0);

    fprintf(stdout, " TEST5\n");
    for (int step = 0; step < 15; step++) {
        b.setCoor(delta, 2.0 * delta, 0.0);
        c.setCoor(-delta, 2.0 * delta, 0.0);
        sav.computeSurfaceAreaAndVolume(av, radii);
        double sa1 = sav.getRadiiSurfaceAreaAndVolume(&b)[1];
        double sa2 = sav.getRadiiSurfaceAreaAndVolume(&c)[1];
        fprintf(stdout, "   %3d %.4f %.5f %.5f\n", step, delta, sa1, sa2);
        delta += 0.1;
    }

    // Test 6
    run_1j4m_Test();
}

void run_1j4m_Test() {
    AtomPointerVector av;
    Atom *tmpAtom;
    vector<double> radii;
    SurfaceAreaAndVolume sav;
    Transforms tf;
        
    tmpAtom = new Atom();
    tmpAtom->setCoor(-6.266, -5.597, -3.441);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-6.000, -5.651, -2.273);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-4.960, -6.501, -2.322);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-4.981, -7.659, -2.274);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-7.215, -5.483, -1.756);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-7.187, -6.802, -0.965);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-8.507, -6.790, -1.117);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-8.191, -7.541, -2.401);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-8.791, -8.077, -2.962);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-9.235, -8.708, -2.983);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-9.172, -8.119, -3.707);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-4.559, -5.582, -1.699);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-3.472, -5.766, -1.187);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.527, -5.183, -1.741);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.708, -5.870, -1.247);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.705, -3.954, -2.050);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.470, -3.261, -1.997);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.781, -1.804, -1.935);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.664, -1.312, -2.627);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-0.769, -3.605, -3.273);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.614, -3.946, -2.865);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.361, -3.764, -4.126);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.645, -5.129, -4.554);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.024, -5.142, -5.895);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.070, -1.160, -1.049);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.247, 0.279, -0.865);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.094, 0.942, -0.883);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.029, 0.516, -0.242);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.987, 0.691, 0.403);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.948, -0.358, 1.458);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.564, -1.464, 1.260);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.441, -0.390, 2.719);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.381, -2.156, 2.351);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.704, -1.565, 3.284);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-0.774, 0.444, 3.503);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.306, -1.898, 4.557);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-0.351, 0.105, 4.759);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-0.609, -1.067, 5.364);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.141, 2.040, -1.576);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.419, 2.727, -1.799);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.193, 4.198, -1.578);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.312, 4.824, -2.165);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.904, 2.521, -3.223);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.895, 1.784, -3.881);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.059, 1.580, -3.307);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.001, 4.711, -0.679);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.873, 6.110, -0.344);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.216, 6.700, -0.009);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.934, 6.173, 0.841);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.908, 5.961, 0.785);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.252, 6.870, 1.789);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.207, 7.936, 2.245);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.739, 6.625, 2.784);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.465, 8.663, 3.223);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.010, 7.332, 3.777);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.929, 8.421, 4.162);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.163, 9.167, 5.122);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.439, 7.840, -0.641);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.701, 8.546, -0.475);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.847, 7.725, -1.049);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(7.010, 7.965, -0.750);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.940, 8.702, 1.023);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.575, 10.027, 1.197);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.647, 10.874, 0.360);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.059, 10.374, 2.280);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.456, 6.738, -1.854);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.429, 5.848, -2.487);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.669, 4.618, -1.622);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(7.445, 3.744, -1.986);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.963, 4.546, -0.510);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.127, 3.393, 0.380);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.021, 2.401, 0.090);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.844, 2.636, 0.366);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.935, 3.838, 1.801);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.726, 5.095, 2.012);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(6.213, 2.738, 2.807);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(8.101, 5.033, 2.615);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.466, 1.289, -0.446);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.549, 0.211, -0.813);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.208, -0.606, 0.415);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.065, -0.954, 1.217);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.192, -0.623, -1.916);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.559, 0.317, -2.922);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.209, -1.607, -2.555);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.941, -0.916, 0.518);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.434, -1.673, 1.649);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.428, -2.693, 1.207);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.290, -2.398, 0.862);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.683, -0.719, 2.495);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.587, 0.054, 3.352);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.770, 1.265, 3.464);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.417, -0.228, 4.238);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.575, 1.960, 4.259);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.230, 0.464, 5.027);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(4.407, 1.674, 5.134);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(5.227, 2.352, 5.905);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.861, -3.918, 1.240);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.926, -4.963, 0.865);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.287, -5.612, 2.057);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(0.670, -6.659, 2.593);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.766, -6.013, 0.407);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(1.789, -6.002, -1.025);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.623, -7.135, -0.979);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(3.642, -7.541, -0.601);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(2.319, -7.747, -1.776);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-0.742, -4.888, 2.443);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.532, -5.260, 3.612);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.826, -6.732, 3.572);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-1.292, -7.573, 4.285);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.796, -6.984, 2.777);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-3.037, -8.365, 2.620);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.406, -8.835, 1.380);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.314, -10.039, 1.192);
    av.push_back(tmpAtom);
    radii.push_back(1.400);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-4.384, -8.417, 2.186);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-4.293, -9.470, 3.224);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-5.410, -9.860, 2.501);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-5.840, -9.768, 3.865);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-6.891, -9.872, 4.145);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-7.500, -9.915, 3.789);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-7.232, -10.084, 4.852);
    av.push_back(tmpAtom);
    radii.push_back(1.550);

    tmpAtom = new Atom();
    tmpAtom->setCoor(-2.241, -7.938, 0.577);
    av.push_back(tmpAtom);
    radii.push_back(1.400);
    
    sav.computeSurfaceAreaAndVolume(av, radii);
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 1361.1362734952666\n", "TEST6a", sav.getSurfaceArea());
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 1042.3727637262789\n", "TEST6a", sav.getVolume());
    
    // Do some rotation and we should get the same result.
    tf.Xrotate(av, 23.3);
    tf.Yrotate(av, 196.2);
    tf.Zrotate(av, 270.6);
    sav.computeSurfaceAreaAndVolume(av, radii);
    fprintf(stdout, " %-20s: Surface Area       %8.3f , expected 1361.1362734952666\n", "TEST6b", sav.getSurfaceArea());
    fprintf(stdout, " %-20s: Volume             %8.3f , expected 1042.3727637262789\n", "TEST6b", sav.getVolume());
    
    av.deletePointers();
}
