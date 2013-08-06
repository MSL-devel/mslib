/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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
#include <bitset>
#include <math.h>
#include <algorithm>
#include <functional>
#include "RandomNumberGenerator.h"
#include "CartesianGeometry.h"

using namespace MSL;
using namespace std;

SurfaceAreaAndVolume::SurfaceAreaAndVolume() {
    surfaceArea = 0.0;
    volume = 0.0;
    probeRadius = 1.4;
    debug = false;
}

SurfaceAreaAndVolume::~SurfaceAreaAndVolume() {
    atoms.clear();
}

void SurfaceAreaAndVolume::computeSurfaceAreaOccludingPoints(AtomPointerVector &_atoms, vector<double> &_radii) {

    for (uint i = 0; i < _atoms.size(); i++) {


        // Create a set of points around sphere...
        vector<pair<double, CartesianPoint> > hiddenPoints; // inner vector holds an angle and a direction vector


        for (uint j = 0; j < _atoms.size(); j++) {
            if (i == j) continue;
            double dist = _atoms(i).distance(_atoms(j));

            if (dist < _radii[i] + _radii[j]) {

                CartesianPoint direction = (_atoms(j).getCoor() - _atoms(i).getCoor()).getUnit();
                double alpha = asin(_radii[j] / dist);

                fprintf(stdout, " DIRECTION %s ALPHA %8.3f\n", direction.toString().c_str(), alpha);
                hiddenPoints.push_back(pair<double, CartesianPoint > (alpha, direction));
            }
        }

        double exposedFraction = getExposedFraction(hiddenPoints, _radii[i]);

        double atomicSurface = 4 * M_PI * _radii[i] * _radii[i] * exposedFraction;

        surfaceArea += atomicSurface;
    }


}

double SurfaceAreaAndVolume::getExposedFraction(vector<pair<double, CartesianPoint> > &_hiddenPoints, double _radii) {
    /*
    double angularGridResolution = 10 * M_PI / 180;

    int totalNumberOfSurfacePoints = int(2*M_PI/angularGridResolution)*int(2*M_PI/angularGridResolution);
    bitset<totalNumberOfSurfacePoints> bs(true);
     */
    /*
      Spherical coordinate transformation

          x = r*cos(angle1)*cos(angle2);
          y = r*cos(angle1)*sin(angle2);
          z = r*sin(angle1);s

          AND

          r = sqrt(x*x + y*y + z*z);
          theta1 = atan(y/x)
          theta2 = acos(z/r)
	  
     */
    /*
    for (uint i = 0; i < _hiddenPoints.size();i++){
		
            // Get Theta1,Theta2,R for this hidden point (R=_radii)
            double theta1 = atan(_hiddenPoints[i].second.getY() / _hiddenPoints[i].second.getX());
            double theta2 = acos(_hiddenPoints[i].second.getZ() / _radii);

            // Convert Theta1-alpha and Theta1+alpha to indices in bs
            int startBitIndex1 = getBitIndex(theta1-_hiddenPoints[i].first,angularGridResolution);
            int endBitIndex1   = getBitIndex(theta1+_hiddenPoints[i].first,angularGridResolution);

            // Convert Theta2-alpha and Theta2+alpha to indices in bs
            int startBitIndex2 = getBitIndex(theta2-_hiddenPoints[i].first,angularGridResolution) + int(totalNumberOfSurfacePoints/2);
            int endBitIndex2   = getBitIndex(theta2+_hiddenPoints[i].first,angularGridResolution) + int(totalNumberOfSurfacePoints/2);

            //  Loop over changing bits to 0
            for (uint j = startBitIndex1; j <= endBitIndex1;j++){

                    for (uint k = startBinIndex2; k <= endBitIndex2;k++){
                            int bitIndex = j * numPoints + k;
                            bs.set(bitIndex,false);
                    }
            }



            //  Loop over changing bits to 0
            for (uint j = startBitIndex; j <= endBitIndex;j++){
                    bs.set(j,false);
            }
    }

    // DEBUG
    // loop over each angle set
    //   Get x,y,z for each bit
    //   Create Spheres using PyMolViz
    //   Color accoring to bit value
    // Write out PyMOLVizObj.
    if (debug){

            for (double angle1 = 0; angle1 < 2*M_PI; angle1+=angularGridResolution){
                    for (double angle2 = 0; angle2 < 2*M_PI; angle2+=angularGridResolution){
				
				
                    }
            }

			

    }


	

    return bs.count() / bs.size();

     */
    return 0.0;
}

/*
double SurfaceAreaAndVolume::getExposedFraction(vector<pair<double,CartesianPoint> > &_hiddenPoints,double _radii){

        double angularGridResolution = 10 * M_PI / 180;

        double totalCount    = 0.0;
        double occludedCount = 0.0;


        for (double angle1 = 0; angle1 < 2*M_PI; angle1+=angularGridResolution){
                for (double angle2 = 0; angle2 < 2*M_PI; angle2+=angularGridResolution){
                        totalCount++;


                        CartesianPoint pointOnSurface;

                        // x = r*cos(angle1)*cos(angle2);
                        // y = r*cos(angle1)*sin(angle2);
                        // z = r*sin(angle1);

                        pointOnSurface.setCoor(_radii*cos(angle1)*cos(angle2),_radii*cos(angle1)*sin(angle2),_radii*sin(angle1));
                        fprintf(stdout,"Point On Surface: %s -----> [ %8.3f %8.3f ]\n",pointOnSurface.toString().c_str(),angle1,angle2);

                        // Test this point agains all hiddenPoints..
                        bool hidden = false;
                        for (uint i = 0; i < _hiddenPoints.size();i++){


                                // beta = angle between this point and direction.
                                double beta = pointOnSurface.angle(_hiddenPoints[i].second) * M_PI/180;
                                fprintf(stdout, "ALPHA,BETA: %8.3f,%8.3f\n",_hiddenPoints[i].first,beta);
                                if (beta < _hiddenPoints[i].first){
                                        hidden = true;
                                        break;
                                }
				
                        }

                        if (hidden){
                                occludedCount++;
                        }
                }
        }

        fprintf(stdout," Number Points Total: %8.3f Occluded: %8.3f\n",totalCount,occludedCount);
        return (1 - occludedCount / totalCount);
}
 */
void SurfaceAreaAndVolume::computeSurfaceAreaAndVolume(AtomPointerVector &_atoms, vector<double> &_radii) {
    /*
      AVRO algorithm : 
      Jan Busa, Jozef Dzurina, Edik Hayryan, Shura Hayryan, Chin-Kun Hu, Jan Plavka, Imrich Pokorny, Jaroslav Skrivanek, Ming-Chya Wu, ARVO: A Fortran package for computing the solvent accessible surface area and the excluded volume of overlapping spheres via analytic equations, Computer Physics Communications, Volume 165, Issue 1, 1 January 2005, Pages 59-96, ISSN 0010-4655, DOI: 10.1016/j.cpc.2004.08.002.
     */
    _atoms.saveCoor("preSurfaceArea");
    bool rotatedMolecule;
    rotatedMolecule = false;

    /*
      Paper talks about whole molecule rotation such that the NorthPoles of each sphere aren't inside any other sphere.
     */



    if (debug) {
        cout << endl << "-->  DISCOVER COLLIDING SPHERES  <--" << endl << endl;
    }

    /****************************************************
         STEP 1: Compute neighbors for each atom

     *****************************************************/
    atomicRadiiSurfaceAreaAndVolume.clear();
    insideSphere.assign(_atoms.size(), false); // store spheres that are completely engulfed in other spheres..
    collidingSpheres.clear();
    collidingSpheres.resize(_atoms.size());
    for (uint i = 0; i < _atoms.size(); i++) {


        atomicRadiiSurfaceAreaAndVolume[_atoms[i]].push_back(_radii[i]);
        atomicRadiiSurfaceAreaAndVolume[_atoms[i]].push_back(0.0);
        atomicRadiiSurfaceAreaAndVolume[_atoms[i]].push_back(0.0);

        // Skip "ENGULFED" spheres
        if (insideSphere[i]) continue;


        // Check for neighbors...
        for (uint j = i + 1; j < _atoms.size(); j++) {
            double dist = _atoms(i).distance(_atoms(j));

            if (dist < _radii[i] + _radii[j]) {


                if (debug) {
                    fprintf(stdout, "Atom %4d is colliding with Atom %4d (and visa-versa)\n", i, j);
                }

                // STORE WHICH SPHERES ARE COLLIDING
                collidingSpheres[i].push_back(j);
                collidingSpheres[j].push_back(i);


                // ENGULFMENT CHECK on both i and j
                if (dist + _radii[i] < _radii[j]) {
                    insideSphere[i] = true;
                }
                if (dist + _radii[j] < _radii[i]) {
                    insideSphere[j] = true;
                }
            }

        }


    }

    filterEngulfedAtoms(_atoms);

    rotatedMolecule = rotateMolecule(_atoms);


    createStereographicProjectedCircles(_atoms);
    getIntersectingAngles(_atoms);
    getArcs(_atoms);
    integrateArcs(_atoms);

    if (rotatedMolecule) {
        _atoms.applySavedCoor("preSurfaceArea");
    }

}

void SurfaceAreaAndVolume::computeSurfaceAreaAndVolume() {
    /*
      AVRO algorithm
     */
    atoms.saveCoor("preSurfaceArea");
    bool rotatedMolecule;
    rotatedMolecule = false;

    /*
      Paper talks about whole molecule rotation such that the NorthPoles of each sphere aren't inside any other sphere.
     */



    if (debug) {
        cout << endl << "-->  DISCOVER COLLIDING SPHERES  <--" << endl << endl;
    }

    /****************************************************
         STEP 1: Compute neighbors for each atom

     *****************************************************/
    insideSphere.assign(atoms.size(), false); // store spheres that are completely engulfed in other spheres..
    collidingSpheres.clear();
    collidingSpheres.resize(atoms.size());
    for (uint i = 0; i < atoms.size(); i++) {


        if (!atoms(i).getActive()) continue;

        // Skip "ENGULFED" spheres
        if (insideSphere[i]) continue;


        double radiiI = atomicRadiiSurfaceAreaAndVolume[atoms[i]][0];

        // Check for neighbors...
        for (uint j = i + 1; j < atoms.size(); j++) {
            if (!atoms(j).getActive()) continue;

            double radiiJ = atomicRadiiSurfaceAreaAndVolume[atoms[j]][0];

            double dist = atoms(i).distance(atoms(j));

            if (dist < radiiI + radiiJ) {


                if (debug) {
                    fprintf(stdout, "Atom %4d is colliding with Atom %4d (and visa-versa)\n", i, j);
                }

                // STORE WHICH SPHERES ARE COLLIDING
                collidingSpheres[i].push_back(j);
                collidingSpheres[j].push_back(i);


                // ENGULFMENT CHECK on both i and j
                if (dist + radiiI < radiiJ) {
                    insideSphere[i] = true;
                }
                if (dist + radiiJ < radiiI) {
                    insideSphere[j] = true;
                }
            }

        }
    }

    filterEngulfedAtoms(atoms);

    rotatedMolecule = rotateMolecule(atoms);

    createStereographicProjectedCircles(atoms);
    getIntersectingAngles(atoms);
    getArcs(atoms);
    integrateArcs(atoms);

    if (rotatedMolecule) {
        atoms.applySavedCoor("preSurfaceArea");
    }

}

void SurfaceAreaAndVolume::addAtomsAndCharmmRadii(AtomPointerVector _atoms, CharmmParameterReader &_par) {

    atoms.clear();
    atomicRadiiSurfaceAreaAndVolume.clear();
    for (AtomPointerVector::iterator atomI = _atoms.begin(); atomI < _atoms.end(); atomI++) {
        string atomItype = (*atomI)->getType();
        //vector<double> vdwParam   = _par.vdwParam(atomItype);
        vector<double> vdwParam;
        // add a check for the return of the following bool function if  the param exist???
        _par.vdwParam(vdwParam, atomItype);
        atomicRadiiSurfaceAreaAndVolume[(*atomI)].push_back(vdwParam[1] + probeRadius); // RADIUS
        atomicRadiiSurfaceAreaAndVolume[(*atomI)].push_back(0.0); // SURFACE AREA
        atomicRadiiSurfaceAreaAndVolume[(*atomI)].push_back(0.0); // VOLUME

        atoms.push_back((*atomI));
    }

}

string SurfaceAreaAndVolume::toString(AtomPointerVector &_atoms) {

    stringstream ss;
    ss << "from pymol.cgo import *" << endl;
    ss << "from pymol import cmd" << endl;
    ss << "from pymol.vfont import plain" << endl;
    ss << "from glob import glob" << endl;

    double atomRadius = 2;
    for (uint i = 0; i < _atoms.size(); i++) {

        CartesianPoint southPole = _atoms(i).getCoor();
        southPole[2] -= atomRadius;
        CartesianPoint northPole = _atoms(i).getCoor();
        northPole[2] += atomRadius;

        double x = _atoms(i).getCoor()[0];
        double y = _atoms(i).getCoor()[1];
        double z = _atoms(i).getCoor()[2];

        stringstream s;
        s << "NP-" << i;
        //ss << northPole.toPyMol(s.str())<<endl;
        s.str("");
        s << "SP-" << i;
        //ss << southPole.toPyMol(s.str())<<endl;


        for (uint j = 0; j < tsCircles[i].size(); j++) {
            //double x2 = _atoms(collidingSpheres[i][j]).getCoor()[0];
            //double y2 = _atoms(collidingSpheres[i][j]).getCoor()[1];
            //double z2 = _atoms(collidingSpheres[i][j]).getCoor()[2];			


            for (uint k = 0; k < 360; k += 10) {
                double pt = tsCircles[i][j][0] + tsCircles[i][j][2] * cos(k);
                double ps = tsCircles[i][j][1] + tsCircles[i][j][2] * sin(k);

                CartesianPoint surfacePoint;
                surfacePoint[0] = x + (4 * atomRadius * atomRadius * pt) / (pt * pt + ps * ps + 4 * atomRadius * atomRadius);
                surfacePoint[1] = y + (4 * atomRadius * atomRadius * ps) / (pt * pt + ps * ps + 4 * atomRadius * atomRadius);
                surfacePoint[2] = z + atomRadius - (8 * atomRadius * atomRadius * atomRadius / (pt * pt + ps * ps + 4 * atomRadius * atomRadius));

                s.str("");
                s << "Surface" << i << "j" << j << "k" << k;
                //ss << surfacePoint.toPyMol(s.str())<<endl;


                double dist = sqrt(pt * pt + ps * ps);

                CartesianPoint vecToCircle = (surfacePoint - northPole).getUnit();
                vecToCircle *= dist;

                CartesianPoint circleCenter;
                circleCenter = northPole;
                circleCenter += vecToCircle;


                CartesianPoint circleNormal;
                circleNormal = (northPole - southPole).getUnit();
                circleNormal *= 0.05;


                CartesianPoint diskBottom;
                diskBottom = circleCenter - circleNormal;

                CartesianPoint diskTop;
                diskTop = circleCenter + circleNormal;


                s.str("");
                s << "centerCircle" << i << "j" << j << "k" << k;
                //ss <<circleCenter.toPyMol(s.str())<<endl;


                ss << "obj" << i << "j" << j << "k" << k << " = [ CYLINDER, " << diskBottom[0] << "," << diskBottom[1] << "," << diskBottom[2] << " ," << diskTop[0] << "," << diskTop[1] << "," << diskTop[2] << ",2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 ]" << endl;
                ss << "cmd.load_cgo(obj" << i << "j" << j << "k" << k << ",\"Circle" << i << "j" << j << "k" << k << "\")" << endl;

            }

            /*
            pt = 0;
            ps = 0;
				
            CartesianPoint surfacePoint;
            surfacePoint[0] = x + (4 * atomRadius * atomRadius * pt) / ( pt * pt + ps *ps + 4 * atomRadius * atomRadius);
            surfacePoint[1] = y + (4 * atomRadius * atomRadius * ps) / ( pt * pt + ps *ps + 4 * atomRadius * atomRadius);
            surfacePoint[2] = z + atomRadius - (8*atomRadius*atomRadius*atomRadius / ( pt * pt + ps *ps + 4 * atomRadius * atomRadius));
			
            double dist = sqrt( tsCircles[i][j][0] * tsCircles[i][j][0] + tsCircles[i][j][1] * tsCircles[i][j][1]);
            dist = sqrt(dist*dist +  (2*atomRadius)*(2*atomRadius));
            CartesianPoint vecToCircle = (surfacePoint - northPole).getUnit();
            vecToCircle *= dist;

            s.str("");
            s << "VecToCir"<<i<<"j"<<j;
            ss <<vecToCircle.toPyMol(s.str())<<endl;
            CartesianPoint circleCenter;
            circleCenter = northPole;
            circleCenter += vecToCircle;

            s.str("");
            s << "centerCircle"<<i<<"j"<<j;
            ss <<circleCenter.toPyMol(s.str())<<endl;

            CartesianPoint circleNormal;
            circleNormal = (northPole - southPole).getUnit();
            circleNormal *= 0.05;
            s.str("");
            s << "cNorm"<<i<<"j"<<j;
            ss <<circleNormal.toPyMol(s.str())<<endl;

            CartesianPoint diskBottom;
            diskBottom = circleCenter - circleNormal;

            CartesianPoint diskTop;
            diskTop   = circleCenter + circleNormal;

            //cout << "Points: "<<endl<<circleCenter<<endl<<circleNormal<<endl<<diskBottom<<endl<<diskTop<<endl<<surfacePoint<<endl;

            ss << "obj"<<i<<"c"<<j<<" = [ CYLINDER, "<<diskBottom[0]<<","<<diskBottom[1]<<","<<diskBottom[2]<<" ,"<<diskTop[0]<<","<<diskTop[1]<<","<<diskTop[2]<<",2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 ]"<<endl;
			
            ss << "cmd.load_cgo(obj"<<i<<"c"<<j<<",\"Circle"<<i<<"c"<<j<<"\")"<<endl;
             */


        }
    }


    return ss.str();
}

void SurfaceAreaAndVolume::filterEngulfedAtoms(AtomPointerVector &_atoms) {

    for (uint i = 0; i < _atoms.size(); i++) {
        int preSize = collidingSpheres[i].size();
        if (insideSphere[i]) {
            collidingSpheres[i].clear();
            continue;
        }

        vector<int>::iterator it;
        for (it = collidingSpheres[i].begin(); it != collidingSpheres[i].end(); it++) {

            if (insideSphere[(*it)]) {
                collidingSpheres[i].erase(it);
                it--;
            }
        }


        if (debug) {
            cout << "FILTER ATOM[" << i << "] for inside spheres: " << preSize << " colliding spheres has been reduced to " << collidingSpheres[i].size() << endl;
        }
    }
}

bool SurfaceAreaAndVolume::rotateMolecule(AtomPointerVector &_atoms) {


    // Rotate molecules
    bool rotatedMolecule = false;
    bool northPoleInside = true;

    RandomNumberGenerator rng;
    //rng.setRNGTimeBasedSeed();
    while (northPoleInside) {

        double minDist = MslTools::doubleMax;
        for (uint i = 0; i < _atoms.size(); i++) {

            double radiiI = atomicRadiiSurfaceAreaAndVolume[_atoms[i]][0];

            // SKIP ENGULFED SPHERES
            if (insideSphere[i]) continue;

            // CHECK NORTH POLE TO ALL OTHER ATOM CENTERS
            CartesianPoint northPole = _atoms(i).getCoor();
            northPole[2] += radiiI;
            for (uint j = 0; j < collidingSpheres[i].size(); j++) {

                double radiiJ = atomicRadiiSurfaceAreaAndVolume[_atoms[collidingSpheres[i][j]]][0];

                double dist = northPole.distance(_atoms[collidingSpheres[i][j]]->getCoor());

                if (fabs(dist - radiiJ) < 0.0001) {
                    if (debug) {
                        cout << "DIST of NP[" << i << "] to atom[" << collidingSpheres[i][j] << "]: " << dist << " and RADII: " << radiiJ << endl;
                        cout << _atoms(i) << endl << _atoms(collidingSpheres[i][j]) << endl << "RADII: " << radiiI << " " << radiiJ << endl;
                    }
                    minDist = dist;
                }

            }


        }

        // FOUND A BAD NORTH POLE...MUST ROTATE WHOLE MOLECULE
        if (minDist != MslTools::doubleMax) {


            CartesianPoint axis(rng.getRandomDouble(), rng.getRandomDouble(), rng.getRandomDouble());
            double angle = MslTools::mod(rng.getRandomDouble()*1000, 360);
            Matrix randRotMatrix = CartesianGeometry::getRotationMatrix(angle, axis);

            //cout << "MOLECULE NEEDS ROTATING \t "<<axis<<" degrees "<<angle<<"\n";			
            Transforms tr;
            tr.rotate(_atoms, randRotMatrix);


            rotatedMolecule = true;


        } else {
            northPoleInside = false;
        }



    }

    return rotatedMolecule;
}

void SurfaceAreaAndVolume::createStereographicProjectedCircles(AtomPointerVector &_atoms, bool _debug) {

    tsCircles.clear();
    tsCircles.resize(_atoms.size());

    // PROCESS EACH ATOM
    for (uint i = 0; i < _atoms.size(); i++) {

        double radiiI = atomicRadiiSurfaceAreaAndVolume[_atoms[i]][0];

        // SKIP ENGULFED
        if (insideSphere[i]) continue;

        if (_debug || debug) {
            fprintf(stdout, "ATOM %4d [%p] has %4d colliding spheres\n", i, _atoms[i], int(collidingSpheres[i].size()));
        }


        // FOR EACH SPHERE COLLIDING WITH ATOM i
        CartesianPoint &a1 = _atoms(i).getCoor();
        tsCircles[i].resize(collidingSpheres[i].size());

        for (uint j = 0; j < collidingSpheres[i].size(); j++) {

            double radiiJ = atomicRadiiSurfaceAreaAndVolume[_atoms[collidingSpheres[i][j]]][0];

            // SKIP ENGULFED
            if (insideSphere[collidingSpheres[i][j]]) continue;

            if (_debug || debug) {
                fprintf(stdout, "\tChecking collision with atom %4d ; %8.3f , %8.3f\n\n", collidingSpheres[i][j], radiiI, radiiJ);
            }

            // STEROGRAPHIC PROJECTION MATH...
            CartesianPoint &a2 = _atoms(collidingSpheres[i][j]).getCoor();
            double dx = a1[0] - a2[0];
            double dy = a1[1] - a2[1];
            double zp = a1[2] + radiiI - a2[2];
            double zm = a1[2] - radiiI - a2[2];
            double a = dx * dx + dy * dy + zp * zp - radiiJ*radiiJ;
            double b = 8 * radiiI * radiiI*dx;
            double c = 8 * radiiI * radiiI*dy;
            double d = 4 * radiiI * radiiI * (dx * dx + dy * dy + zm * zm - radiiJ * radiiJ);


            /*
              STORE CIRCLE DEFINITION
              t = t0+r0*cos(theta)
              s = s0+r0*sin(theta)
             */
            tsCircles[i][j].resize(4);
            tsCircles[i][j][0] = -b / (2 * a); // t0
            tsCircles[i][j][1] = -c / (2 * a); // s0
            tsCircles[i][j][2] = sqrt((b * b + c * c - 4 * a * d) / (4 * a * a)); // r0
            tsCircles[i][j][3] = 1; // orientation
            if (a > 0) {
                tsCircles[i][j][3] = -1;
            }
        }
    }
}

void SurfaceAreaAndVolume::getIntersectingAngles(AtomPointerVector &_atoms, bool _debug) {

    intersectionAngles.clear();
    intersectionAngles.resize(_atoms.size());


    // FOR EACH ATOM
    for (uint i = 0; i < _atoms.size(); i++) {


        // SKIP ENGULFED ATOMS
        if (insideSphere[i]) continue;
        if (_debug || debug) {
            fprintf(stdout, "ATOM %4d has %4d circles\n", i, int(tsCircles[i].size()));
        }
        intersectionAngles[i].resize(tsCircles[i].size());

        // Only 1 circle, full circumference intersecting angles are 0 and 2PI ( sign change for orientation of given circle)
        if (tsCircles[i].size() == 1) {

            intersectionAngles[i][0].push_back(0.0);
            intersectionAngles[i][0].push_back(2 * M_PI * tsCircles[i][0][3]);

            continue;
        }


        // PROCESS EACH CIRCLE FOR INTERSECTIONS
        for (uint j = 0; j < tsCircles[i].size(); j++) {

            double t1 = tsCircles[i][j][0];
            double s1 = tsCircles[i][j][1];
            double r1 = tsCircles[i][j][2];
            if (_debug || debug) {
                fprintf(stdout, "\tCIRCLE1 %4d [ %8.3f , %8.3f , %8.3f, %8.3f ]\n", j, t1, s1, r1, tsCircles[i][j][3]);
            }

            // CHECK VS ALL OTHER CIRCLES
            for (uint k = 0; k < tsCircles[i].size(); k++) {

                if (j == k) continue;

                double t2 = tsCircles[i][k][0];
                double s2 = tsCircles[i][k][1];
                double r2 = tsCircles[i][k][2];
                if (_debug || debug) {
                    fprintf(stdout, "\tCIRCLE2 %4d [ %8.3f , %8.3f , %8.3f, %8.3f ]\n", k, t2, s2, r2, tsCircles[i][k][3]);
                }

                /*
				 
                  /|\
              r1 / | \ r2
                /  |  \
               /   |   \
     (t1,s1)  /____|____\  (t2,s2)

                x   dist-x

             [---dist----]


                 */

                double dist = sqrt((t1 - t2) * (t1 - t2) + (s1 - s2) * (s1 - s2));
                double x = (r1 * r1 - r2 * r2 + dist * dist) / (2 * dist);
                double alpha = acos(x / r1);
                //double beta  = acos( (dist - x) / r2);


                // Test for 2 points of intersection...

                // TEST 1 : Non-intersecting circles
                if (r1 + r2 < dist) {
                    if (_debug || debug) {
                        //cout << "NON-INTESECTING CIRCLES\n";
                    }
                    continue;
                }

                // TEST 2: Enclosed circle
                if (r1 + r2 > dist && (dist + r1 <= r2 || dist + r2 <= r1)) {
                    if (_debug || debug) {
                        //cout << "ENCLOSED CIRCLE\n";
                    }
                    continue;
                }

                // TEST 3: Single Point
                if (r1 + r2 == dist) {
                    if (_debug || debug) {
                        //cout << "SINGLE POINT\n";
                    }
                    intersectionAngles[i][j].push_back(0.0);

                    continue;
                }

                // TEST 4: Two points of intersection
                if ((r1 + r2 > dist)) {
                    //double t = t1 +r1 *cos(alpha);
                    //double s = s1 +r1 *sin(alpha);

                    // Find relative angle...  between (t1+1,s1), (t1,s1), (t2,s2)
                    CartesianPoint p1(t1 + 1, s1, 0);
                    CartesianPoint p2(t1, s1, 0);
                    CartesianPoint p3(t2, s2, 0);

                    double angle = (M_PI / 180) * p1.angle(p2, p3);
                    if (s1 > s2) {
                        angle = 2 * M_PI - angle;
                    }

                    double tmpAngle = MslTools::mod((alpha + angle), (2.0 * M_PI));
                    intersectionAngles[i][j].push_back(tmpAngle);
                    if (_debug || debug) {
                        fprintf(stdout, "\t\tAngle is %8.3f, mod is %8.3f", (alpha + angle), tmpAngle);
                    }

                    tmpAngle = MslTools::mod(2 * M_PI - alpha + angle, (2.0 * M_PI));
                    intersectionAngles[i][j].push_back(tmpAngle);
                    if (_debug || debug) {
                        fprintf(stdout, "\t\tAngle is %8.3f, mod is %8.3f\n", (2 * M_PI - alpha + angle), tmpAngle);
                    }
                }
            }
        }
    }
}

bool SurfaceAreaAndVolume::insideOutsideTest(uint i, uint j) {
    bool inCircle = true;

    for (uint k = 0; k < tsCircles[i].size(); k++) {
        if (j == k)
            continue;

        double a = tsCircles[i][j][0] + tsCircles[i][j][2] - tsCircles[i][k][0];
        double b = tsCircles[i][j][1] - tsCircles[i][k][1];
        double d = sqrt(a * a + b * b);

        if (d < tsCircles[i][k][2]) {
            if (tsCircles[i][k][3] > 0)
                inCircle = inCircle && true;
            else
                inCircle = inCircle && false;
        } else if (d > tsCircles[i][k][2]) {
            if (tsCircles[i][k][3] > 0)
                inCircle = inCircle && false;
            else
                inCircle = inCircle && true;
        } else {
            a = tsCircles[i][j][0] - tsCircles[i][k][0];
            b = tsCircles[i][j][1] - tsCircles[i][k][1];
            d = sqrt(a * a + b * b);
            if (d < tsCircles[i][k][2]) {
                if (tsCircles[i][k][3] > 0)
                    inCircle = inCircle && true;
                else
                    inCircle = inCircle && false;
            } else {
                if (tsCircles[i][k][3] > 0)
                    inCircle = inCircle && false;
                else
                    inCircle = inCircle && true;
            }
        }
    }
    return inCircle;
}

void SurfaceAreaAndVolume::getArcs(AtomPointerVector &_atoms, bool _debug) {
    arcs.clear();
    arcs.resize(_atoms.size());

    // FOR EACH ATOM
    for (uint i = 0; i < _atoms.size(); i++) {

        // SKIP ENGULFED ATOMS
        if (insideSphere[i]) continue;

        // SKIP ATOMS WITH NO CIRCLES
        // Should this really be the case?  
        // We seem to take care of no intersection circles
        // below (see FULL CIRCLE ARC comment).
        if (tsCircles[i].size() == 0) continue;

        arcs[i].resize(tsCircles[i].size());

        if (_debug || debug) {
            fprintf(stdout, "ATOM %4d has %4d circles\n", i, int(tsCircles[i].size()));
        }


        // FULL CIRCLE ARC WHEN ONLY ONE CIRCLE
        // NOTE, we should never get here due to the continue
        // statement a few lines above!
        if (tsCircles[i].size() + 1 == 1) {
            arcs[i][0].push_back(0);
            arcs[i][0].push_back(2 * M_PI * tsCircles[i][0][3]);
            continue;
        }

        for (uint j = 0; j < tsCircles[i].size(); j++) {

            // PARAMETERS FOR THIS CIRCLE
            double t0 = tsCircles[i][j][0];
            double s0 = tsCircles[i][j][1];
            double r0 = tsCircles[i][j][2];

            if (_debug || debug) {
                fprintf(stdout, "\tCIRCLE %4d [ %8.3f , %8.3f , %8.3f, %8.3f ]\n", j, t0, s0, r0, tsCircles[i][j][3]);
            }

            if (intersectionAngles[i][j].size() == 0) {
                if (_debug || debug) {
                    fprintf(stdout, "\t\tARC FULL CIRCLE [ -2pi ; 0 ]\n");
                }
                // BTH
                bool passed = insideOutsideTest(i, j);

                if (passed) {
                    arcs[i][j].push_back(0);
                    arcs[i][j].push_back(2 * M_PI * tsCircles[i][j][3]);
                }
                continue;
            }

            // Sort depenendent on circle orientation..
            if (tsCircles[i][j][3] > 0) {
                std::sort(intersectionAngles[i][j].begin(), intersectionAngles[i][j].end());
            } else {
                std::sort(intersectionAngles[i][j].begin(), intersectionAngles[i][j].end(), std::greater<double>());
            }

            int validAngleCount = 0; // KEEP TRACK OF NUMBER OF VALID ANGLES IN THIS CIRCLE
            //double firstAngle = 0.0;
            //double lastAngle = 2 * M_PI;

            for (uint k = 0; k < intersectionAngles[i][j].size(); k++) {
                if (_debug || debug) {
                    cout << "\t\tIntersecting angle: " << k << " : " << intersectionAngles[i][j][k] << endl;
                }

                double t = t0 + r0 * cos(intersectionAngles[i][j][k]);
                double s = s0 + r0 * sin(intersectionAngles[i][j][k]);

                // Check to make sure not inside any other circle.
                bool midPointInside = false;
                for (uint jj = 0; jj < tsCircles[i].size(); jj++) {
                    if (j == jj) continue;

                    // Check to see if arc midpoint is inside a circle...
                    if (k != 0) {
                        double tmid = t0 + r0 * cos((intersectionAngles[i][j][k] + intersectionAngles[i][j][k - 1]) / 2.0);
                        double smid = s0 + r0 * sin((intersectionAngles[i][j][k] + intersectionAngles[i][j][k - 1]) / 2.0);
                        double dist = sqrt((tmid - tsCircles[i][jj][0]) * (tmid - tsCircles[i][jj][0]) + (smid - tsCircles[i][jj][1]) * (smid - tsCircles[i][jj][1]));
                        //cout  << "****** ANGLES: "<<intersectionAngles[i][j][k]<<" , "<<intersectionAngles[i][j][k-1]<<" t0,s0,r0 "<<t0<<","<<s0<<","<<r0<<endl;
                        if ((dist < tsCircles[i][jj][2]) && (fabs(dist - tsCircles[i][jj][2]) > 0.00001)) {

                            if (tsCircles[i][jj][3] < 0) {
                                midPointInside = true;


                                if (_debug || debug) {
                                    fprintf(stdout, "\t\t\tMIDPOINT INSIDE CIRCLE(%4d) WITH DIST %8.3f AND RADIUS %8.3f\n", jj, dist, tsCircles[i][jj][2]);
                                }
                            } else {
                                if (_debug || debug) {
                                    fprintf(stdout, "\t\t\tMIDPOINT2 INSIDE CIRCLE(%4d) WITH DIST %8.3f AND RADIUS %8.3f but POS ORIEN\n", jj, dist, tsCircles[i][jj][2]);
                                }
                            }
                        } else {
                            if (tsCircles[i][jj][3] > 0) {
                                midPointInside = true;
                            }
                        }
                    }
                }


                if (_debug || debug) {
                    cout << endl;
                }

                if (midPointInside) {
                    continue;
                }

                //firstAngle = intersectionAngles[i][j][k];
                //lastAngle = intersectionAngles[i][j][k];
                validAngleCount++;
                if (validAngleCount > 1) {

                    // Portion of arc that is not overlapping:
                    //  intersectionAngles[i][j][k-1] -> intersectionAngles[i][j][k];
                    arcs[i][j].push_back(intersectionAngles[i][j][k - 1]);
                    arcs[i][j].push_back(intersectionAngles[i][j][k]);

                    if (_debug || debug) {
                        double tm1 = t0 + r0 * cos(intersectionAngles[i][j][k - 1]);
                        double sm1 = s0 + r0 * sin(intersectionAngles[i][j][k - 1]);
                        fprintf(stdout, "\t\tARC [ %8.3f , %8.3f ] ; [ %8.3f , %8.3f ] --> t,s (%8.3f, %8.3f) (%8.3f, %8.3f)\n", intersectionAngles[i][j][k - 1], intersectionAngles[i][j][k], intersectionAngles[i][j][k - 1], intersectionAngles[i][j][k], tm1, sm1, t, s);
                    }
                }
            }

            // BTH
            double angle = intersectionAngles[i][j][0] + 2 * M_PI + intersectionAngles[i][j][intersectionAngles[i][j].size() - 1];
            // BTH I added the /2.0 below.  It seems to make sense to me...
            double tmid = t0 + r0 * cos(angle / 2.0);
            double smid = s0 + r0 * sin(angle / 2.0);
            bool midPointInside = false;
            for (uint jj = 0; jj < tsCircles[i].size(); jj++) {
                if (j == jj) continue;

                double dist = sqrt((tmid - tsCircles[i][jj][0]) * (tmid - tsCircles[i][jj][0]) + (smid - tsCircles[i][jj][1]) * (smid - tsCircles[i][jj][1]));
                if ((dist < tsCircles[i][jj][2]) && (fabs(dist - tsCircles[i][jj][2]) > 0.00001)) {

                    if (tsCircles[i][jj][3] < 0) {
                        midPointInside = true;

                        if (_debug || debug) {
                            fprintf(stdout, "\t\t\tMIDPOINT2 INSIDE CIRCLE(%4d) WITH DIST %8.3f AND RADIUS %8.3f\n", jj, dist, tsCircles[i][jj][2]);
                        }
                    }
                } else {
                    if (tsCircles[i][jj][3] > 0) {
                        midPointInside = true;
                    }
                }
            }

            if (midPointInside) continue;

            arcs[i][j].push_back(intersectionAngles[i][j][intersectionAngles[i][j].size() - 1]);
            arcs[i][j].push_back(intersectionAngles[i][j][0] + 2 * M_PI * tsCircles[i][j][3]);

            if (_debug || debug) {
                fprintf(stdout, "\t\tARC0 [ %8.3f %8.3f ] %8.3f \n", intersectionAngles[i][j][intersectionAngles[i][j].size() - 1], intersectionAngles[i][j][0] + 2 * M_PI * tsCircles[i][j][3], intersectionAngles[i][j][0]);
            }

            //arcs[i][j].push_back(-2*M_PI);
            //arcs[i][j].push_back(0);

            //arcs[i][j].push_back(lastAngle);
            //arcs[i][j].push_back(firstAngle);

            if (_debug || debug) {
                //fprintf(stdout,"\tARC1 [ %8.3f %8.3f ] \n",firstAngle,lastAngle);
            }
        }
    }
}

void SurfaceAreaAndVolume::integrateArcs(AtomPointerVector &_atoms, bool _debug) {
    double totalVolume = 0.0;
    double totalSurfaceArea = 0.0;
    for (uint i = 0; i < _atoms.size(); i++) {

        double radiiI = atomicRadiiSurfaceAreaAndVolume[_atoms[i]][0];
        if (insideSphere[i]) continue;
        // if (tsCircles[i].size() == 0) continue;
        // BTH
        if(tsCircles[i].size() == 0) {
            double volume = 4.0/3.0*M_PI * radiiI*radiiI*radiiI;
            double area = 4.0*M_PI * radiiI*radiiI;
            totalVolume += volume;
            totalSurfaceArea += area;
            atomicRadiiSurfaceAreaAndVolume[_atoms[i]][1] = area;
            atomicRadiiSurfaceAreaAndVolume[_atoms[i]][2] = volume;
            continue;
        }
        if (_debug || debug) {

            fprintf(stdout, "ATOM %4d has %4d circles\n", i, int(tsCircles[i].size()));
        }


        double volume = 0.0;
        double surfaceArea = 0.0;
        //double z = _atoms(i).getZ();
        bool positiveOrientation = false;
        for (uint j = 0; j < tsCircles[i].size(); j++) {


            if (tsCircles[i][j][3] > 0) {
                positiveOrientation = true;
            }

            if (_debug || debug) {
                fprintf(stdout, "\tCIRCLE %4d [ %8.3f , %8.3f , %8.3f, %8.3f ]\n", j, tsCircles[i][j][0], tsCircles[i][j][1], tsCircles[i][j][2], tsCircles[i][j][3]);
            }
            double A = (4 * radiiI * radiiI + tsCircles[i][j][0] * tsCircles[i][j][0] + tsCircles[i][j][1] * tsCircles[i][j][1] + tsCircles[i][j][2] * tsCircles[i][j][2]) / 2;
            double B = tsCircles[i][j][0] * tsCircles[i][j][2];
            double C = tsCircles[i][j][1] * tsCircles[i][j][2];
            double D = A * A - B * B - C*C;
            double Dsqrt = sqrt(D);

            for (uint a = 0; a < arcs[i][j].size(); a += 2) {

                double alpha = arcs[i][j][a];
                double beta = arcs[i][j][a + 1];

                bool flip = false;
                if (beta - alpha < 0) {
                    flip = true;
                    alpha = arcs[i][j][a + 1];
                    beta = arcs[i][j][a];

                }

                if (_debug || debug) {
                    fprintf(stdout, "\t\tARC %4d [ %8.3f , %8.3f ] \n", a, alpha, beta);
                }

                double sum = (beta + alpha) / 2.0;
                double diff = (beta - alpha) / 2.0;


                double I1 = 0.0;
                double I2 = 0.0;
                double I3 = 0.0;
                double J1 = 0.0;
                double J2 = 0.0;
                double J3 = 0.0;
                if (fabs(fabs(beta - alpha) - 2.0 * M_PI) > 0.000000001) {

                    I1 = (2.0 / Dsqrt) * ((M_PI / 2.0) - atan((A * cos(diff) + B * cos(sum) + C * sin(sum)) / (Dsqrt * sin(diff))));

                    I2 = (1.0 / D) * (A * I1 +  \
							         ((-B * sin(beta) + C * cos(beta)) / (A + B * cos(beta) + C * sin(beta))) - \
							         ((-B * sin(alpha) + C * cos(alpha)) / (A + B * cos(alpha) + C * sin(alpha))));

                    I3 = (1 / (2.0 * D)) * (\

                            ((-B * sin(beta) + C * cos(beta)) / pow(A + B * cos(beta) + C * sin(beta), 2)) - \
						              ((-B * sin(alpha) + C * cos(alpha)) / pow(A + B * cos(alpha) + C * sin(alpha), 2)) + \
							      ((((-B / A) * sin(beta) + (C / A) * cos(beta)) / (A + B * cos(beta) + C * sin(beta))) - \
								  (((-B / A) * sin(alpha) + (C / A) * cos(alpha)) / (A + B * cos(alpha) + C * sin(alpha)))) + \
						(((2 * A * A + B * B + C * C) * I2) / A));

                    J1 = (fabs(beta - alpha) + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) * I1) / 2.0;
                    J2 = (I1 + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) * I2) / 4.0;
                    J3 = (I2 + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) * I3) / 8.0;

                } else {
                    //cout << "FULL CIRCLE"<<endl;
                    I1 = 2 * M_PI / Dsqrt;
                    I2 = (2 * M_PI * A) / (Dsqrt * Dsqrt * Dsqrt);
                    I3 = (M_PI * (2 * A * A + B * B + C * C)) / (Dsqrt * Dsqrt * Dsqrt * Dsqrt * Dsqrt);

                    J1 = M_PI + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) / 2.0 * I1;
                    J2 = (I1 + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) * I2) / 4.0;
                    J3 = (I2 + (tsCircles[i][j][2] * tsCircles[i][j][2] - A) * I3) / 8.0;
                }

                double radii3 = radiiI * radiiI*radiiI;

                // DIRECTLY FROM THE PAPER
                //double volumeIntegral      = (   ( (128*J3*radii3*radii3*radiiI ) / 3) + 
                //                                 ( (2  *J1*radii3                  ) / 3) - 
                //                                 ( (8  *J2*radii3*radiiI        ) * (3*_atoms(i).getZ()+2*radiiI) / 3));

                //WORKS ALMOST: double volumeIntegral = ((128*J3*radii3*radii3*radiiI + 8*J2*radii3*radiiI*radiiI + 2 *J1*radii3) / 3) -8*radii3*radiiI*J2*(_atoms(i).getZ()+radiiI);
                double volumeIntegral = (128.0 * J3 * radii3 * radii3 * radiiI + 8.0 * J2 * radii3 * radiiI * radiiI + 2.0 * J1 * radii3) / 3.0 - 8.0 * radii3 * radiiI * J2 * (_atoms(i).getZ() + radiiI);


                if (_debug || debug) {
                    cout << "J1,J2,J3: " << J1 << "," << J2 << "," << J3 << "," << radiiI << "," << _atoms(i).getZ() << "," << volumeIntegral << endl;
                    cout << "  T1: " << 128 * J3 * radii3 * radii3 * radiiI << endl;
                    cout << "  T2: " << 8 * J2 * radii3 * radiiI * radiiI << endl;
                    cout << "  T3: " << 2 * J1*radii3;
                    cout << "  T4: " << 8 * radii3 * radiiI * J2 * (_atoms(i).getZ() + radiiI) << endl;
                    cout << "  Sum: " << ((128.0 * J3 * radii3 * radii3 * radiiI)+(8.0 * J2 * radii3 * radiiI * radiiI)+(2.0 * J1 * radii3)) / 3.0;
                    cout << "  ALL: " << ((128.0 * J3 * radii3 * radii3 * radiiI)+(8.0 * J2 * radii3 * radiiI * radiiI)+(2.0 * J1 * radii3)) / 3.0 - 8.0 * radii3 * radiiI * J2 * (_atoms(i).getZ() + radiiI);
                }
                double surfaceAreaIntegral = J1;

                // -1,+1 
                int sign = 1;
                //cout << "Alpha,Beta: "<<alpha<<","<<beta<<","<<(beta-alpha)<<endl;
                if (fabs(fabs(beta - alpha) - 2.0 * M_PI) < 0.0001) {
                    if (tsCircles[i][j][3] < 0.0) {
                        if (_debug || debug) {
                            cout << "SIGN CHANGE" << endl;
                        }
                        sign = -1;
                    }

                } else if (flip) {
                    if (_debug || debug) {
                        cout << "SIGN CHANGE" << endl;
                    }
                    sign = -1;
                }
                if (_debug || debug) {
                    fprintf(stdout, "\tPOS? %4d\n\tI: %8.3f, %8.3f, %8.3f\n\tJ: %8.3f, %8.3f, %8.3f\n\tdelta: %8.3f, %8.3f\n\trr: %8.3f\n", positiveOrientation, I1, I2, I3, J1, J2, J3, sign * 2 * radiiI * radiiI*surfaceAreaIntegral, sign*volumeIntegral, (tsCircles[i][j][2] * tsCircles[i][j][2] - A));
                }
                volume += sign*volumeIntegral;
                surfaceArea += sign * 2.0 * radiiI * radiiI*surfaceAreaIntegral;
                if (_debug || debug) {
                    cout << "V[" << j << "," << a << "] = " << volume << endl;
                }
            }
        }


        double atomicSurfaceArea = 0.0;
        double atomicVolume = 0.0;

        if (positiveOrientation) {

            atomicSurfaceArea = surfaceArea;
            atomicVolume = volume;


            if (_debug || debug) {
                fprintf(stdout, "\tVolume subtotal       : %16.3f  %16.3f ===> %16.3f\n", volume, 0.0, totalVolume);
                fprintf(stdout, "\tSurface Area subtotal : %16.3f  %16.3f ===> %16.3f\n", surfaceArea, 0.0, totalSurfaceArea);
            }

        } else {
            atomicSurfaceArea = 4.0 * M_PI * radiiI * radiiI + surfaceArea;
            atomicVolume = (4.0 / 3.0) * M_PI * radiiI * radiiI * radiiI + volume;

            if (_debug || debug) {
                fprintf(stdout, "\tVolume subtotal       : %16.3f  %16.3f ===> %16.3f\n", volume, (4.0 / 3.0) * M_PI * radiiI * radiiI * radiiI + volume, totalVolume);
                fprintf(stdout, "\tSurface Area subtotal : %16.3f  %16.3f ===> %16.3f\n", surfaceArea, 4.0 * M_PI * radiiI * radiiI + surfaceArea, totalSurfaceArea);
            }

        }

        atomicRadiiSurfaceAreaAndVolume[_atoms[i]][1] = atomicSurfaceArea;
        atomicRadiiSurfaceAreaAndVolume[_atoms[i]][2] = atomicVolume;


        totalSurfaceArea += atomicSurfaceArea;
        totalVolume += atomicVolume;
    }


    if (_debug || debug) {
        fprintf(stdout, "TOTALS : VOLUME = %8.3f Angs^3 ; SURFACE AREA = %8.3f Angs^2\n", totalVolume, totalSurfaceArea);
    }

    surfaceArea = totalSurfaceArea;
    volume = totalVolume;
}
