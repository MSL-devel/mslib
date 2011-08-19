/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
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

#ifndef SURFACEAREAANDVOLUME_H
#define SURFACEAREAANDVOLUME_H

#include "AtomPointerVector.h"
#include "CharmmParameterReader.h"
#include "Transforms.h"
#include <vector>

namespace MSL {

    class SurfaceAreaAndVolume {
    public:
        SurfaceAreaAndVolume();
        ~SurfaceAreaAndVolume();

        // Standalone functions
        void computeSurfaceAreaOccludingPoints(AtomPointerVector &_atoms, std::vector<double> &_radii);
        void computeSurfaceAreaAndVolume(AtomPointerVector &_atoms, std::vector<double> &_radii);

        // Member-Variable Dependent functions
        void addAtomsAndCharmmRadii(AtomPointerVector _atoms, CharmmParameterReader &_par); // Store std::map of AllAtoms and a Atomic+Probe radii, getting atomic radii from CharmmParameterReader
        void computeSurfaceAreaAndVolume(); // Compute All Atom Surface Area and Volume      ; does     set surfaceArea, Volume member variables


        std::vector<double> getRadiiSurfaceAreaAndVolume(Atom *_atom); // Get entry to atomicRadiiSurfaceAreaAndVolume

        double getSurfaceArea();
        double getVolume();

        double getProbeRadius();
        void setProbeRadius(double _probeRadius);

        bool getDebug();
        void setDebug(bool _debug);
        std::string toString(AtomPointerVector &_atoms);
    private:

        void setup();
        void copy(const SurfaceAreaAndVolume &_sav);


        /*

          Occluded Points Surface Area Helper Functions

         */
        double getExposedFraction(std::vector<std::pair<double, CartesianPoint> > &_hiddenPoints, double _radii);

        /*

          Stereographic Surface Area Helper Functions

         */

        void filterEngulfedAtoms(AtomPointerVector &_atoms);
        bool rotateMolecule(AtomPointerVector &_atoms);
        void createStereographicProjectedCircles(AtomPointerVector &_atoms, bool _debug = false);
        void getIntersectingAngles(AtomPointerVector &_atoms, bool _debug = false);
        void getArcs(AtomPointerVector &_atoms, bool _debug = false);
        bool insideOutsideTest(uint i, uint j);
        void integrateArcs(AtomPointerVector &_atoms, bool _debug = false);

        std::vector<bool> insideSphere;
        std::vector<std::vector<int> > collidingSpheres;
        std::vector<std::vector<std::vector<double> > > tsCircles;
        std::vector<std::vector<std::vector<double> > > arcs;
        std::vector<std::vector<std::vector< double> > > intersectionAngles;

        // Contains All Atoms of a system....
        AtomPointerVector atoms;
        std::map<Atom *, std::vector<double> > atomicRadiiSurfaceAreaAndVolume;



        //
        double probeRadius;
        double surfaceArea;
        double volume;
        bool debug;


    };

    inline double SurfaceAreaAndVolume::getSurfaceArea() {
        return surfaceArea;
    }

    inline double SurfaceAreaAndVolume::getVolume() {
        return volume;
    }

    inline bool SurfaceAreaAndVolume::getDebug() {
        return debug;
    }

    inline void SurfaceAreaAndVolume::setDebug(bool _debug) {
        debug = _debug;
    }

    inline double SurfaceAreaAndVolume::getProbeRadius() {
        return probeRadius;
    }

    inline void SurfaceAreaAndVolume::setProbeRadius(double _probeRadius) {
        probeRadius = _probeRadius;
    }

    inline std::vector<double> SurfaceAreaAndVolume::getRadiiSurfaceAreaAndVolume(Atom *_at) {
        return atomicRadiiSurfaceAreaAndVolume[_at];
    }
}

#endif
