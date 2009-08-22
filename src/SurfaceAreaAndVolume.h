/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#include "AtomVector.h"
#include "CharmmParameterReader.h"
#include <vector>

class SurfaceAreaAndVolume{
	public:
		SurfaceAreaAndVolume();
		~SurfaceAreaAndVolume();

		// Standalone functions
		void computeSurfaceAreaOccludingPoints(AtomVector &_atoms, vector<double> &_radii);
		void computeSurfaceAreaAndVolumeStereographicProjectIntegration(AtomVector &_atoms,vector<double> &_radii);
		void computeTest(AtomVector &_atoms,vector<double> &_radii);
		
		// Member-Variable Dependent functions
		void addAtomsAndCharmmRadii(AtomVector _atoms, CharmmParameterReader &_par);  // Store map of AllAtoms and a Atomic+Probe radii, getting atomic radii from CharmmParameterReader
		void computeSurfaceAreaAndVolumeStereographicProjectIntegration();                                 // Compute All Atom Surface Area and Volume      ; does     set surfaceArea, Volume member variables
		//void computeSurfaceAreaAndVolumeStereographicProjectIntegration(Atom &_atom);                      // Compute a single Atom Surface Area and Volume ; does NOT set surfaceArea, Volume member variables


		vector<double> getRadiiSurfaceAreaAndVolume(Atom *_atom); // Get entry to atomicRadiiSurfaceAreaAndVolume

		double getSurfaceArea();
		double getVolume();

		double getProbeRadius();
		void   setProbeRadius(double _probeRadius);

		bool getDebug();
		void setDebug(bool _debug);
		string toString(AtomVector &_atoms);
	private:

		void setup();
		void copy(const SurfaceAreaAndVolume &_sav);


		/*

		  Occluded Points Surface Area Helper Functions

		 */
		double getExposedFraction(vector<pair<double,CartesianPoint> > &_hiddenPoints,double _radii);

		/*

		  Stereographic Surface Area Helper Functions

		 */

		void filterEngulfedAtoms(AtomVector &_atoms);
		bool rotateMolecule(AtomVector &_atoms);
		void createStereographicProjectedCircles(AtomVector &_atoms,bool _debug=false);
		void getIntersectingAngles(AtomVector &_atoms,bool _debug=false);
		void getArcs(AtomVector &_atoms,bool _debug=false);
		void integrateArcs(AtomVector &_atoms, bool _debug=false);

		vector<bool> insideSphere;
		vector<vector<int> > collidingSpheres;
		vector<vector<vector<double> > > tsCircles;
		vector<vector<vector<double> > > arcs;
		vector<vector<vector< double> > > intersectionAngles;

		// Contains All Atoms of a system....
		AtomVector atoms;
		map<Atom *, vector<double> > atomicRadiiSurfaceAreaAndVolume;
		


		//
		double probeRadius;
		double surfaceArea;
		double volume;
		bool debug;

		
};

inline double SurfaceAreaAndVolume::getSurfaceArea(){ return surfaceArea; }
inline double SurfaceAreaAndVolume::getVolume()     { return volume;      }

inline bool SurfaceAreaAndVolume::getDebug()        { return debug;       }
inline void SurfaceAreaAndVolume::setDebug(bool _debug) { debug = _debug; }

inline double SurfaceAreaAndVolume::getProbeRadius()        { return probeRadius;   }
inline void SurfaceAreaAndVolume::setProbeRadius(double _probeRadius) { probeRadius = _probeRadius; }

inline vector<double> SurfaceAreaAndVolume::getRadiiSurfaceAreaAndVolume(Atom *_at) { return atomicRadiiSurfaceAreaAndVolume[_at]; }
#endif
