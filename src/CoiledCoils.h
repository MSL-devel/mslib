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

#ifndef COILEDCOILS_H
#define COILEDCOILS_H



// STL Includes


// MSL Includes
#include "CartesianPoint.h"
#include "AtomPointerVector.h"
//#include "BBQTableReader.h"
//#include "BBQTableWriter.h"
//#include "BBQTable.h"
#include "System.h"
#include "Atom.h"
#include "Symmetry.h"
//#include "System.h"
#include "MslTools.h"

// BOOST Includes

/*************************************************************************
 *   NOTE:  REQUIRES COPY CONSTRUCTOR AND = OPERATOR!!!
 *
 *   It will be difficult until we have a system copy constructor.
 *   Can we generate a sys on demand or do we need the getSystem()
 *   call (same goes for Chain)
 ************************************************************************/

/*************************************************************************
 * Conversion from Gevorg/Crick to Dan/North notation:
 *
 * Gevorg/Crick:
 * r0 = super helical radius (Angstroms)
 * r1 = alpha helical radius (Angstroms)
 * w0 = super helical frequency (radians/residue)
 * w1 = alpha helical frequency (radians/residue)
 * a = alpha = helix crossing angle (radians)
 * p1 = alpha helical phase (radians)
 * dZ = Z offset (initially at 0) (angstroms)
 *
 * Dan/North:
 * r0 = super helical radius (Angstroms)
 * r1 = alpha helical radius (Angstroms)
 * theta = alpha helical frequency (radians/residue)
 * helicalPhase =  phi1 = alpha helical phase (radians)
 * h = rise/residue in alpha helix (Angstroms)
 * pitch = super helical pitch, distance between repeats in superhelix (Angstroms)
 *
 * Conversions:
 *
 * North		Crick	
 * ================================================
 * r0		=	r0
 * pitch	~	a	=   -atan(2*PI*r0/pitch)
 * risePerRes	~	w0	=   ( sin(a)*(risePerRes) ) / r0
 * r1		=	r1
 * w1           =	w1
 * phi1 	=	phi1
 * nRes		=	nRes
 * dZ		=	dZ
 *
 * 
 * wa * t2 	=       w1
 * t2           =       w0
 *
 ************************************************************************/

namespace MSL { 
class CoiledCoils {

	public:
		CoiledCoils();
		~CoiledCoils();

		bool primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  std::string _symmetry, int _N, System& _sys, std::vector<std::string> _startingPositions);
		AtomPointerVector& getCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes);
		AtomPointerVector& getCoiledCoilCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes);
		AtomPointerVector& getCoiledCoilBundle(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  std::string _symmetry, int _N);
		AtomPointerVector& getCoiledCoilBundleCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes, std::string _symmetry, int _N);
		//void applyCoiledCoil(AtomPointerVector &_av, double _p0);
		//void offersCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1);
		//void sotoCoiledCoils(double _r0, double _risePerRes, double _r1, int _nRes, double _resPerTurn, double _alpha, double _helicalPhase);

		AtomPointerVector& getAtomPointers();
		//Chain& getChain(); // include all backbone atoms
		//System& getSystem(); // include all backbone atoms

		//void useBBQTable(std::string _bbqTable);

		/*
		  Varing parameters as a function of residue number
		  
		  Need to specify:
		     1. Starting residue
		     2. Ending Residue
		     3. Start Params
		     4. End Params
		     5. Param Function
		 */

		//void addNorthParamSchedule(std::vector<int> _residues, std::vector<double> _t, std::vector<double> _r1, std::vector<double> _wAlpha, std::vector<double> _helicalPhase, std::vector<double> _r0, std::vector<double> _p0, std::vector<int> _paramFunc);

	private:
		
		System *sys;

		void radCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes);

		// North Parameter Schedule Variables
		//std::vector<int> npStartStopResidues;
		//std::vector<double> npPhase;
		//std::vector<double> npR1;
		//std::vector<double> npWAlpha;
		//std::vector<double> npHelicalPhase;
		//std::vector<double> npR0;
		//std::vector<double> npP0;
		//std::vector<int> npParamFuncs;
		
		void generateNorthCoor(Atom *_pAtom, unsigned int i, double _r1, double _wa, double _phi1, double _r0, double _pitch, double _risePerRes, double _L, double _corr, double _zCenter); 
		AtomPointerVector atoms;
		//void generateGevorgCoor(Atom *_pAtom, unsigned int i, double _r0, double _w0, double _a, double _r1, double _w1, double _p1, double _dZ, int _nRes);
		//void generateGevorgCoor(Atom *_pAtom, unsigned int i, double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, double _corr, double _wa, double _zCenter);
		//BBQTable bbqTable;


};
//INLINES
inline AtomPointerVector& CoiledCoils::getAtomPointers() { return atoms; }

}

#endif
