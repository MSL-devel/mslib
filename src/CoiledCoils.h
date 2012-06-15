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

#ifndef COILEDCOILS_H
#define COILEDCOILS_H



// STL Includes


// MSL Includes
#include "System.h"
#include "Symmetry.h"
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

		bool setSystemToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  std::string _symmetry, int _N, System& _sys, std::vector<std::string> _startingPositions);
		// DEPRECATED
		bool primarySequenceToCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  std::string _symmetry, int _N, System& _sys, std::vector<std::string> _startingPositions);

		AtomPointerVector& getCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes);
		AtomPointerVector& getCoiledCoilCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes);
		AtomPointerVector& getCoiledCoilBundle(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes,  std::string _symmetry, int _N);
		AtomPointerVector& getCoiledCoilBundleCricks(double _r0, double _w0, double _a, double _r1, double _w1, double _phi1, double _dZ, int _nRes, std::string _symmetry, int _N);

		AtomPointerVector& getAtomPointers();

		void setBackboneAtomNames(std::string _CAname, std::string _Nname, std::string _Cname, std::string _Oname); // Defaults CA N C O

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
		
		//System *sys;

		void radCoiledCoil(double _r0, double _risePerRes, double _pitch, double _r1, double _w1, double _phi1, double _dZ, int _nRes);

		std::string CAname;
		std::string Nname;
		std::string Cname;
		std::string Oname;

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
inline void CoiledCoils::setBackboneAtomNames(std::string _CAname, std::string _Nname, std::string _Cname, std::string _Oname) {
	// Defaults are CA N C O
	CAname = _CAname;
	Nname = _Nname;
	Cname = _Cname;
	Oname = _Oname;

}

}

#endif
