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

#ifndef COILEDCOILS_H
#define COILEDCOILS_H



// STL Includes


// MSL Includes
#include "CartesianPoint.h"
#include "AtomPointerVector.h"
#include "BBQTableReader.h"
#include "BBQTableWriter.h"
#include "BBQTable.h"
#include "System.h"
#include "Atom.h"

// BOOST Includes

/*******************************************************************
 *   NOTE:  REQUIRES COPY CONSTRUCTOR AND = OPERATOR!!!
 *
 *   It will be difficult until we have a system copy constructor.
 *   Can we generate a sys on demand or do we need the getSystem()
 *   call (same goes for Chain)
 *******************************************************************/

class CoiledCoils {

	public:
		CoiledCoils();
		~CoiledCoils();

		void cricksCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1);
		void offersCoiledCoil(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _w1);
		void northCoiledCoils(double _r0, double _risePerRes, double _p0, double _r1, int _nRes, double _theta, double _helicalPhase);
		void applyCoiledCoil(AtomPointerVector &_av, double _p0);
		void sotoCoiledCoils(double _r0, double _risePerRes, double _r1, int _nRes, double _resPerTurn, double _alpha, double _helicalPhase);
				     //int NumberResidues, double Rsuper, double Ralpha, double ResiduesPerTurnMinorHelix, double Alpha, double RisePerResidue, double Phase){

		AtomPointerVector& getAtoms();
		Chain& getChain(); // include all backbone atoms
		System& getSystem(); // include all backbone atoms

		void useBBQTable(string _bbqTable);

		/*
		  Varing parameters as a function of residue number
		  
		  Need to specify:
		     1. Starting residue
		     2. Ending Residue
		     3. Start Params
		     4. End Params
		     5. Param Function
		 */

		void addNorthParamSchedule(vector<int> _residues, vector<double> _t, vector<double> _r1, vector<double> _wAlpha, vector<double> _helicalPhase, vector<double> _r0, vector<double> _p0, vector<int> _paramFunc);

	private:
		
		// North Parameter Schedule Variables
		vector<int> npStartStopResidues;
		vector<double> npPhase;
		vector<double> npR1;
		vector<double> npWAlpha;
		vector<double> npHelicalPhase;
		vector<double> npR0;
		vector<double> npP0;
		vector<int> npParamFuncs;
		
		

		CartesianPoint generateNorthCoor(double _t, double _r1, double _wAlpha, double _helicalPhase, double _r0, double _p0);
		AtomPointerVector atoms;
		BBQTable bbqTable;
		System *sys;


};
//INLINES
inline AtomPointerVector& CoiledCoils::getAtoms() { return atoms; }

#endif
