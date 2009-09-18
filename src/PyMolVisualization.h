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

#ifndef PYMOLVISUALIZATION_H
#define PYMOLVISUALIZATION_H

/*
  PyMolVisualization Object
  
  Goals:
    Represent (named) Spheres at given CartesianPoints
    Represent (named) cylinders between 2 CartesianPoints (width optional)
    Represent (named) frame (set of 3 cylinders) from a Frame object
    Represent (named) plane object from 3 CartesianPoints

    Write out any representations.

 */

#include <iostream>
#include <map>

//#include "Grid.h"
#include "CartesianPoint.h"
#include "Frame.h"

class PyMolVisualization {

	public:
		PyMolVisualization();
		~PyMolVisualization();


		bool createAtom(CartesianPoint &_a, string _name="", double _radiusScaling=1.0);               // create a single object for a single atom
		bool createAtom(Atom &_a, string _name="", double _radiusScaling=1.0);               // create a single object for a single atom
		//bool createAtoms(AtomVector &_av, string _name="", double _radiusScaling=1.0);      // create a single object for multiple atoms.
		bool createSphere(CartesianPoint &_cp, string _name="", double _radius=2.0,int rgbR=1,int rgbG=0,int rgbB=0); 
		bool createCylinder(CartesianPoint &_start,CartesianPoint &_end, string _name="",double _radius=2.0,int rgbR=1,int rgbG=0,int rgbB=0);
		bool createCone(CartesianPoint &_start,CartesianPoint &_end, string _name="",double _radius=2.0,int rgbR=1,int rgbG=0,int rgbB=0);
		//bool createFrame(Frame &_frame, string _name="");
		//bool createPlane(CartesianPoint &_cp1, CartesianPoint &_cp2, CartesianPoint &_cp3); 
		
		bool createArrow(CartesianPoint &_start, CartesianPoint &_vector, string _name,double _cylinderRadius=0.1, double _coneRadius=0.2, int rgbR=0,int rgbG=0,int rgbB=1); 


		// Create a openDX style grid
		//bool createGrid(Grid &_grid, string _name="");

		
		string toString();
		friend ostream & operator<<(ostream &_os, PyMolVisualization & _pyObj)  {_os << _pyObj.toString(); return _os;};


		bool   existPyMolObject(string _name);
		string getPyMolObject(string _name);
		


	private:
		map<string,string> pymolObjectStrings;
		map<string,string> pymolAtomStrings;



};
#endif
