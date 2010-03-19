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

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "AtomPointerVector.h"
#include "Transforms.h"

namespace MSL { 
class Symmetry {
	public:
		Symmetry();
		~Symmetry();


		// C2 Symmetry
		//void applyC2(AtomPointerVector &_ats);

		// C3 Symmetry
		//void applyC3(AtomPointerVector &_ats);

		// General CN Symmetry
		void applyCN(AtomPointerVector &_ats, int _N); // N is number of symmetry mates

		// General DN Symmetry
		void applyDN(AtomPointerVector &_ats, int _N); // N is number of symmetry mates


		void applyD2(AtomPointerVector &_ats);
		

		AtomPointerVector& getAtomPointers();
	private:

		AtomPointerVector atoms;
	
};
inline AtomPointerVector& Symmetry::getAtomPointers() { return atoms; }

}

#endif
