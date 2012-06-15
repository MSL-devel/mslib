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

#ifndef OPTIMAL_RMSD_CALCULATOR_H_
#define OPTIMAL_RMSD_CALCULATOR_H_

#include <vector>
#include "AtomPointerVector.h"
using namespace std;

namespace MSL {

class OptimalRMSDCalculator {

 public:
  
	OptimalRMSDCalculator() {}
	~OptimalRMSDCalculator() {}

	// getters
	double lastRMSD() { return rmsd; }
	vector<double> lastTranslation();
	vector<vector<double> > lastRotation();

	// calculate optimal superposition and the resulting RMSD, applying transformation to given atoms
	bool align(AtomPointerVector &_align, AtomPointerVector &_ref, AtomPointerVector& _moveable);
	// quickly calculate RMSD upon optimal superposition without generating the rotation matrix
	double bestRMSD(AtomPointerVector &_align, AtomPointerVector &_ref, bool* _suc = NULL);

 protected:
	// implemetation of Kabsch algoritm for optimal superposition
	bool Kabsch(AtomPointerVector &_align, AtomPointerVector &_ref, int mode);

 private:
	double rmsd;
	double t[3];    // translation vector
	double u[3][3]; // rotation matrix

};

}

#endif // OPTIMAL_RMSD_CALCULATOR_H_

