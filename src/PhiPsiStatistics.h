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

#ifndef PHIPSISTATISTICS
#define PHIPSISTATISTICS


//MSL Includes
#include "Residue.h"

// STL Includes
using namespace std;

class PhiPsiStatistics {
	public:
		PhiPsiStatistics();
		PhiPsiStatistics(const PhiPsiStatistics &_phiPsiStat);
		~PhiPsiStatistics();

		void operator=(const PhiPsiStatistics &_phiPsiStat);
		void addStatisitics(string _residueType, string _phiBin, string _psiBin,int _count);

		int operator()(string _key);

		int getCounts(Residue nMinus1, Residue n, Residue nPlus1);
		double getProbability(Residue nMinus1, Residue n, Residue nPlus1);
		double getProbabilityAll(Residue nMinus1, Residue n, Residue nPlus1);
		double getPropensity(Residue nMinus1, Residue n, Residue nPlus1);

		static double getPhi(Residue nMinus1, Residue n);
		static double getPsi(Residue n, Residue nPlus1);
		void computeTotalCounts();

		map<string,int>  getPhiPsiCounts() const { return phiPsiTable; }
	private:
		
		void copy(const PhiPsiStatistics &_phiPsiStat);
		double getPhiPsiBin(double _angle);
        double gridSize;

		map<string,int> phiPsiTable;
	
	
	
};
#endif
