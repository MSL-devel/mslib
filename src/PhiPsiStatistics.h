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

namespace MSL { 
class PhiPsiStatistics {
	public:
		PhiPsiStatistics();
		PhiPsiStatistics(const PhiPsiStatistics &_phiPsiStat);
		~PhiPsiStatistics();

		void operator=(const PhiPsiStatistics &_phiPsiStat);
		void addStatisitics(std::string _residueType, std::string _phiBin, std::string _psiBin,int _count);

		int operator()(std::string _key);

		int getCounts(std::string &resName, double phi, double psi);
		int getCounts(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getProbability(std::string &resName, double phi, double psi);
		double getProbability(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getProbabilityAll(double phi, double psi);
		double getProbabilityAll(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);
		double getPropensity(std::string &resName, double phi, double psi);
		double getPropensity(const Residue &nMinus1, const Residue &n, const Residue &nPlus1);

		static double getPhi(const Residue &nMinus1, const Residue &n);
		static double getPsi(const Residue &n, const Residue &nPlus1);
		void computeTotalCounts();

		std::map<std::string,int>  getPhiPsiCounts() const { return phiPsiTable; }
	private:
		
		void copy(const PhiPsiStatistics &_phiPsiStat);
		double getPhiPsiBin(double _angle);
        double gridSize;

		std::map<std::string,int> phiPsiTable;
	
	
	
};
}

#endif
