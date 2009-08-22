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

#ifndef CHISTATISTICS
#define CHISTATISTICS


//MSL Includes
#include "Residue.h"

// STL Includes
using namespace std;

class ChiStatistics {
	public:
		ChiStatistics();
		ChiStatistics(const ChiStatistics &_chiStat);
		~ChiStatistics();

		void operator=(const ChiStatistics &_chiStat);
		/*
		void addStatisitics(string _residueType, int _chiNumber, string _chiBin, int _count);

		int operator()(string _key);
		
		int getCounts(Residue &_n, int _chiNumber);
		double getProbability(Residue &_n,int _chiNumber);
		double getProbabilityAll(Residue &_n, int _chiNumber);
		double getPropensity(Residue &_n, int _chiNumber);
		void computeTotalCounts();
		map<string,int>  getChiCounts() const { return chiTable; }
		*/

		bool   atomsExist(Residue &_n, int _chiNumber);
		int    getNumberChis(Residue &_n);
		double getChi(Residue &_n, int _chiNumber,bool _angleInRadians=false);



	private:
		
		void copy(const ChiStatistics &_phiPsiStat);

		/*
		double getChiBin(double _angle);
		map<string,vector<int> >  chiTable;
		*/

		// Stores atom names for each residue, each chi angle
		map<string, vector< vector<string> > > chis;
	
	
};
#endif
