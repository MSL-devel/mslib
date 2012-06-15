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
#ifndef CHISTATISTICS_H
#define CHISTATISTICS_H


//MSL Includes
#include "Residue.h"
#include "DegreeOfFreedomReader.h"
#include "SysEnv.h"

// STL Includes

namespace MSL { 
class ChiStatistics {



	public:
		ChiStatistics();
		ChiStatistics(const ChiStatistics &_chiStat);

		// read the degreeoffreedom file
		bool read(std::string _dofFile);
		~ChiStatistics();

		void operator=(const ChiStatistics &_chiStat);
		/*
		void addStatisitics(std::string _residueType, int _chiNumber, std::string _chiBin, int _count);

		int operator()(std::string _key);
		
		int getCounts(Residue &_n, int _chiNumber);
		double getProbability(Residue &_n,int _chiNumber);
		double getProbabilityAll(Residue &_n, int _chiNumber);
		double getPropensity(Residue &_n, int _chiNumber);
		void computeTotalCounts();
		std::map<std::string,int>  getChiCounts() const { return chiTable; }
		*/

		bool   atomsExist(Residue &_n, int _chiNumber);
		int    getNumberChis(Residue &_n);
		double getChi(Residue &_n, int _chiNumber,bool _angleInRadians=false);
		std::vector<double>  getChis(Residue &_n,bool _angleInRadians=false);

	private:
		
		void copy(const ChiStatistics &_phiPsiStat);
		DegreeOfFreedomReader dofReader;

		/*
		double getChiBin(double _angle);
		std::map<std::string,std::vector<int> >  chiTable;
		*/

};
}

#endif
