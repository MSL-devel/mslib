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

#ifndef RESIDUEPAIRTABLE_H
#define RESIDUEPAIRTABLE_H

#include "Residue.h"

#include <string>
#include <map>

namespace MSL { 
class ResiduePairTable {
	
	public:
		ResiduePairTable();
		ResiduePairTable(const ResiduePairTable &_rpt);
		~ResiduePairTable();

		void operator=(const ResiduePairTable &_rpt);
		void addResiduePair(std::string _res1, std::string _res2, double _value);
	
		double getValue(std::string _r1, std::string _r2);

		std::map<std::string,double> getPairTable() const { return pairTable; }
	private:

		void copy(const ResiduePairTable &_rpt);
		std::map<std::string,double> pairTable;
};
}

#endif
