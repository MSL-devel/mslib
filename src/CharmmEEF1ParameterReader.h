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


#ifndef CHARMMEFF1PARAMETERREADER_H
#define CHARMMEFF1PARAMETERREADER_H
// STL Includes
#include <vector>
#include <map>
#include <iostream>
#include <math.h>

//MSL Includes
#include "Reader.h"
#include "MslTools.h"


namespace MSL { 

	class CharmmEEF1ParameterReader : public Reader {

		public:
			CharmmEEF1ParameterReader();
			CharmmEEF1ParameterReader(const std::string & _filename);
			CharmmEEF1ParameterReader(const CharmmEEF1ParameterReader & _par);
			~CharmmEEF1ParameterReader();

			void operator=(const CharmmEEF1ParameterReader & _par);

			bool read();
			void reset();

			bool solventExists(std::string _solvent) const;

			//VdwParam will return a std::vector with 8 values _V_i, Gfree_i, Sigw_i, rmin_i, V_j, Gfree_j, Sigw_j, rmin_j
			bool EEF1Param(std::vector<double> & _params, std::string type, std::string solvent="") const;
			//To be implemented
			//std::vector<double> IMM1Param(std::string _type, std::string _solvent1, std::string _solvent2) const;
			std::string getDefaultSolvent() const;

		private:
			void addEEF1(std::string _solvent, std::string _atomType, double _V, double _Gref, double _Gfree, double _Href, double _CPref, double _Sigw);
			
			//void setup();
			void copy(const CharmmEEF1ParameterReader & _par);
			
			//EEF1Params will contain a std::vector with 4 values (Eps,Rmin,Esp14,Rmin14)
			// ["solvent"]["atomType"][0] = V    
			// ["solvent"]["atomType"][1] = Gref 
			// ["solvent"]["atomType"][2] = Gfree
			// ["solvent"]["atomType"][3] = Href 
			// ["solvent"]["atomType"][4] = CPref
			// ["solvent"]["atomType"][5] = Sigw 
			std::map<std::string, std::map<std::string,std::vector<double> > > EEF1Map;
		//	std::map<std::string, std::map<std::string,std::map<std::string, std::vector<double> > > EEF1PairMap;
			static const std::string defaultSolvent;


	};

	inline bool CharmmEEF1ParameterReader::solventExists(std::string _solvent) const {
		if (_solvent == "") {
			_solvent = defaultSolvent;
		}
		std::map<std::string, std::map<std::string, std::vector<double> > >::const_iterator found1;
		return EEF1Map.find(_solvent) != EEF1Map.end();
	}
	inline std::string CharmmEEF1ParameterReader::getDefaultSolvent() const {return defaultSolvent;}

}

#endif

