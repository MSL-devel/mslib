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



#ifndef EZPOTENTIALBUILDER_H
#define EZPOTENTIALBUILDER_H

#include <iostream>
#include <vector>
#include <map>

#include "System.h"
#include "PolymerSequence.h"
#include "EZpotentialInteraction.h"

/*******************************************************
 *  This is an implementation of the EZ empirical membrane 
 *  insertion potential.  If used please cite:
 *
 *     Senes A, Chadi DC, Law PB, Walters RF, Nanda V, DeGrado WF. 
 *     "E(z), a depth-dependent potential for assessing the energies 
 *     of insertion of amino acid side-chains into membranes, 
 *     derivation and applications to determining the orientation 
 *     of transmembrane and interfacial helices."
 *     J Mol Biol. 2007 366(2), 436-48
 *
 *******************************************************/

namespace MSL { 
	class EZpotentialBuilder {
		public:
			EZpotentialBuilder();
			EZpotentialBuilder(System & _system);
			EZpotentialBuilder( EZpotentialBuilder & _sysBuild);
			~EZpotentialBuilder();

			void operator=( EZpotentialBuilder & _sysBuild);

			void setSystem(System & _system);

			//bool readParameters(std::string _parameterFile);

			bool buildInteractions();

			void setUseCB(bool _flag);
			bool getUseCB() const;

			void setAddTermini(bool _flag);
			bool getAddTermini() const;

			//void printParameters();

			
		private:
			void setup();
			void copy(EZpotentialBuilder & _sysBuild);
			void deletePointers();
			void reset();
			void setParams();

			System * pSystem;

			struct ResidueParams {
				std::vector<double> params;
				bool sigmoidalFunc;
				bool Nterminal;
				bool Cterminal;
			};

			map<string, map<std::string, ResidueParams> > parameters;

			bool useCB_flag;
			bool addTermini_flag;

	};
	inline void EZpotentialBuilder::setSystem(System & _system) {
		pSystem = &_system;
	}

inline void EZpotentialBuilder::setUseCB(bool _flag) {
	useCB_flag = _flag;
}
inline bool EZpotentialBuilder::getUseCB() const {
	return useCB_flag;
}
inline void EZpotentialBuilder::setAddTermini(bool _flag) {
	addTermini_flag = _flag;
}
inline bool EZpotentialBuilder::getAddTermini() const {
	return addTermini_flag;
}

}

#endif
