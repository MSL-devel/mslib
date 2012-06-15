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


#ifndef HYDROGENBONDBUILDER_H
#define HYDROGENBONDBUILDER_H

#include <iostream>
#include <vector>
#include <map>

#include "System.h"
#include "Reader.h"
#include "PolymerSequence.h"
#include "Scwrl4HBondInteraction.h"

namespace MSL { 
	class HydrogenBondBuilder {
		public:
			HydrogenBondBuilder();
			HydrogenBondBuilder(System & _system, std::string _scwrl4ParameterFile);
			HydrogenBondBuilder( HydrogenBondBuilder & _sysBuild);
			~HydrogenBondBuilder();

			void operator=( HydrogenBondBuilder & _sysBuild);

			void setSystem(System & _system);

			bool readParameters(std::string _scwrl4ParameterFile);

			bool buildInteractions(double _cutoff = -1.0); // called the first time interactions are built, a negative cutoff means no distance cutoff will be applied while creating interactions
			bool update(double _cutoff = -1.0); // works on the donor and acceptor list created by buildInteractions

			void printParameters();

			/**************************************************
			 * Add one or more new identities to a position,
			 * if a series of backbone atoms are specified
			 * (_bbAtoms) then the coordinates are passed to
			 * the atoms with the same name in the new residues
			 * to preserve the position of the backbone
			 **************************************************/
			
		private:
			void setup();
			void copy(HydrogenBondBuilder & _sysBuild);
			void collectDonorsAndAcceptors();
			void deletePointers();
			void reset();

			System * pSystem;

			std::vector<std::vector<Atom*> > acceptors; // 0 - acceptor(O), 1 - lonePairAngleAtom(C), 2 - lonePairDihedralAtom  (CA)
			std::vector<std::vector<Atom*> > donors; // 0- hydrogen(HN), 1 - donorAtom (N)

			std::map<std::string,std::string> donorAtom; // [ARG HN] = N
			std::map<std::string,std::string> lonePairAngleAtom; // [ASP O] = C
			std::map<std::string,std::string> lonePairDihedralAtom; // [ASP O] =[CA]

			std::map<std::string,std::vector<double> > donorData; // [ARG HN ] [2.08 0.67 35 37 49 ] 

			std::map<std::string,std::vector<double> > acceptorData; // [ASP O ] [1 120 0 180]

	};
	inline void HydrogenBondBuilder::setSystem(System & _system) {
		reset();
		pSystem = &_system;
	}
}

#endif
