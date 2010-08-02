/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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


#ifndef HYDROGENBONDBUILDER_H
#define HYDROGENBONDBUILDER_H

#include <iostream>
#include <vector>

#include "System.h"
#include "Reader.h"
#include "PolymerSequence.h"
#include "Scrwl4HBondInteraction.h"

namespace MSL { 
	class HydrogenBondBuilder {
		public:
			HydrogenBondBuilder();
			HydrogenBondBuilder(System & _system, std::string _scrwl4ParameterFile);
			HydrogenBondBuilder( HydrogenBondBuilder & _sysBuild);
			~HydrogenBondBuilder();

			void operator=( HydrogenBondBuilder & _sysBuild);

			void setSystem(System & _system);

			bool readParameters(std::string _scrwl4ParameterFile);

			bool buildInteractions(double _cutoff = -1.0);

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
			void deletePointers();
			void reset();

			System * pSystem;

			std::map<std::string,std::map<std::string,std::string> > donorData; // [ARG][HH11][NH1] [ARG][HH12][NH1] ..

			std::map<std::string,std::map<std::string,std::vector<std::string> > > acceptorData; // [ASP] [O ] [C CA 1 120 0 180]

	};
	inline void HydrogenBondBuilder::setSystem(System & _system) {
		reset();
		pSystem = &_system;
	}
}

#endif
