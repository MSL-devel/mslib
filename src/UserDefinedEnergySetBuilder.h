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

#ifndef USERDEFINEDENERGYSETBUILDER_H
#define USERDEFINEDENERGYSETBUILDER_H

#include "System.h"


/*
  This class takes the System and the Energy set. The system atoms are looped over and 
  the singleton object UserDefinedEnergy is used to build/add UserDefinedInteractions into the energy set.
 */
class UserDefinedEnergySetBuilder {
	public:
		UserDefinedEnergySetBuilder();

		void addEnergyInteractions(System &_sys);
		void addTwoBodyInteractions(System &_sys, string _filename);
	private:

		void setup();

};
#endif
