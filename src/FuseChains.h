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
#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "Chain.h"
using namespace std;

namespace MSL {
  class FuseChains {
	public:
		FuseChains();
		FuseChains(const FuseChains &_fuse);
		~FuseChains();

		/*
		  Fuse a insert into a template chain, meaning 2 fuse points. Use this for loop insertions.
		  This assumes that the insert chain has been pre-aligned into the proper spot in the template.
		  The function will renumber all the atoms from 1 - n , where n is _template.size()-2 + _insert.size().
		*/
		AtomPointerVector & fuseInsert(Chain &_template, Chain &_insert, string templateStem1posId="", string templateStem2posId="", bool _includeTemplateStems=false);


	private:
		AtomContainer fusedChains;

	
  };
}
