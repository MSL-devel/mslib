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

#ifndef SYSTEMROTAMERLOADER_H
#define SYSTEMROTAMERLOADER_H

#include <iostream>

#include "RotamerLibrary.h"
#include "RotamerLibraryReader.h"
#include "System.h"

using namespace std;

class SystemRotamerLoader {
	public:
		SystemRotamerLoader();
		SystemRotamerLoader(System & _sys);
		SystemRotamerLoader(System & _sys, string _libraryFile);
		SystemRotamerLoader(const SystemRotamerLoader & _sysrotload);
		~SystemRotamerLoader();

		bool readRotamerLibraryFile(string _libraryFile);

		void setSystem(System & _sys);
		void setRotamerLibrary(RotamerLibrary * _pRotlib);

		RotamerLibrary * getRotamerLibrary() const;

		bool loadRotamers(unsigned int _resIndex, string _rotLib, string _residue, int _start, int _end, bool _keepOldRotamers=false);
		bool loadRotamers(string _chainId, string _resNumAndIcode, string _rotLib, string _residue, int _start, int _end, bool _keepOldRotamers=false);
		bool loadRotamers(Position * _pos, string _rotLib, string _residue, int _start, int _end, bool _keepOldRotamers=false);


		// the next functions add rotamers, preserving the old one
		// it is basically wrapper functions to load rotamers with _keepOldRotamers=true
		bool addRotamers(unsigned int _resIndex, string _rotLib, string _residue, int _start, int _end);
		bool addRotamers(string _chainId, string _resNumAndIcode, string _rotLib, string _residue, int _start, int _end);
		bool addRotamers(Position * _pos, string _rotLib, string _residue, int _start, int _end);

	private:
		void setup(System * _pSys, string _libraryFile);
		void deletePointers();

		bool deleteRotLib_flag;
		
		RotamerLibrary * pRotLib;
		RotamerLibraryReader * pRotRead;
		System * pSystem;

};

inline void SystemRotamerLoader::setSystem(System & _sys) {pSystem = &_sys;}
inline void SystemRotamerLoader::setRotamerLibrary(RotamerLibrary * _pRotlib) {if (deleteRotLib_flag) {delete pRotLib;} pRotLib = _pRotlib; deleteRotLib_flag=false;}
inline RotamerLibrary * SystemRotamerLoader::getRotamerLibrary() const {return pRotLib;}
inline bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, string _rotLib, string _residue, int _start, int _end) {
	// add rotamers, preserving the old one
	// it is basically wrapper functions to load rotamers with _keepOldRotamers=true
	return loadRotamers(_resIndex, _rotLib, _residue, _start, _end, true);
}
inline bool SystemRotamerLoader::addRotamers(string _chainId, string _resNumAndIcode, string _rotLib, string _residue, int _start, int _end) {
	return loadRotamers(_chainId, _resNumAndIcode, _rotLib, _residue, _start, _end, true);
}
inline bool SystemRotamerLoader::addRotamers(Position * _pos, string _rotLib, string _residue, int _start, int _end) {
	return loadRotamers(_pos, _rotLib, _residue, _start, _end, true);
}

#endif
