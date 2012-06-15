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

#ifndef SYSTEMROTAMERLOADER_H
#define SYSTEMROTAMERLOADER_H

#include <iostream>

#include "RotamerLibrary.h"
#include "RotamerLibraryReader.h"
#include "System.h"


namespace MSL { 
class SystemRotamerLoader {
	public:
		SystemRotamerLoader();
		SystemRotamerLoader(System & _sys);
		SystemRotamerLoader(System & _sys, std::string _libraryFile);
		SystemRotamerLoader(const SystemRotamerLoader & _sysrotload);
		~SystemRotamerLoader();

		bool readRotamerLibraryFile(std::string _libraryFile);

		void setSystem(System & _sys);
		void setRotamerLibrary(RotamerLibrary * _pRotlib);

		bool defineRotamerSamplingLevels();

		RotamerLibrary * getRotamerLibrary() const;

		// load rotamers using levellabels
		bool loadRotamers(unsigned int _resIndex, std::string _residue, std::string _level, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(std::string _positionId, std::string _residue, std::string _level, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(Position * _pos, std::string _residue, std::string _level, std::string _rotLib="", bool _keepOldRotamers=false);

		// load n rotamers
		bool loadRotamers(unsigned int _resIndex, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(std::string _positionId, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(Position * _pos, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="", bool _keepOldRotamers=false);
		// load a specific range of rotamers
		bool loadRotamers(unsigned int _resIndex, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(std::string _positionId, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="", bool _keepOldRotamers=false);
		bool loadRotamers(Position * _pos, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="", bool _keepOldRotamers=false);
		// DEPRECATED FUNCTIONS
		bool loadRotamers(unsigned int _resIndex, std::string _rotLib, std::string _residue, int _start, int _end, bool _keepOldRotamers=false); // DEPRECATED
		bool loadRotamers(std::string _positionId, std::string _rotLib, std::string _residue, int _start, int _end, bool _keepOldRotamers=false); // DEPERCATED
		bool loadRotamers(Position * _pos, std::string _rotLib, std::string _residue, int _start, int _end, bool _keepOldRotamers=false); // DEPRECATED
		bool loadRotamers(std::string _chainId, std::string _resNumAndIcode, std::string _rotLib, std::string _residue, int _start, int _end, bool _keepOldRotamers=false);


		// the next functions add rotamers, preserving the old one
		// it is basically wrapper functions to load rotamers with _keepOldRotamers=true
		bool addRotamers(unsigned int _resIndex, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="");
		bool addRotamers(std::string _positionId, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="");
		bool addRotamers(Position * _pos, std::string _residue, unsigned int _numberOfRots, std::string _rotLib="");
		bool addRotamers(unsigned int _resIndex, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="");
		bool addRotamers(std::string _positionId, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="");
		bool addRotamers(Position * _pos, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib="");
		//DEPRECATED
		bool addRotamers(unsigned int _resIndex, std::string _rotLib, std::string _residue, int _start, int _end);
		bool addRotamers(std::string _positionId, std::string _rotLib, std::string _residue, int _start, int _end);
		bool addRotamers(Position * _pos, std::string _rotLib, std::string _residue, int _start, int _end);
		bool addRotamers(std::string _chainId, std::string _resNumAndIcode, std::string _rotLib, std::string _residue, int _start, int _end);


		// getter
		std::string getRotamerLibraryFileName() const;
		
	private:
		void setup(System * _pSys, std::string _libraryFile);
		void deletePointers();
	
		bool deleteRotLib_flag;
		
		RotamerLibrary * pRotLib;
		RotamerLibraryReader * pRotRead;

		std::string rotamerLibraryFile;

		System * pSystem;

};
inline std::string SystemRotamerLoader::getRotamerLibraryFileName() const { return rotamerLibraryFile;}
inline void SystemRotamerLoader::setSystem(System & _sys) {pSystem = &_sys;}
inline void SystemRotamerLoader::setRotamerLibrary(RotamerLibrary * _pRotlib) {if (deleteRotLib_flag) {delete pRotLib;} pRotLib = _pRotlib; deleteRotLib_flag=false;}
inline bool SystemRotamerLoader::defineRotamerSamplingLevels() {
	if(pSystem) {
		if(pRotLib) {
			std::map<std::string,std::map<std::string,unsigned int> > levels = pRotLib->getAllLevels();
			return pSystem->defineRotamerSamplingLevels(levels);
		} else {
			std::cerr << "ERROR 12466: RotamerLibrary not set in SystemRotamerLoader::defineRotamerSamplingLevels()" << std::endl;
			return false;
		}
	} else {
		std::cerr << "ERROR 12466: System not set in SystemRotamerLoader::defineRotamerSamplingLevels()" << std::endl;
		return false;
	}
}
inline RotamerLibrary * SystemRotamerLoader::getRotamerLibrary() const {return pRotLib;}

// the following add rotamers, preserving the old one
// they are wrapper functions to load rotamers with _keepOldRotamers=true
inline bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, std::string _residue, unsigned int _numberOfRots, std::string _rotLib) {
	return loadRotamers(_resIndex, _residue, _numberOfRots, _rotLib, true);
}
inline bool SystemRotamerLoader::addRotamers(std::string _positionId, std::string _residue, unsigned int _numberOfRots, std::string _rotLib) {
	return loadRotamers(_positionId, _residue, _numberOfRots, _rotLib, true);
}
inline bool SystemRotamerLoader::addRotamers(Position * _pos, std::string _residue, unsigned int _numberOfRots, std::string _rotLib) {
	return loadRotamers(_pos, _residue, _numberOfRots, _rotLib, true);
}

inline bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib) {
	return loadRotamers(_resIndex, _residue, _start, _end, _rotLib, true);
}
inline bool SystemRotamerLoader::addRotamers(std::string _positionId, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib) {
	return loadRotamers(_positionId, _residue, _start, _end, _rotLib, true);
}
inline bool SystemRotamerLoader::addRotamers(Position * _pos, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib) {
	return loadRotamers(_pos, _residue, _start, _end, _rotLib, true);
}

inline bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, std::string _rotLib, std::string _residue, int _start, int _end) {
	std::cerr << "DEPRECATED bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, std::string _rotLib, std::string _residue, int _start, int _end), use bool SystemRotamerLoader::addRotamers(unsigned int _resIndex, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib) instead" << std::endl; 
	return loadRotamers(_resIndex, _rotLib, _residue, _start, _end, true);
}
inline bool SystemRotamerLoader::addRotamers(std::string _positionId, std::string _rotLib, std::string _residue, int _start, int _end) {
	std::cerr << "DEPRECATED bool SystemRotamerLoader::addRotamers(std::string _positionId, std::string _rotLib, std::string _residue, int _start, int _end), use bool SystemRotamerLoader::addRotamers(std::string _positionId, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib) instead" << std::endl; 
	return loadRotamers(_positionId, _rotLib, _residue, _start, _end, true);
}
inline bool SystemRotamerLoader::addRotamers(Position * _pos, std::string _rotLib, std::string _residue, int _start, int _end) {
	std::cerr << "DEPRECATED bool SystemRotamerLoader::addRotamers(Position * _pos, std::string _rotLib, std::string _residue, int _start, int _end), use bool SystemRotamerLoader::addRotamers(Position * _pos, std::string _residue, unsigned int _start, unsigned int _end, std::string _rotLib instead" << std::endl; 
	return loadRotamers(_pos, _rotLib, _residue, _start, _end, true);
}

}

#endif
