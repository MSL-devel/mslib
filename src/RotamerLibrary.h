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

#ifndef ROTAMERLIBRARY_H
#define ROTAMERLIBRARY_H

#include "MslTools.h"

// STL Includes
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>


namespace MSL { 
class RotamerLibraryReader;
class RotamerLibraryWriter;

class RotamerLibrary {

	public:
		// Constructors/Destructors
		RotamerLibrary();
		RotamerLibrary(const RotamerLibrary & _rotlib);
		~RotamerLibrary();

		/*********************************************************
		 *  Setup functions
		 *********************************************************/
		bool addLibrary(std::string _libName);
		bool addResidue(std::string _libName, std::string _resName);
		bool addInitAtoms(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms); // DEPRECATED
		bool addMobileAtoms(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms);
		bool addInternalCoorDefinition(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms); // Assumes _atoms is always of size 4.  blank std::strings "" added if bond,angle.
		bool addConformation(std::string _libName, std::string _resName, const std::vector<double> & _values, unsigned int _bin = 0); 

		/*********************************************************
		 *  Queries
		 *********************************************************/
		unsigned int getNumberOfLibraries() const;
		bool libraryExists(std::string _libName) const;
		bool residueExists(std::string _libName, std::string _resName);
		unsigned int size(std::string _libName, std::string _resName);
		// backbone dependent version
		unsigned int size(std::string _libName, std::string _resName, double phi, double psi);

		/* I/O */
		bool readFile(std::string _filename, std::string _beblFile="", bool _append=false);
		bool writeFile(std::string _filename, std::string beblFile="");

		/*********************************************************
		 *  Getters
		 *********************************************************/
		struct InternalCoorDefi {
			unsigned int type;
			std::vector<std::string> atomNames;
			//std::vector<int> resnumCorrectors;
		};

		struct RotamerBuildingIC {
			std::vector<std::string> atomNames; // the name of the IC entry that could build the atom
			std::vector<unsigned int> defiIndeces; // the indexes of the entry of the defi that correspond to the values
			bool improper_flag;
		};

		std::vector<std::string> getInitAtoms(std::string _libName, std::string _resName); // DEPRECATED
		std::vector<std::string> getMobileAtoms(std::string _libName, std::string _resName);
		std::vector<InternalCoorDefi> getInternalCoorDefinition(std::string _libName, std::string _resName);
		std::vector<std::vector<double> > getInternalCoor(std::string _libName, std::string _resName);

		std::vector<unsigned int> getRotamerBins(std::string _libName, std::string _resName);
		std::map<std::string, RotamerBuildingIC> getRotamerBuildingIC(std::string _libName, std::string _resName);
		std::string getInitAtomsLine(std::string _libName,std::string _resName); // DEPRECATED
		std::string getMobileAtomsLine(std::string _libName,std::string _resName);
		std::vector<std::string> getInternalCoorDefinitionLines(std::string _libName, std::string _resName) ;
		std::vector<std::string> getAllInternalCoorLines(std::string _libName, std::string _resName) ;
		std::string getInternalCoorLine(std::string _libName, std::string _resName, unsigned int _num) ;

		void setLevel(std::string _levelName, std::string _resName, unsigned int _numRots);
		unsigned int getLevel(std::string _levelName, std::string _resName);
	
		// TODO How about the backbone-dependent version for getAllLevels?
		std::map<std::string, std::map<std::string,unsigned int> > getAllLevels();

		std::string getDefaultLibrary();
		std::set<std::string> getResList (std::string _libName); 
		std::set<std::string>  getAllResList(); 
		std::vector<std::string> getLibraryNames() const;

		bool calculateBuildingICentries();
	//	num is 0-based
		void removeAllConformations();

		// remove the _num'th rotamer from the library
		bool removeRotamer(std::string _libName,std::string _resName,int _num);

		// removes _numOfRotamers starting from _startIdx (if _numOfRotamers = -1 removes till the end)
		bool removeRotamers(std::string _libName,std::string _resName,unsigned _startIdx, unsigned _numOfRotamers = -1);

		bool trimToLevel(std::string _libName,std::string _level);

		void reset();

		std::string toString();

		bool isBackboneDependent() const;
		// Methods for bebl

		void setBeblLevelLabels(std::vector<std::string>& _levelLabels);
		void setPhiBinSize(int _phiBinSize);
		void setPsiBinSize(int _psiBinSize);

		std::vector<std::string>& getBeblLevelLabels();
		int getPhiBinSize();
		int getPsiBinSize();

		bool addBeblBin(std::string _resName, int _phiBin, int _psiBin, std::vector<unsigned> _numConfsPerLevel, std::vector<unsigned> _confIndices);
		std::vector<std::vector<double> > getInternalCoor(std::string _libName, std::string _resName, double phi, double psi);
		unsigned int getLevel(std::string _levelName, std::string _resName, double phi, double psi);
		int getDefaultBin();

	protected:		
	private:
		void setup();
		void copy(const RotamerLibrary & _rotlib);
		void deletePointers();

		std::pair<int, int> getPhiPsiBin(std::string _resName, double _phi, double _psi);
		bool beblResidueExists(std::string resName);

		// A residue, contains the name of the residue, 
		// the degrees of freedom, the conformations
		// the
		//
		// the mobileAtoms std::vector contains all atoms that 
		// need to be initalized
		struct Res {
			std::vector<std::string> mobileAtoms;
			std::vector<InternalCoorDefi> defi; 
			std::vector<std::vector<double> > internalCoor;
			std::vector<unsigned int > rotamerBins; // bin 0 is the default bin .. actual bin numbers start from 1
			std::map<std::string, RotamerBuildingIC> buildingInstructions;
		};

		// the first level std::map std::string is the library name, the second level is the residue name
		std::map<std::string, std::map<std::string, Res> > libraries;  
		std::map<std::string,std::map<std::string,unsigned int> > levels;

		std::string defaultLibrary;
		std::map<std::string, Res>::iterator lastFoundRes;
		RotamerLibraryReader *rotReader;
		RotamerLibraryWriter *rotWriter;

		// Members to store  Bebl info

		bool isBbdep; // is this library backbonedependent
		
		struct BeblBinInfo {
			std::vector<unsigned> numConfsPerLevel; // number of conformers at each level 
			std::vector<unsigned> confIndices; 		
		};	
	
		std::vector<std::string> levelLabels; // defines the correspondence between levellabels and numConfs
		// map from resName to phiBin to psiBin to BeblBinInfo
		std::map<std::string,std::map<int,std::map<int,BeblBinInfo> > > bebl;
		int phiBinSize;
		int psiBinSize;

};

// INLINE FUNCTIONS

inline bool RotamerLibrary::beblResidueExists(std::string resName) {return (bebl.find(resName) != bebl.end());}
inline int RotamerLibrary::getDefaultBin() {return MslTools::intMax;}
inline unsigned int RotamerLibrary::getNumberOfLibraries() const {return libraries.size();}
inline std::vector<std::string> RotamerLibrary::getLibraryNames() const {
	std::vector<std::string> libraryNames;
	for (std::map<std::string, std::map<std::string, Res> >::const_iterator l = libraries.begin(); l!= libraries.end(); l++) {
		libraryNames.push_back(l->first);
	}
	return libraryNames;
}
inline bool RotamerLibrary::addLibrary(std::string _libName) {if (libraries.size() == 0) {defaultLibrary=_libName;} libraries[_libName]; return true;}
inline bool RotamerLibrary::addResidue(std::string _libName, std::string _resName) {if (_libName == "") {_libName = defaultLibrary;} if (libraryExists(_libName)) {libraries[_libName][_resName]; return true;} else {return false;}}
inline bool RotamerLibrary::addInitAtoms(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms) {
	// DEPRECATED 
	std::cerr << "WARNING: using deprecated function inline bool RotamerLibrary::addInitAtoms(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms), use addMobileAtoms instead" << std::endl;
	return addMobileAtoms(_libName, _resName, _atoms);
}
inline bool RotamerLibrary::addMobileAtoms(std::string _libName, std::string _resName, const std::vector<std::string> & _atoms) {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	if (residueExists(_libName, _resName)) {
		libraries[_libName][_resName].mobileAtoms = _atoms;
		return true;
	} else {
		return false;
	}
}

inline bool RotamerLibrary::addConformation(std::string _libName, std::string _resName, const std::vector<double> & _values,unsigned int _bin) {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	if (residueExists(_libName, _resName)) {
		if (_values.size() != libraries[_libName][_resName].defi.size()) {
			// the number of numbers must match the number of internal coor definitions
			return false;
		}

		libraries[_libName][_resName].internalCoor.push_back(_values);
		libraries[_libName][_resName].rotamerBins.push_back(_bin);
		return true;
	} else {
		return false;
	}
}
inline bool RotamerLibrary::libraryExists(std::string _libName) const {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	return libraries.find(_libName) != libraries.end();
}
inline bool RotamerLibrary::residueExists(std::string _libName, std::string _resName) {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	if (libraryExists(_libName)) {
		//return libraries[_libName].find(_resName) != libraries[_libName].end();
		lastFoundRes = libraries.find(_libName)->second.find(_resName);
		return lastFoundRes != libraries.find(_libName)->second.end();
		//return libraries.find(_libName)->second.find(_resName) != libraries.find(_libName)->second.end();
	} else {
		return false;
	}
}
inline std::vector<std::string> RotamerLibrary::getInitAtoms(std::string _libName, std::string _resName) {
	// DEPRECATED 
	std::cerr << "WARNING: using deprecated function inline std::vector<std::string> RotamerLibrary::getInitAtoms(std::string _libName, std::string _resName), use getMobileAtoms instead" << std::endl;
	return getMobileAtoms(_libName, _resName);
}
inline std::vector<std::string> RotamerLibrary::getMobileAtoms(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.mobileAtoms;
	}
	return std::vector<std::string>();

}
inline std::vector<RotamerLibrary::InternalCoorDefi> RotamerLibrary::getInternalCoorDefinition(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.defi;
	} else {
		return std::vector<InternalCoorDefi>();
	}
}
inline std::vector<unsigned int> RotamerLibrary::getRotamerBins(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.rotamerBins;
	} else {
		return std::vector<unsigned int>();
	}
}
/*
inline void RotamerLibrary::getInternalCoorDefinition(std::string _libName, std::string _resName, std::vector<unsigned int> & _type, std::vector<std::vector<std::string> > & _atomNames, std::vector<std::vector<int> > & _resnumCorrectors) {
	_type.clear();
	_atomNames.clear();
	_resnumCorrectors.clear();
	if (residueExists(_libName, _resName)) {
		for (std::vector<InternalCoorDefi>::iterator k=lastFoundRes->second.defi.begin(); k!=lastFoundRes->second.defi.end(); k++) {
			_type.push_back(k->type);
			_atomNames.push_back(k->atomNames);
			_resnumCorrectors.push_back(k->resnumCorrectors);
		}
	}
}
*/
inline std::map<std::string, RotamerLibrary::RotamerBuildingIC> RotamerLibrary::getRotamerBuildingIC(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.buildingInstructions;
	} else {
		return std::map<std::string, RotamerBuildingIC>();
	}
}
inline std::vector<std::vector<double> > RotamerLibrary::getInternalCoor(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.internalCoor;
	}
	return std::vector<std::vector<double> >();
}

inline std::string RotamerLibrary::getDefaultLibrary() {
	return defaultLibrary;
}
inline std::set<std::string> RotamerLibrary::getResList(std::string _libName) {

	if(libraryExists(_libName)) {
		std::set<std::string> list;
		for(std::map<std::string,Res>::iterator i = libraries[_libName].begin(); i != libraries[_libName].end(); i++) {
			list.insert(i->first);
			//std::cout << "UUUUU inserting " << i->first << std::endl;
		}
		return list;
	} else {
		return std::set<std::string> ();
	}
}


inline std::set<std::string> RotamerLibrary::getAllResList() {
	std::set<std::string> list;
	 //std::cout << "UUUUU resList.size" << resList.size() << std::endl;
	for (std::map <std::string, std::map<std::string,Res> >::iterator i = libraries.begin(); i != libraries.end(); i++) {
		 // std::cout << "UUUUU size of std::vector" << (*i).size() << std::endl;
		for(std::map<std::string,Res>::iterator j = i->second.begin(); j != i->second.end(); j++) {
				
			list.insert(j->first);
		}
	}
	return list;
}



inline void RotamerLibrary::setLevel(std::string _levelName, std::string _resName, unsigned int _numRots) {
	levels[_levelName][_resName] = _numRots;
}
inline unsigned int RotamerLibrary::getLevel(std::string _levelName, std::string _resName) {
	if(levels.find(_levelName) != levels.end()) {
		if(levels[_levelName].find(_resName) != levels[_levelName].end()) {
			return levels[_levelName][_resName];
		} else {
			std::cerr << "WARNING 10324: LEVEL " << _levelName << " not defined for " << _resName << std::endl;
			return 0;
		}
	} else {
		std::cerr << "WARNING 10324: LEVEL " << _levelName << " does not exist." << std::endl;
		return 0;
	}
}

inline std::map<std::string,std::map<std::string,unsigned int> > RotamerLibrary::getAllLevels() {
	return levels;
}

inline unsigned int RotamerLibrary::size(std::string _libName, std::string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.internalCoor.size();
	} else {
		return 0;
	}
}

inline bool RotamerLibrary::isBackboneDependent() const {return isBbdep;} 
inline void RotamerLibrary::setBeblLevelLabels(std::vector<std::string>& _levelLabels) {levelLabels = _levelLabels;}
inline void RotamerLibrary::setPhiBinSize(int _phiBinSize) {phiBinSize = _phiBinSize;}
inline void RotamerLibrary::setPsiBinSize(int _psiBinSize) {psiBinSize = _psiBinSize;}
inline std::vector<std::string>& RotamerLibrary::getBeblLevelLabels() {return levelLabels;}
inline int RotamerLibrary::getPhiBinSize() {return phiBinSize;}
inline int RotamerLibrary::getPsiBinSize() {return psiBinSize;}
}

#endif
