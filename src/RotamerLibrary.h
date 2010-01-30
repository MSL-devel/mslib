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

#ifndef ROTAMERLIBRARY_H
#define ROTAMERLIBRARY_H

// STL Includes
#include <vector>
#include <map>
#include <set>
#include <iostream>
using namespace std;

class RotamerLibrary {

	public:
		// Constructors/Destructors
		RotamerLibrary();
		RotamerLibrary(const RotamerLibrary & _rotlib);
		~RotamerLibrary();

		/*********************************************************
		 *  Setup functions
		 *********************************************************/
		bool addLibrary(string _libName);
		bool addResidue(string _libName, string _resName);
		bool addInitAtoms(string _libName, string _resName, const vector<string> & _atoms);
		bool addInternalCoorDefinition(string _libName, string _resName, const vector<string> & _atoms); // Assumes _atoms is always of size 4.  blank strings "" added if bond,angle.
		bool addConformation(string _libName, string _resName, const vector<double> & _values); 

		/*********************************************************
		 *  Queries
		 *********************************************************/
		unsigned int getNumberOfLibraries() const;
		bool libraryExists(string _libName) const;
		bool residueExists(string _libName, string _resName);
		unsigned int size(string _libName, string _resName);

		/*********************************************************
		 *  Getters
		 *********************************************************/
		struct InternalCoorDefi {
			unsigned int type;
			vector<string> atomNames;
			//vector<int> resnumCorrectors;
		};

		struct RotamerBuildingIC {
			vector<string> atomNames; // the name of the IC entry that could build the atom
			vector<unsigned int> defiIndeces; // the indexes of the entry of the defi that correspond to the values
			bool improper_flag;
		};

		vector<string> getInitAtoms(string _libName, string _resName);
		vector<InternalCoorDefi> getInternalCoorDefinition(string _libName, string _resName);
		vector<vector<double> > getInternalCoor(string _libName, string _resName);
		map<string, RotamerBuildingIC> getRotamerBuildingIC(string _libName, string _resName);
		string getInitAtomsLine(string _libName,string _resName);
		vector<string> getInternalCoorDefinitionLines(string _libName, string _resName) ;
		vector<string> getAllInternalCoorLines(string _libName, string _resName) ;
		string getInternalCoorLine(string _libName, string _resName, unsigned int _num) ;

		string getDefaultLibrary();
		set<string> getResList (string _libName); 
		set<string>  getAllResList(); 
		vector<string> getLibraryNames() const;

		bool calculateBuildingICentries();
	//	num is 0-based
		void removeAllConformations();
		bool removeRotamer(string _libName,string _resName,int _num);
		void reset();

		string toString();

	protected:		
	private:
		void setup();
		void copy(const RotamerLibrary & _rotlib);

		// A residue, contains the name of the residue, 
		// the degrees of freedom, the conformations
		// the
		//
		// the initAtoms vector contains all atoms that 
		// need to be initalized
		struct Res {
			vector<string> initAtoms;
			vector<InternalCoorDefi> defi; 
			vector<vector<double> > internalCoor;
			map<string, RotamerBuildingIC> buildingInstructions;
		};

		// the first level map string is the library name, the second level is the residue name
		map<string, map<string, Res> > libraries;  

		string defaultLibrary;
		map<string, Res>::iterator lastFoundRes;

};

// INLINE FUNCTIONS

inline unsigned int RotamerLibrary::getNumberOfLibraries() const {return libraries.size();}
inline vector<string> RotamerLibrary::getLibraryNames() const {
	vector<string> libraryNames;
	for (map<string, map<string, Res> >::const_iterator l = libraries.begin(); l!= libraries.end(); l++) {
		libraryNames.push_back(l->first);
	}
	return libraryNames;
}
inline bool RotamerLibrary::addLibrary(string _libName) {if (libraries.size() == 0) {defaultLibrary=_libName;} libraries[_libName]; return true;}
inline bool RotamerLibrary::addResidue(string _libName, string _resName) {if (_libName == "") {_libName = defaultLibrary;} if (libraryExists(_libName)) {libraries[_libName][_resName]; return true;} else {return false;}}
inline bool RotamerLibrary::addInitAtoms(string _libName, string _resName, const vector<string> & _atoms) {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	if (residueExists(_libName, _resName)) {
		libraries[_libName][_resName].initAtoms = _atoms;
		return true;
	} else {
		return false;
	}
}
inline bool RotamerLibrary::addConformation(string _libName, string _resName, const vector<double> & _values) {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	if (residueExists(_libName, _resName)) {
		if (_values.size() != libraries[_libName][_resName].defi.size()) {
			// the number of numbers must match the number of internal coor definitions
			return false;
		}

		libraries[_libName][_resName].internalCoor.push_back(_values);
		return true;
	} else {
		return false;
	}
}
inline void RotamerLibrary::reset() { libraries.clear(); setup();}
inline bool RotamerLibrary::libraryExists(string _libName) const {
	if (_libName == "") {
		_libName = defaultLibrary;
	}
	return libraries.find(_libName) != libraries.end();
}
inline bool RotamerLibrary::residueExists(string _libName, string _resName) {
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
inline vector<string> RotamerLibrary::getInitAtoms(string _libName, string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.initAtoms;
	}
	return vector<string>();

}
inline vector<RotamerLibrary::InternalCoorDefi> RotamerLibrary::getInternalCoorDefinition(string _libName, string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.defi;
	} else {
		return vector<InternalCoorDefi>();
	}
}
/*
inline void RotamerLibrary::getInternalCoorDefinition(string _libName, string _resName, vector<unsigned int> & _type, vector<vector<string> > & _atomNames, vector<vector<int> > & _resnumCorrectors) {
	_type.clear();
	_atomNames.clear();
	_resnumCorrectors.clear();
	if (residueExists(_libName, _resName)) {
		for (vector<InternalCoorDefi>::iterator k=lastFoundRes->second.defi.begin(); k!=lastFoundRes->second.defi.end(); k++) {
			_type.push_back(k->type);
			_atomNames.push_back(k->atomNames);
			_resnumCorrectors.push_back(k->resnumCorrectors);
		}
	}
}
*/
inline map<string, RotamerLibrary::RotamerBuildingIC> RotamerLibrary::getRotamerBuildingIC(string _libName, string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.buildingInstructions;
	} else {
		return map<string, RotamerBuildingIC>();
	}
}
inline vector<vector<double> > RotamerLibrary::getInternalCoor(string _libName, string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.internalCoor;
	}
	return vector<vector<double> >();
}
inline string RotamerLibrary::getDefaultLibrary() {
	return defaultLibrary;
}
inline set<string> RotamerLibrary::getResList(string _libName) {

	if(libraryExists(_libName)) {
		set<string> list;
		for(map<string,Res>::iterator i = libraries[_libName].begin(); i != libraries[_libName].end(); i++) {
			list.insert(i->first);
			//cout << "UUUUU inserting " << i->first << endl;
		}
		return list;
	} else {
		return set<string> ();
	}
}


inline set<string> RotamerLibrary::getAllResList() {
	set<string> list;
	 //cout << "UUUUU resList.size" << resList.size() << endl;
	for (map <string, map<string,Res> >::iterator i = libraries.begin(); i != libraries.end(); i++) {
		 // cout << "UUUUU size of vector" << (*i).size() << endl;
		for(map<string,Res>::iterator j = i->second.begin(); j != i->second.end(); j++) {
				
			list.insert(j->first);
		}
	}
	return list;
}


inline unsigned int RotamerLibrary::size(string _libName, string _resName) {
	if (residueExists(_libName, _resName)) {
		return lastFoundRes->second.internalCoor.size();
	} else {
		return 0;
	}
}



#endif
