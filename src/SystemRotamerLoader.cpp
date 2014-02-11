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

#include "SystemRotamerLoader.h"

using namespace MSL;
using namespace std;
#include "MslOut.h"

static MslOut MSLOUT("SystemRotamerLoader");

SystemRotamerLoader::SystemRotamerLoader() {
	setup(NULL, "");
}

SystemRotamerLoader::SystemRotamerLoader(System & _sys) {
	setup(&_sys, "");
}

SystemRotamerLoader::SystemRotamerLoader(System & _sys, string _libraryFile) {
	setup(&_sys, _libraryFile);
}

SystemRotamerLoader::SystemRotamerLoader(const SystemRotamerLoader & _sysrotload) {
}

SystemRotamerLoader::~SystemRotamerLoader() {
	deletePointers();
}

void SystemRotamerLoader::setup(System * _pSys, string _libraryFile) {
	pRotLib = new RotamerLibrary;
	pRotRead = new RotamerLibraryReader(pRotLib);
	rotamerLibraryFile = _libraryFile;
	pSystem = _pSys;
	deleteRotLib_flag = true;
	if (_libraryFile != "") {
		if (!readRotamerLibraryFile(_libraryFile)) {
			cerr << "WARNING 3840: error reading the rotamer library file " << _libraryFile << " in void SystemRotamerLoader::setup(System * _pSys, string _libraryFile)" << endl;
		}
	}
}

void SystemRotamerLoader::deletePointers() {
	if (deleteRotLib_flag) {
		delete pRotLib;
	}
	delete pRotRead;
}

bool SystemRotamerLoader::readRotamerLibraryFile(string _libraryFile) {
	rotamerLibraryFile = _libraryFile;
	if (!pRotRead->open(_libraryFile)) {
		cerr << "WARNING 3831: cannot open rotamer library " << _libraryFile << " file in void SystemRotamerLoader::setup(System * _pSys, string _libraryFile)" << endl;
		return false;
	} else {
		if (!pRotRead->read()) {
			cerr << "WARNING 3836: cannot read rotamer library file " << _libraryFile << " in void SystemRotamerLoader::setup(System * _pSys, string _libraryFile)" << endl;
			return false;
		}
		pRotRead->close();
	}
	return true;
}


/*
bool SystemRotamerLoader::loadRotamers(string _chainId, string _resNumAndIcode, string _rotLib, string _resName, int _start, int _end, bool _keepOldRotamers) {
	cerr << "DEPRECATED bool SystemRotamerLoader::loadRotamers(string _chainId, string _resNumAndIcode, string _rotLib, string _resName, int _start, int _end, bool _keepOldRotamers), use bool loadRotamers(std::string _positionId, std::string _rotLib, std::string _residue, int _start, int _end, bool _keepOldRotamers=false) instead" << endl; 
	string positionId = _chainId + (string)"," + _resNumAndIcode;
	return loadRotamers(positionId, _rotLib, _resName, _start, _end, _keepOldRotamers);
}
*/

bool SystemRotamerLoader::loadRotamers(string _positionId, string _rotLib, string _resName, int _start, int _end, bool _keepOldRotamers) {
	cerr << "DEPRECATED bool SystemRotamerLoader::loadRotamers(string _positionId, string _rotLib, string _resName, unsigned int _start, unsigned int _end, bool _keepOldRotamers), use bool SystemRotamerLoader::loadRotamers(string _positionId, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) instead" << endl; 
	return loadRotamers(_positionId, _rotLib, _start,  _end, _resName, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _rotLib, string _resName, int _start, int _end, bool _keepOldRotamers) {
	cerr << "DEPRECATED bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _rotLib, string _resName, unsigned int _start, unsigned int _end, bool _keepOldRotamers), use bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) instead" << endl; 
	return loadRotamers(_resIndex, _resName, _start, _end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _rotLib, string _resName, int _start, int _end, bool _keepOldRotamers) {
	cerr << "DEPRECATED bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _rotLib, string _resName, unsigned int _start, unsigned int _end, bool _keepOldRotamers), use bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) instead" << endl; 
	return loadRotamers(_pPos, _resName, _start, _end, _rotLib, _keepOldRotamers);
}


bool SystemRotamerLoader::loadRotamers(string _positionId, string _resName, string _levelName, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = pRotLib->getLevel(_levelName,_resName) - 1;
	return loadRotamers(_positionId, _resName, start,  end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _resName, string _levelName, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = pRotLib->getLevel(_levelName,_resName) - 1;
	return loadRotamers(_resIndex, _resName, start, end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, string _levelName, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = pRotLib->getLevel(_levelName,_resName) - 1;
	return loadRotamers(_pPos, _resName, start, end, _rotLib, _keepOldRotamers);
}


bool SystemRotamerLoader::loadRotamers(string _positionId, string _resName, unsigned int _numberOfRots, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = _numberOfRots - 1;
	return loadRotamers(_positionId, _resName, start,  end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _resName, unsigned int _numberOfRots, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = _numberOfRots - 1;
	return loadRotamers(_resIndex, _resName, start, end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _numberOfRots, string _rotLib, bool _keepOldRotamers) {
	unsigned int start = 0;
	unsigned int end = _numberOfRots - 1;
	return loadRotamers(_pPos, _resName, start, end, _rotLib, _keepOldRotamers);
}


bool SystemRotamerLoader::loadRotamers(string _positionId, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) {
	if (pSystem->positionExists(_positionId)) {
		// get the residue and find the index
		Position * pPos = &(pSystem->getLastFoundPosition());
		return loadRotamers(pPos, _resName, _start, _end, _rotLib, _keepOldRotamers);
	} else {
		cerr << "WARNING 58229: Position " << _positionId << " does not exist, bool SystemRotamerLoader::loadRotamers(string _positionId, string _rotLib, string _resName, unsigned int _start, unsigned int _end, bool _keepOldRotamers)" << endl;
		return false;
	}
		
}

bool SystemRotamerLoader::loadRotamers(unsigned int _resIndex, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) {
	return loadRotamers(&(pSystem->getPosition(_resIndex)), _resName, _start, _end, _rotLib, _keepOldRotamers);
}

bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers) {


	/*
	  Check that position is NOT defined.
	 */
	if (_pPos == NULL) {
		return false;
	}

	/*
	  Insure identity exists
	 */
	if (!_pPos->identityExists(_resName)){
		return false;
	}


	/*
	  Insure that the residue exists in specific rotamer library
	 */
	if (!pRotLib->residueExists(_rotLib, _resName)) {
		return false;
	}

	// Extract DEFI for this rotamer library/residue type
	vector<RotamerLibrary::InternalCoorDefi> defi = pRotLib->getInternalCoorDefinition(_rotLib, _resName);

	// Extract ICVALUES for this rotamer library/residue type
	vector<vector<double> > icValues = pRotLib->getInternalCoor(_rotLib, _resName);

	// Make sure we have enough icValues to load from _start to _end.
	if (_start > _end || _end >= icValues.size()) {
		cerr << "WARNING 58229: Indeces " << _start << " and " << _end << " are out of range for residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << "; icValues.size(): "<<icValues.size()<<" in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
		return false;
	}


	// remember what was the current identity to set it back later
	string currentIdentity = _pPos->getResidueName();


	// get the residue and find the index
	_pPos->setActiveIdentity(_resName,false); // Set to proper identity. Don't apply to linked positions

	AtomPointerVector atoms = _pPos->getAtomPointers();

	map<string, Atom*> atomMap;
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		atomMap[(*k)->getName()] = *k;
	}

	IcTable * pIcTable = &(pSystem->getIcTable());
	vector<vector<Atom*> > defiAtoms;
	for (vector<RotamerLibrary::InternalCoorDefi>::iterator defiItr=defi.begin(); defiItr != defi.end(); defiItr++) {
		defiAtoms.push_back(vector<Atom*>());
		for (unsigned int j=0; j<defiItr->atomNames.size(); j++) {
			map<string, Atom*>::iterator found = atomMap.find(defiItr->atomNames[j]);
			if (found == atomMap.end()) {
				cerr << "WARNING 58239: Definition atom " << defiItr->atomNames[j] << " not found in residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
				return false;
			}
			defiAtoms.back().push_back(found->second);

		}
	}


	// make sure that _resName is defined at the position
	if (_pPos->getResidueName() != _resName) {
		if (_pPos->identityExists(_resName)) {
			_pPos->setActiveIdentity(_resName);
		} else {
			cerr << "WARNING 58234: Identity " << _resName << " does not exists at position " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
			return false;
		}
	}
	if(!_keepOldRotamers) {
		_pPos->getCurrentIdentity().removeAllAltConformations();
	}

	// get the atoms that need to be rebuilt
	vector<string> mobileAtoms = pRotLib->getMobileAtoms(_rotLib, _resName);
// 	cout << "UUU " << _rotLib << "/" << _resName << " has " << mobileAtoms.size() << " mobileAtoms atoms" << endl;
// 	for (vector<string>::iterator k=mobileAtoms.begin(); k<mobileAtoms.end(); k++) {
// 		cout << "UUU " << *k << endl;
//	}

	map<string, RotamerLibrary::RotamerBuildingIC> buildIC = pRotLib->getRotamerBuildingIC(_rotLib, _resName);
	vector<vector<Atom*> > rotamerBuildingICAtoms;
	for (map<string, RotamerLibrary::RotamerBuildingIC>::iterator bldItr=buildIC.begin(); bldItr!=buildIC.end(); bldItr++) {
		rotamerBuildingICAtoms.push_back(vector<Atom*>());
		for (unsigned int i=0; i<bldItr->second.atomNames.size(); i++) {
			map<string, Atom*>::iterator found = atomMap.find(bldItr->second.atomNames[i]);
			if (found != atomMap.end()) {
				rotamerBuildingICAtoms.back().push_back(found->second);
			} else {
				cerr << "WARNING 58239: Atom " << mobileAtoms[i] << " not found in residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
				return false;
			}
		}
	}


	vector<Atom*> initAtomPointers;
	for (unsigned int i=0; i<mobileAtoms.size(); i++) {
		//cout << "Working on InitAtom: "<<mobileAtoms[i]<<endl;
		// for each atom to be rebuilt...
		map<string, Atom*>::iterator found = atomMap.find(mobileAtoms[i]);
		if (found == atomMap.end()) {
			cerr << "WARNING 58239: Atom " << mobileAtoms[i] << " not found in residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
			return false;
		}
		initAtomPointers.push_back(found->second);

		/****************************************************************
		 * Check if the atom has an IC entry that can build...
		 ****************************************************************/
		bool icFoundForAtom = false;
		vector<IcEntry*> & icEntries = found->second->getIcEntries();
		for (vector<IcEntry*>::iterator k=icEntries.begin(); k!=icEntries.end(); k++) {

			(*k)->clearAllBuffers(); // remove any saved IC values
			for (vector<vector<Atom*> >::iterator bldAtomsItr=rotamerBuildingICAtoms.begin(); bldAtomsItr!=rotamerBuildingICAtoms.end(); bldAtomsItr++) {
				if ((*k)->match((*bldAtomsItr)[0], (*bldAtomsItr)[1], (*bldAtomsItr)[2], (*bldAtomsItr)[3])) {
					icFoundForAtom = true;
				}
			}
			if (icFoundForAtom) {
				break;
			}
		}
		if (!icFoundForAtom) {

			/************************************************************
			 *   ...if the IC entry does not exist:
			 *   add a new IC entry that can build the atom from the
			 *   rotamer library definitions
			 ************************************************************/
			for (vector<RotamerLibrary::InternalCoorDefi>::iterator defiItr=defi.begin(); defiItr != defi.end(); defiItr++) {

				// Skip non-dihedral entries..
				if (defiItr->atomNames.size() != 4){
					continue;
				}
				if (defiItr->atomNames[3] == mobileAtoms[i]) {
					// this IC can build the init atom
					Atom * pAtom1 = NULL;
					Atom * pAtom2 = NULL;
					Atom * pAtom3 = NULL;
					Atom * pAtom4 = NULL;
					if (_pPos->atomExists(defiItr->atomNames[0])) {
						pAtom1 = &(_pPos->getAtom(defiItr->atomNames[0]));
					} else {
						continue;
					}
					if (_pPos->atomExists(defiItr->atomNames[1])) {
						pAtom2 = &(_pPos->getAtom(defiItr->atomNames[1]));
					} else {
						continue;
					}
					if (_pPos->atomExists(defiItr->atomNames[2])) {
						pAtom3 = &(_pPos->getAtom(defiItr->atomNames[2]));
					} else {
						continue;
					}
					if (_pPos->atomExists(defiItr->atomNames[3])) {
						pAtom4 = &(_pPos->getAtom(defiItr->atomNames[3]));
					} else {
						continue;
					}
					bool improper_flag = false;
					if (defiItr->type == 3) {
						improper_flag = true;
					}

					pSystem->addIcEntry(pAtom1, pAtom2, pAtom3, pAtom4, 0.0, 0.0, 0.0, 0.0, 0.0, improper_flag);
					pSystem->addIcEntry(pAtom4, pAtom2, pAtom3, pAtom1, 0.0, 0.0, 0.0, 0.0, 0.0, improper_flag);
					icFoundForAtom = true;
					break;
				}
			}
		}
		if (!icFoundForAtom) {
			// we could not add the IC either, we fail
			cerr << "WARNING 58249: could not create building IC entry for atom " << mobileAtoms[i] << " in residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
			return false;
		}

	}
	/***********************************************************************
	 *
	 *   Now apply the values for the desired rotamers and rebuild,
	 *   saving alternate conformations
	 *   
	 ************************************************************************/

	for (unsigned int i=_start; i<=_end; i++) {

		if (icValues[i].size() != defiAtoms.size()) {
			cerr << "WARNING 58244: Mismatching number of definitions and values in residue " << _pPos->getChainId() << " " << _pPos->getResidueNumber() << " " << _resName << " in bool SystemRotamerLoader::loadRotamers(Position * _pPos, string _resName, unsigned int _start, unsigned int _end, string _rotLib, bool _keepOldRotamers)" << endl;
			return false;
		}
		
		for (unsigned int j=0; j<icValues[i].size(); j++) {
			if (defiAtoms[j].size() == 2) {
				//cout << "IC BOND["<<i<<"]: "<<(*defiAtoms[j][0]).getName()<<" "<<(*defiAtoms[j][1]).getName()<<" "<<icValues[i][j]<<endl;
				pIcTable->editBond(defiAtoms[j][0], defiAtoms[j][1], icValues[i][j]);
				pIcTable->editBond(defiAtoms[j][1], defiAtoms[j][0], icValues[i][j]);
			} else if (defiAtoms[j].size() == 3) {
				//cout << "IC ANGLE["<<i<<"]: "<<(*defiAtoms[j][0]).getName()<<" "<<(*defiAtoms[j][1]).getName()<<" "<<(*defiAtoms[j][2]).getName()<<" "<<icValues[i][j]<<endl;
				pIcTable->editAngle(defiAtoms[j][0], defiAtoms[j][1], defiAtoms[j][2], icValues[i][j]);
				pIcTable->editAngle(defiAtoms[j][2], defiAtoms[j][1], defiAtoms[j][0], icValues[i][j]);
			} else if (defiAtoms[j].size() == 4) {
				//cout << "IC DIHEDRAL["<<i<<"]: "<<(*defiAtoms[j][0]).getName()<<" "<<(*defiAtoms[j][1]).getName()<<" "<<(*defiAtoms[j][2]).getName()<<" "<<(*defiAtoms[j][3]).getName()<<" "<<icValues[i][j]<<endl;
				pIcTable->editDihedral(defiAtoms[j][0], defiAtoms[j][1], defiAtoms[j][2], defiAtoms[j][3], icValues[i][j]);
				pIcTable->editDihedral(defiAtoms[j][3], defiAtoms[j][1], defiAtoms[j][2], defiAtoms[j][0], -icValues[i][j]);
			}

		}
		for (vector<Atom*>::iterator k=initAtomPointers.begin(); k!=initAtomPointers.end(); k++) {
			//if (i>0 && (*k)->getNumberOfAltConformations() != 0) {
			if(i > _start || _keepOldRotamers) {
				(*k)->addAltConformation();
				(*k)->setActiveConformation((*k)->getNumberOfAltConformations()-1);
			} 
			//(*k)->setActiveConformation(i);
			(*k)->wipeCoordinates();
		}

		

		for (vector<Atom*>::iterator k=initAtomPointers.begin(); k!=initAtomPointers.end(); k++) {
			//cout << "Trying to build from IC: "<<(*k)->getName()<<endl;
			(*k)->buildFromIc();
		}
		//pSystem->printIcTable();
	}

	// restore the original identity
	_pPos->setActiveIdentity(currentIdentity);
	pSystem->updateVariablePositions();
	return true;


}
