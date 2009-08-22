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

#ifndef SYSTEM_H
#define SYSTEM_H

#include <string>
#include <vector>

#include "Chain.h"
#include "IcTable.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "EnergySet.h"

using namespace std;

class System {
	public:
		System();
		System(const Chain & _chain);
		System(const vector<Chain> & _chains);
		System(const AtomVector & _atoms);
		System(const System & _system);
		~System();

		void operator=(const System & _system); // assignment

		void addChain(const Chain & _chain, string _chainId="");
		bool removeChain(string _chainId);
		void removeAllChains();
		bool duplicateChain(string _chainId, string _newChainId="");
		bool duplicateChain(size_t _n, string _newChainId="");

		void addAtoms(const AtomVector & _atoms);
		
		EnergySet* getEnergySet();
		double calcEnergy(bool _activeOnly=true);
		double calcEnergy(string _selection, bool _activeOnly=true);
		double calcEnergy(string _selection1, string _selection2, bool _activeOnly=true);

		void setNameSpace(string _nameSpace);
		string getNameSpace() const;

		Chain & operator()(string _chainId); // return chain by chainId ("A")
		unsigned int size() const; // number of chains
		unsigned int atomSize() const; // number of active atoms
		unsigned int allAtomSize() const; // number of atoms active and inactive
		unsigned int positionSize() const; // number of positions
		unsigned int linkedPositionSize() const; // number of positions linked, labeled "MASTER"
		unsigned int slavePositionSize() const; // number of positions linked, labeled "SLAVE"
		
		vector<Chain*> & getChains();
		Chain & getChain(size_t _n);
		Chain & getChain(string _chainId);

		unsigned int residueSize() const; // number of positions
		vector<Position*> & getPositions();
		Position & getPosition(size_t _n);
		Position & getPosition(string _chainId, int _resNum);
		Position & getPosition(string _chainId, string _resNumAndIcode);
		Residue & getResidue(size_t _n);
		unsigned int getPositionIndex(string _chainId, int _resNum);	
		unsigned int getPositionIndex(string _chainId, string _resNumAndIcode);	
		unsigned int getPositionIndex(const Position * _pPos) const;	

		AtomVector & getAtoms();
		AtomVector & getAllAtoms();
		Atom & operator[](size_t _n);
		Atom & getAtom(size_t _n);
		
		// check the existance of chains, residues, atoms
		bool exists(string _chainId); // chain
		bool exists(string _chainId, int _resNum); // residue, by int
		bool exists(string _chainId, string _resNumAndIcode); // residue by string (possibly with insertion code
		bool exists(string _chainId, int _resNum, string _name); // atom
		bool exists(string _chainId, string _resNumAndIcode, string _name); // atom
		bool exists(string _chainId, int _resNum, string _name, string _identity); // atom specifying identity (i.e. ALA)
		bool exists(string _chainId, string _resNumAndIcode, string _name, string _identity); // atom specifying identity (i.e. ALA)

		Chain & getLastFoundChain();
		Position & getLastFoundPosition();
		Residue & getLastFoundResidue();
		Atom & getLastFoundAtom();

		/* IC TABLE FOR BUILDING AND EDITING CONFORMATION FROM INTERNAL COORDINATES */
		//vector<IcEntry*> & getIcTable();
		IcTable & getIcTable();
		bool addIcEntry(string _1_chain_resNumIcode_name, string _2_chain_resNumIcode_name, string _3_chain_resNumIcode_name, string _4_chain_resNumIcode_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag=false); 
		bool addIcEntry(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag=false); 
		bool addIcEntry(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag);
		bool seed(string _1_chain_resNumIcode_name, string _2_chain_resNumIcode_name, string _3_chain_resNumIcode_name);
		bool seed(string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name); // removes all coordinates and finds 3 atoms to seed in cartesian space

		void fillIcFromCoor();
		void buildAtoms(); // build all possible coordinates of the active atoms from the IC table (active atoms only)
		void buildAllAtoms(); // build all possible coordinates of the active atoms from the IC table (active and inactive atoms)
		void printIcTable() const;
		void saveIcToBuffer(string _name);
		void restoreIcFromBuffer(string _name);
		void clearAllIcBuffers();

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		/* I/O */
		bool readPdb(string _filename); // add atoms or alt coor
		bool writePdb(string _filename);

		unsigned int assignCoordinates(const AtomVector & _atoms,bool checkIdentity=true); // only set coordinates for existing matching atoms, return the number assigned

		/* UPDATES REQUESTED BY POSITIONS */
		void updateChainMap(Chain * _chain);
		void updateIndexing();
		void updateAllAtomIndexing();
		//void swapInActiveList(Position * _position, AtomVector & _atoms);

		// copy coordinates for specified atoms from the current identity of each position to all other identities
		void copyCoordinatesOfAtomsInPosition(vector<string> _sourcePosNames=vector<string>());
		
		void updateVariablePositions(); // either multiple identities or rotamers
		vector<unsigned int> getVariablePositions() const;  // get the index of the variable positions, need to run updateVariablePositions() first
		void setActiveRotamers(vector<unsigned int> _rots); // set the active rotamers for all variable positions

		
		// takes vector<vector<string> > which is:
		//  [0]  A_19 B_19 C_19   
		//  [1]  A_22 B_22 C_22
		// ..
		// Meaning A_19 is linked as a MASTER to B_19,C_19
		//         A_22 is linekd as a MASTER to B_22,C_22
		void setLinkedPositions(vector<vector<string> > &_linkedPositions);

	private:
		void setup();
		void copy(const System & _system);
		void reset();
		void deletePointers();
		bool findIcAtoms(Atom *& _pAtom1, Atom *& _pAtom2, Atom *& _pAtom3, Atom *& _pAtom4, string _1_chain, string _1_resNumIcode, string _1_name, string _2_chain, string _2_resNumIcode, string _2_name, string _3_chain, string _3_resNumIcode, string _3_name, string _4_chain, string _4_resNumIcode, string _4_name);

		vector<Chain*> chains;
		vector<Position*> positions;
		map<string, Chain*> chainMap;

		//vector<IcEntry*> icTable;
		IcTable icTable;

		string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		/*********************************************
		 *  We keep a list of active atoms
		 *
		 *  We also keep an index of where the chain
		 *  atoms start and end so that we can quickly
		 *  swap them when a position changes identity
		 *********************************************/
		AtomVector activeAtoms;
		AtomVector activeAndInactiveAtoms;
		bool noUpdateIndex_flag;

		map<string, Chain*>::iterator foundChain;
		EnergySet* ESet;
		PDBReader * pdbReader;
		PDBWriter * pdbWriter;
	
		vector<unsigned int> variablePositions; // this needs to be updated with updateVariablePositions()


};

// INLINED FUNCTIONS
inline void System::setNameSpace(string _nameSpace) {nameSpace = _nameSpace;}
inline string System::getNameSpace() const {return nameSpace;}
inline Chain & System::operator()(string _chainId) {return *chainMap[_chainId];}
inline unsigned int System::size() const {return chains.size();}
inline unsigned int System::atomSize() const {return activeAtoms.size();}
inline unsigned int System::allAtomSize() const {return activeAndInactiveAtoms.size();}
inline unsigned int System::positionSize() const {return positions.size();}
inline vector<Chain*> & System::getChains() {return chains;}
inline Chain & System::getChain(size_t _n) {return *(chains[_n]);}
inline Chain & System::getChain(string _chainId) {exists(_chainId); return *(foundChain->second);}
inline unsigned int System::residueSize() const {return positions.size();}
inline vector<Position*> & System::getPositions() {return positions;}
inline Position & System::getPosition(size_t _n) {return *(positions[_n]);}
inline Position & System::getPosition(string _chainId, int _resNum) {
	if (!exists(_chainId, _resNum)) {
		cerr << "ERROR 49129: Position " << _chainId << " " << _resNum << " not found in inline Position & System::getPosition(string _chainId, int _resNum)" << endl;
	}
	return foundChain->second->getLastFoundPosition();
}
inline Position & System::getPosition(string _chainId, string _resNumAndIcode) {
	if (!exists(_chainId, _resNumAndIcode)) {
		cerr << "ERROR 49134: Position " << _chainId << " " << _resNumAndIcode << " not found in inline Position & System::getPosition(string _chainId, string _resNumAndIcode)" << endl;
	}
	return foundChain->second->getLastFoundPosition();
}

inline Residue & System::getResidue(size_t _n) {return positions[_n]->getCurrentIdentity();}
inline AtomVector & System::getAtoms() {return activeAtoms;}
inline AtomVector & System::getAllAtoms() {return activeAndInactiveAtoms;}
inline Atom & System::operator[](size_t _n) {return *(activeAtoms[_n]);}
inline Atom & System::getAtom(size_t _n) {return *(activeAtoms[_n]);}
inline bool System::exists(string _chainId) {foundChain = chainMap.find(_chainId); return foundChain != chainMap.end();}
inline bool System::exists(string _chainId, int _resNum) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum);} return false;}
inline bool System::exists(string _chainId, string _resNumAndIcode) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode);} return false;}
inline bool System::exists(string _chainId, int _resNum, string _name) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum, _name);} return false;}
inline bool System::exists(string _chainId, string _resNumAndIcode, string _name) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode, _name);} return false;}
inline bool System::exists(string _chainId, int _resNum, string _name, string _identity) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum, _name, _identity);} return false;}
inline bool System::exists(string _chainId, string _resNumAndIcode, string _name, string _identity) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode, _name, _identity);} return false;}
inline EnergySet* System::getEnergySet() { return(ESet); }
inline double System::calcEnergy(bool _activeOnly) { return (ESet->calcEnergy(_activeOnly));}
inline double System::calcEnergy(string _selection, bool _activeOnly) { return (ESet->calcEnergy(_selection, _activeOnly));}
inline double System::calcEnergy(string _selection1, string _selection2, bool _activeOnly) { return (ESet->calcEnergy(_selection1, _selection2, _activeOnly));}

//inline vector<IcEntry*> & System::getIcTable() {return icTable;}
inline IcTable & System::getIcTable() {return icTable;}
inline Chain & System::getLastFoundChain() {return *(foundChain->second);}
inline Position & System::getLastFoundPosition() {return foundChain->second->getLastFoundPosition();}
inline Residue & System::getLastFoundResidue() {return foundChain->second->getLastFoundResidue();}
inline Atom & System::getLastFoundAtom() {return foundChain->second->getLastFoundAtom();}
inline void System::wipeAllCoordinates() {for (vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {(*k)->wipeAllCoordinates();}}
inline void System::buildAtoms() {
	// build only the active atoms
	for (AtomVector::iterator k=activeAtoms.begin(); k!=activeAtoms.end(); k++) {\
		(*k)->buildFromIc(true); // build only from active atoms = true
	}
}
inline void System::buildAllAtoms() {
	// build active and inactive atoms
	for (AtomVector::iterator k=activeAndInactiveAtoms.begin(); k!=activeAndInactiveAtoms.end(); k++) {
		(*k)->buildFromIc(false); // build only from active atoms = true
	}
}
inline void System::fillIcFromCoor() {for (IcTable::iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->fillFromCoor();}}
inline void System::printIcTable() const {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {cout << *(*k) << endl;}}
inline void System::saveIcToBuffer(string _name) {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->saveBuffer(_name);}}
inline void System::restoreIcFromBuffer(string _name) {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->restoreFromBuffer(_name);}}
inline void System::clearAllIcBuffers() {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->clearAllBuffers();}}
inline bool System::readPdb(string _filename) {reset(); if (!pdbReader->open(_filename) || !pdbReader->read()) return false; addAtoms(pdbReader->getAtoms()); return true;}
inline bool System::writePdb(string _filename) {if (!pdbWriter->open(_filename)) return false; return pdbWriter->write(activeAtoms);}

inline unsigned int System::getPositionIndex(string _chainId, int _resNum) {
	if (exists(_chainId, _resNum)) {
		return getPositionIndex(&getLastFoundPosition());
	} else {
		cerr << "ERROR 44198: Position " << _chainId << " " << _resNum << " not found in unsigned inline unsigned int System::getPositionIndex(string _chainId, int _resNum)" << endl;
		exit(44193);
	}
}
inline unsigned int System::getPositionIndex(string _chainId, string _resNumAndIcode) { 
	if (exists(_chainId, _resNumAndIcode)) {
		return getPositionIndex(&getLastFoundPosition());
	} else {
		cerr << "ERROR 44203: Position " << _chainId << " " << _resNumAndIcode << " not found in unsigned inline unsigned inline unsigned int System::getPositionIndex(string _chainId, string _resNumAndIcode)" << endl;
		exit(44193);
	}
}
inline void System::copyCoordinatesOfAtomsInPosition(vector<string> _sourcePosNames) {for (vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {(*k)->copyCoordinatesOfAtoms(_sourcePosNames);} }

inline void System::updateVariablePositions() {
	// update the list of positions with either multiple identities or rotamers
	variablePositions.clear();
	for (unsigned int i=0; i<positions.size(); i++) {
		if (positions[i]->size() > 1 || positions[i]->getTotalNumberOfRotamers() > 1) {
			variablePositions.push_back(i);
		}
	}
}
inline vector<unsigned int> System::getVariablePositions() const {
	// get the index of the variable positions, need to run updateVariablePositions() first
	return variablePositions;
}
inline void System::setActiveRotamers(vector<unsigned int> _rots) {
	// set the active rotamers for all variable positions
	for (unsigned int i=0; i<_rots.size(); i++) {
		if (i >= variablePositions.size()) {
			break;
		}
		positions[variablePositions[i]]->setActiveRotamer(_rots[i]);
	}
}

inline unsigned int System::linkedPositionSize() const {
	unsigned int result = 0;
	for (uint i = 0; i < positionSize();i++){
		if (positions[i]->getLinkedPositionType() == Position::MASTER){
			result++;
		}
	}

	return result;
}

inline unsigned int System::slavePositionSize() const {
	unsigned int result = 0;
	for (uint i = 0; i < positionSize();i++){
		if (positions[i]->getLinkedPositionType() == Position::SLAVE){
			result++;
		}
	}

	return result;
}

#endif
