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
#include <sys/types.h>

#include "Chain.h"
#include "IcTable.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "EnergySet.h"


namespace MSL { 
class PolymerSequence;

class System {
	public:
		System();
		System(const Chain & _chain);
		System(const std::vector<Chain> & _chains);
		System(const AtomPointerVector & _atoms);
		System(const System & _system);
		~System();

		void operator=(const System & _system); // assignment

		void addChain(const Chain & _chain, std::string _chainId="");
		bool removeChain(std::string _chainId);
		void removeAllChains();
		bool duplicateChain(std::string _chainId, std::string _newChainId="");
		bool duplicateChain(size_t _n, std::string _newChainId="");

		void addAtoms(const AtomPointerVector & _atoms);
		
		EnergySet* getEnergySet();
		/* Calculate the energies */
		double calcEnergy();
		double calcEnergy(std::string _selection);
		double calcEnergy(std::string _selection1, std::string _selection2);

		/* Calculate the energies including the interactions that inlcude atoms that belong to inactive side chains */
		double calcEnergyAllAtoms();
		double calcEnergyAllAtoms(std::string _selection);
		double calcEnergyAllAtoms(std::string _selection1, std::string _selection2);

		std::string getEnergySummary () const;
		void printEnergySummary() const;

		double calcEnergyOfSubset(std::string _subsetName);

		void saveEnergySubset(std::string _subsetName);
		void saveEnergySubset(std::string _subsetName, std::string _selection);
		void saveEnergySubset(std::string _subsetName, std::string _selection1, std::string _selection2);
		void saveEnergySubsetAllAtoms(std::string _subsetName);
		void saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection);
		void saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection1, std::string _selection2);

		void removeEnergySubset(std::string _subsetName);



		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		Chain & operator()(std::string _chainId); // return chain by chainId ("A")
		unsigned int size() const; // number of chains
		unsigned int atomSize() const; // number of active atoms
		unsigned int allAtomSize() const; // number of atoms active and inactive
		unsigned int positionSize() const; // number of positions
		unsigned int linkedPositionSize() const; // number of positions linked, labeled "MASTER"
		unsigned int slavePositionSize() const; // number of positions linked, labeled "SLAVE"
		
		std::vector<Chain*> & getChains();
		Chain & getChain(size_t _n);
		Chain & getChain(std::string _chainId);

		unsigned int residueSize() const; // number of positions
		std::vector<Position*> & getPositions();
		Position & getPosition(size_t _n);
		Position & getPosition(std::string _chainId, int _resNum);
		Position & getPosition(std::string _chainId, std::string _resNumAndIcode);
		Residue & getResidue(size_t _n);
		unsigned int getPositionIndex(std::string _chainId, int _resNum);	
		unsigned int getPositionIndex(std::string _chainId, std::string _resNumAndIcode);	
		unsigned int getPositionIndex(const Position * _pPos) const;	

		AtomPointerVector & getAtoms();
		AtomPointerVector & getAllAtoms();
		Atom & operator[](size_t _n);
		Atom & getAtom(size_t _n);
		
		// check the existance of chains, residues, atoms
		bool exists(std::string _chainId); // chain
		bool exists(std::string _chainId, int _resNum); // residue, by int
		bool exists(std::string _chainId, std::string _resNumAndIcode); // residue by std::string (possibly with insertion code
		bool exists(std::string _chainId, int _resNum, std::string _name); // atom
		bool exists(std::string _chainId, std::string _resNumAndIcode, std::string _name); // atom
		bool exists(std::string _chainId, int _resNum, std::string _name, std::string _identity); // atom specifying identity (i.e. ALA)
		bool exists(std::string _chainId, std::string _resNumAndIcode, std::string _name, std::string _identity); // atom specifying identity (i.e. ALA)

		Chain & getLastFoundChain();
		Position & getLastFoundPosition();
		Residue & getLastFoundResidue();
		Atom & getLastFoundAtom();

		/* IC TABLE FOR BUILDING AND EDITING CONFORMATION FROM INTERNAL COORDINATES */
		//std::vector<IcEntry*> & getIcTable();
		IcTable & getIcTable();
		bool addIcEntry(std::string _1_chain_resNumIcode_name, std::string _2_chain_resNumIcode_name, std::string _3_chain_resNumIcode_name, std::string _4_chain_resNumIcode_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag=false); 
		bool addIcEntry(std::string _1_chain, std::string _1_resNumIcode, std::string _1_name, std::string _2_chain, std::string _2_resNumIcode, std::string _2_name, std::string _3_chain, std::string _3_resNumIcode, std::string _3_name, std::string _4_chain, std::string _4_resNumIcode, std::string _4_name, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag=false); 
		bool addIcEntry(Atom * _pAtom1, Atom * _pAtom2, Atom * _pAtom3, Atom * _pAtom4, double _d1, double _a1, double _dihe, double _a2, double _d2, bool _improperFlag);
		bool seed(std::string _1_chain_resNumIcode_name, std::string _2_chain_resNumIcode_name, std::string _3_chain_resNumIcode_name);
		bool seed(std::string _1_chain, std::string _1_resNumIcode, std::string _1_name, std::string _2_chain, std::string _2_resNumIcode, std::string _2_name, std::string _3_chain, std::string _3_resNumIcode, std::string _3_name); // removes all coordinates and finds 3 atoms to seed in cartesian space

		void fillIcFromCoor();
		void buildAtoms(); // build all possible coordinates of the active atoms from the IC table (active atoms only)
		void buildAllAtoms(); // build all possible coordinates of the active atoms from the IC table (active and inactive atoms)
		void printIcTable() const;
		void saveIcToBuffer(std::string _name);
		void restoreIcFromBuffer(std::string _name);
		void clearAllIcBuffers();

		void wipeAllCoordinates(); // flag all active and inactive atoms as not having cartesian coordinates

		/* I/O */
		bool readPdb(std::string _filename); // add atoms or alt coor
		bool writePdb(std::string _filename);

		unsigned int assignCoordinates(const AtomPointerVector & _atoms,bool checkIdentity=true); // only set coordinates for existing matching atoms, return the number assigned

		/* UPDATES REQUESTED BY POSITIONS */
		void updateChainMap(Chain * _chain);
		void updateIndexing();
		void updateAllAtomIndexing();
		//void swapInActiveList(Position * _position, AtomPointerVector & _atoms);

		// copy coordinates for specified atoms from the current identity of each position to all other identities
		void copyCoordinatesOfAtomsInPosition(std::vector<std::string> _sourcePosNames=std::vector<std::string>());
		
		void updateVariablePositions(); // either multiple identities or rotamers
		std::vector<unsigned int> getVariablePositions() const;  // get the index of the variable positions, need to run updateVariablePositions() first
		void setActiveRotamers(std::vector<unsigned int> _rots); // set the active rotamers for all variable positions

		
		// takes std::vector<std::vector<std::string> > which is:
		//  [0]  A_19 B_19 C_19   
		//  [1]  A_22 B_22 C_22
		// ..
		// Meaning A_19 is linked as a MASTER to B_19,C_19
		//         A_22 is linekd as a MASTER to B_22,C_22
		void setLinkedPositions(std::vector<std::vector<std::string> > &_linkedPositions);


		std::string toString() const;
		friend std::ostream & operator<<(std::ostream &_os, const System & _sys)  {_os << _sys.toString(); return _os;};
		
	private:
		void setup();
		void copy(const System & _system);
		void reset();
		void deletePointers();
		bool findIcAtoms(Atom *& _pAtom1, Atom *& _pAtom2, Atom *& _pAtom3, Atom *& _pAtom4, std::string _1_chain, std::string _1_resNumIcode, std::string _1_name, std::string _2_chain, std::string _2_resNumIcode, std::string _2_name, std::string _3_chain, std::string _3_resNumIcode, std::string _3_name, std::string _4_chain, std::string _4_resNumIcode, std::string _4_name);

		std::vector<Chain*> chains;
		std::vector<Position*> positions;
		std::map<std::string, Chain*> chainMap;

		//std::vector<IcEntry*> icTable;
		IcTable icTable;

		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		/*********************************************
		 *  We keep a list of active atoms
		 *
		 *  We also keep an index of where the chain
		 *  atoms start and end so that we can quickly
		 *  swap them when a position changes identity
		 *********************************************/
		AtomPointerVector activeAtoms;
		AtomPointerVector activeAndInactiveAtoms;
		bool noUpdateIndex_flag;

		std::map<std::string, Chain*>::iterator foundChain;
		EnergySet* ESet;
		PDBReader * pdbReader;
		PDBWriter * pdbWriter;
	
		std::vector<unsigned int> variablePositions; // this needs to be updated with updateVariablePositions()

		PolymerSequence * polSeq;


};

// INLINED FUNCTIONS
inline void System::setNameSpace(std::string _nameSpace) {nameSpace = _nameSpace;}
inline std::string System::getNameSpace() const {return nameSpace;}
inline Chain & System::operator()(std::string _chainId) {return *chainMap[_chainId];}
inline unsigned int System::size() const {return chains.size();}
inline unsigned int System::atomSize() const {return activeAtoms.size();}
inline unsigned int System::allAtomSize() const {return activeAndInactiveAtoms.size();}
inline unsigned int System::positionSize() const {return positions.size();}
inline std::vector<Chain*> & System::getChains() {return chains;}
inline Chain & System::getChain(size_t _n) {return *(chains[_n]);}
inline Chain & System::getChain(std::string _chainId) {exists(_chainId); return *(foundChain->second);}
inline unsigned int System::residueSize() const {return positions.size();}
inline std::vector<Position*> & System::getPositions() {return positions;}
inline Position & System::getPosition(size_t _n) {return *(positions[_n]);}
inline Position & System::getPosition(std::string _chainId, int _resNum) {
	if (!exists(_chainId, _resNum)) {
		std::cerr << "ERROR 49129: Position " << _chainId << " " << _resNum << " not found in inline Position & System::getPosition(std::string _chainId, int _resNum)" << std::endl;
	}
	return foundChain->second->getLastFoundPosition();
}
inline Position & System::getPosition(std::string _chainId, std::string _resNumAndIcode) {
	if (!exists(_chainId, _resNumAndIcode)) {
		std::cerr << "ERROR 49134: Position " << _chainId << " " << _resNumAndIcode << " not found in inline Position & System::getPosition(std::string _chainId, std::string _resNumAndIcode)" << std::endl;
	}
	return foundChain->second->getLastFoundPosition();
}

inline Residue & System::getResidue(size_t _n) {return positions[_n]->getCurrentIdentity();}
inline AtomPointerVector & System::getAtoms() {return activeAtoms;}
inline AtomPointerVector & System::getAllAtoms() {return activeAndInactiveAtoms;}
inline Atom & System::operator[](size_t _n) {return *(activeAtoms[_n]);}
inline Atom & System::getAtom(size_t _n) {return *(activeAtoms[_n]);}
inline bool System::exists(std::string _chainId) {foundChain = chainMap.find(_chainId); return foundChain != chainMap.end();}
inline bool System::exists(std::string _chainId, int _resNum) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum);} return false;}
inline bool System::exists(std::string _chainId, std::string _resNumAndIcode) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode);} return false;}
inline bool System::exists(std::string _chainId, int _resNum, std::string _name) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum, _name);} return false;}
inline bool System::exists(std::string _chainId, std::string _resNumAndIcode, std::string _name) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode, _name);} return false;}
inline bool System::exists(std::string _chainId, int _resNum, std::string _name, std::string _identity) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNum, _name, _identity);} return false;}
inline bool System::exists(std::string _chainId, std::string _resNumAndIcode, std::string _name, std::string _identity) {foundChain=chainMap.find(_chainId); if (foundChain != chainMap.end()) {return foundChain->second->exists(_resNumAndIcode, _name, _identity);} return false;}
inline EnergySet* System::getEnergySet() { return(ESet); }
inline double System::calcEnergy() { return (ESet->calcEnergy());}
inline double System::calcEnergy(std::string _selection) { return (ESet->calcEnergy(_selection));}
inline double System::calcEnergy(std::string _selection1, std::string _selection2) { return (ESet->calcEnergy(_selection1, _selection2));}
inline double System::calcEnergyAllAtoms() { return ESet->calcEnergyAllAtoms(); }
inline double System::calcEnergyAllAtoms(std::string _selection) { return ESet->calcEnergyAllAtoms(_selection); }
inline double System::calcEnergyAllAtoms(std::string _selection1, std::string _selection2) { return ESet->calcEnergyAllAtoms(_selection1, _selection2); }
inline std::string System::getEnergySummary () const {return ESet->getSummary();}
inline void System::printEnergySummary() const {ESet->printSummary();}
inline double System::calcEnergyOfSubset(std::string _subsetName) { return ESet->calcEnergyOfSubset(_subsetName); }
inline void System::saveEnergySubset(std::string _subsetName) { ESet->saveEnergySubset(_subsetName); }
inline void System::saveEnergySubset(std::string _subsetName, std::string _selection) { ESet->saveEnergySubset(_subsetName, _selection); }
inline void System::saveEnergySubset(std::string _subsetName, std::string _selection1, std::string _selection2) { ESet->saveEnergySubset(_subsetName, _selection1, _selection2); }
inline void System::saveEnergySubsetAllAtoms(std::string _subsetName) { ESet->saveEnergySubsetAllAtoms(_subsetName); }
inline void System::saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection) { ESet->saveEnergySubsetAllAtoms(_subsetName, _selection); }
inline void System::saveEnergySubsetAllAtoms(std::string _subsetName, std::string _selection1, std::string _selection2) { ESet->saveEnergySubsetAllAtoms(_subsetName, _selection1, _selection2); }
inline void System::removeEnergySubset(std::string _subsetName) { ESet->removeEnergySubset(_subsetName); }

//inline std::vector<IcEntry*> & System::getIcTable() {return icTable;}
inline IcTable & System::getIcTable() {return icTable;}
inline Chain & System::getLastFoundChain() {return *(foundChain->second);}
inline Position & System::getLastFoundPosition() {return foundChain->second->getLastFoundPosition();}
inline Residue & System::getLastFoundResidue() {return foundChain->second->getLastFoundResidue();}
inline Atom & System::getLastFoundAtom() {return foundChain->second->getLastFoundAtom();}
inline void System::wipeAllCoordinates() {for (std::vector<Chain*>::iterator k=chains.begin(); k!=chains.end(); k++) {(*k)->wipeAllCoordinates();}}
inline void System::buildAtoms() {
	// build only the active atoms
	for (AtomPointerVector::iterator k=activeAtoms.begin(); k!=activeAtoms.end(); k++) {\
		(*k)->buildFromIc(true); // build only from active atoms = true
	}
}
inline void System::buildAllAtoms() {
	// build active and inactive atoms
	for (AtomPointerVector::iterator k=activeAndInactiveAtoms.begin(); k!=activeAndInactiveAtoms.end(); k++) {
		(*k)->buildFromIc(false); // build only from active atoms = true
	}
}
inline void System::fillIcFromCoor() {for (IcTable::iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->fillFromCoor();}}
inline void System::printIcTable() const {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {std::cout << *(*k) << std::endl;}}
inline void System::saveIcToBuffer(std::string _name) {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->saveBuffer(_name);}}
inline void System::restoreIcFromBuffer(std::string _name) {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->restoreFromBuffer(_name);}}
inline void System::clearAllIcBuffers() {for (IcTable::const_iterator k=icTable.begin(); k!=icTable.end(); k++) {(*k)->clearAllBuffers();}}
inline bool System::readPdb(std::string _filename) {reset(); if (!pdbReader->open(_filename) || !pdbReader->read()) return false; addAtoms(pdbReader->getAtoms()); return true;}
inline bool System::writePdb(std::string _filename) {if (!pdbWriter->open(_filename)) return false; bool result = pdbWriter->write(activeAtoms); pdbWriter->close();return result;}

inline unsigned int System::getPositionIndex(std::string _chainId, int _resNum) {
	if (exists(_chainId, _resNum)) {
		return getPositionIndex(&getLastFoundPosition());
	} else {
		std::cerr << "ERROR 44198: Position " << _chainId << " " << _resNum << " not found in unsigned inline unsigned int System::getPositionIndex(std::string _chainId, int _resNum)" << std::endl;
		exit(44193);
	}
}
inline unsigned int System::getPositionIndex(std::string _chainId, std::string _resNumAndIcode) { 
	if (exists(_chainId, _resNumAndIcode)) {
		return getPositionIndex(&getLastFoundPosition());
	} else {
		std::cerr << "ERROR 44203: Position " << _chainId << " " << _resNumAndIcode << " not found in unsigned inline unsigned inline unsigned int System::getPositionIndex(std::string _chainId, std::string _resNumAndIcode)" << std::endl;
		exit(44193);
	}
}
inline void System::copyCoordinatesOfAtomsInPosition(std::vector<std::string> _sourcePosNames) {for (std::vector<Position*>::iterator k=positions.begin(); k!=positions.end(); k++) {(*k)->copyCoordinatesOfAtoms(_sourcePosNames);} }

inline void System::updateVariablePositions() {
	// update the list of positions with either multiple identities or rotamers
	variablePositions.clear();
	for (unsigned int i=0; i<positions.size(); i++) {
		if (positions[i]->size() > 1 || positions[i]->getTotalNumberOfRotamers() > 1) {
			variablePositions.push_back(i);
		}
	}
}
inline std::vector<unsigned int> System::getVariablePositions() const {
	// get the index of the variable positions, need to run updateVariablePositions() first
	return variablePositions;
}
inline void System::setActiveRotamers(std::vector<unsigned int> _rots) {
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

/* Calculate the energies including the interactions that inlcude atoms that belong to inactive side chains */

}

#endif
