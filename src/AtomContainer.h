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


#ifndef ATOMCONTAINER_H
#define ATOMCONTAINER_H

#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBReader.h"
#include "PDBWriter.h"


namespace MSL { 
class Residue;

/**************************************************
 *   TO DO
 *
 *   Reorganize the exist and getAtom(string) and the internal
 *   atom map to use the AtomId and AtomOfIdentityId, with
 *   some support for accepting both (with an internal check)
 *   DONE
 * 
 *   Read multi model PDB with multiple atom coordinates instead
 *   of adding atoms multiple times
 *   DONE
 **************************************************/

class AtomContainer {
	public:
		AtomContainer();
		AtomContainer(const AtomPointerVector & _atoms);
		AtomContainer(const AtomContainer & _AC);
		~AtomContainer();
		void operator=(const AtomContainer & _AC); // assignment
		
		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		/* ADD ATOMS TO THE END */
		void addAtom(const Atom & _atom);
		void addAtom(std::string _atomId, const CartesianPoint & _coor=CartesianPoint(0.0, 0.0, 0.0), std::string _element="");
		void addAtom(std::string _atomId, double _x, double _y, double _z, std::string _element="");
		void addAtoms(const AtomPointerVector & _atoms);

		/* INSERT ATOMS AT A CERTAIN POSITION */
		void insertAtom(const Atom & _atom, unsigned int _skipPositions);
		void insertAtom(std::string _atomId, const CartesianPoint & _coor, unsigned int _skipPositions);
		void insertAtom(std::string _atomId, double _x, double _y, double _z, unsigned int _skipPositions);
		void insertAtoms(const AtomPointerVector & _atoms, unsigned int _skipPositions);

		/* REMOVE ATOMS */
		bool removeAtom(std::string _atomId);
		void removeAtom(unsigned int _n);
		void removeAllAtoms();

		
		unsigned int size() const;
		unsigned int atomSize() const;
		// the () and [] operators are redundant
		Atom & operator[](unsigned int _n);
		Atom & operator()(unsigned int _n);
		Atom & operator[](std::string _atomId); // use argument as ("A,7,CA");
		Atom & operator()(std::string _atomId); // use argument as ("A,7,CA");
		Atom & getAtom(unsigned int _n);
		Atom & getAtom(std::string _atomId);
		AtomPointerVector & getAtomPointers();

		bool atomExists(std::string _atomId);
		Atom & getLastFoundAtom();
		void updateAtomMap(Atom * _atom);

		/* I/O */
		bool readPdb(std::string _filename);
                bool readPdb(std::stringstream& _filename);
		bool writePdb(std::string _filename);

		// print the atom container using the AtomPointerVector toString
		std::string toString() const;
		friend std::ostream & operator<<(std::ostream &_os, const AtomContainer & _atomContainer)  {_os << _atomContainer.toString(); return _os;};

		/************************************************
		 *  If true, when atoms are added, if they already
		 *  exist, an alt coor is added. Otherwise atoms
		 *  are added even if they exist
		 ************************************************/
		void setAddAtomsAsAltCoors(bool _flag);
		bool getAddAtomsAsAltCoors() const;

		/************************************************
		 *  Setting the active conformation (return false if 
		 *  not all atoms have enough conformations
		 ************************************************/
		bool setActiveConformation(unsigned int _i);


		/***************************************************
		 *  Saving coordinates to buffers:
		 *
		 *  coordinates can be saved to named buffers (string _coordName),
		 *  and copied back from them
		 *
		 *  The difference between save coordinates to a buffer, and 
		 *  having multiple alternate coor is that the saved coord 
		 *  are simply a buffer that can be restored
		 *
		 *  Coor can be saved to buffer with two different commands:
		 *    saveCoor:
		 *      - saveCoor saves ONLY the current coor
		 *      - when restored with applySavedCoor, a buffer created with
		 *        saveCoor will replace the CURRENT coorinate only
		 *    saveAltCoor:
		 *      - saveAltCoor saves ALL alternative coordinates and
		 *        also remembers what was the current coordinate
		 *      - when restored with the same applySavedCoor, a buffer
		 *        created with saveAltCoor will wipe off all alternative
		 *        cordinates and recreate the situation that was present
		 *        when the buffer was saved
		 *
		 *  More details in Atom.h
		 ***************************************************/
		void saveCoor(std::string _coordName);
		void saveAltCoor(std::string _coordName);
		bool applySavedCoor(std::string _coordName);
		void clearSavedCoor(std::string _coordName="");		

	private:
		void setup();
                void setup(std::stringstream& _str);
		void copy(const AtomContainer & _AC);
		void deletePointers();		
		void reset();

		AtomPointerVector atoms;
		std::map<std::string, Atom*> atomMap;
		std::map<std::string, Atom*> atomMapWithIdentities;
		std::map<std::string, Atom*>::iterator found; // without Identity
		//std::map<std::string, Atom*>::iterator found2; // with Identity
		
		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd

		PDBReader * pdbReader;
		PDBWriter * pdbWriter;

		bool addAtomsAsAltCoors_flag;
		
};
// inlined functions
inline unsigned int AtomContainer::size() const {return atoms.size();}
inline unsigned int AtomContainer::atomSize() const {return atoms.size();}
inline void AtomContainer::setNameSpace(std::string _nameSpace) {
	nameSpace = _nameSpace;
}

inline std::string AtomContainer::getNameSpace() const{
	return nameSpace;
}

inline Atom & AtomContainer::operator[](unsigned int _n) {return *(atoms[_n]);}
inline Atom & AtomContainer::operator()(unsigned int _n) {return *(atoms[_n]);}
inline Atom & AtomContainer::operator[](std::string _atomId) { atomExists(_atomId); return getLastFoundAtom(); }
inline Atom & AtomContainer::operator()(std::string _atomId) { atomExists(_atomId); return getLastFoundAtom(); }
inline Atom & AtomContainer::getAtom(unsigned int _n) {return *atoms[_n];}
inline AtomPointerVector & AtomContainer::getAtomPointers() {return atoms;}
inline Atom & AtomContainer::getLastFoundAtom() { return *(found->second);}

inline bool AtomContainer::readPdb(std::string _filename) {
	reset();
	if (pdbReader->open(_filename)) {
		if (pdbReader->read()) {
			addAtoms(pdbReader->getAtomPointers());
			return true;
		} else {
			return false;
		}
		pdbReader->close();
	} else {
		return false;
	}
}
inline bool AtomContainer::readPdb(std::stringstream& _str) { deletePointers(); setup(_str); if (!pdbReader->open(_str) || !pdbReader->read()) return false; addAtoms(pdbReader->getAtomPointers()); return true; }
inline bool AtomContainer::writePdb(std::string _filename) {if (!pdbWriter->open(_filename)) return false; bool result = pdbWriter->write(atoms); pdbWriter->close();return result;}
inline std::string AtomContainer::toString() const {return atoms.toString();}
inline void AtomContainer::saveCoor(std::string _coordName) {atoms.saveCoor(_coordName);}
inline void AtomContainer::saveAltCoor(std::string _coordName) {atoms.saveAltCoor(_coordName);}
inline bool AtomContainer::applySavedCoor(std::string _coordName) {return atoms.applySavedCoor(_coordName);}
inline void AtomContainer::clearSavedCoor(std::string _coordName) {atoms.clearSavedCoor(_coordName);}
inline bool AtomContainer::setActiveConformation(unsigned int _i) {
	bool out = true;
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (!(*k)->setActiveConformation(_i)) {
			out = false;
		}
	}
	return out;
}
inline void AtomContainer::setAddAtomsAsAltCoors(bool _flag) { addAtomsAsAltCoors_flag = _flag;}
inline bool AtomContainer::getAddAtomsAsAltCoors() const {return addAtomsAsAltCoors_flag;}


}

#endif
