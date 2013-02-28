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


#ifndef ATOMGROUP_H
#define ATOMGROUP_H

#include "AtomPointerVector.h"


namespace MSL { 
class Residue;
class Position;
class Chain;
class System;

class AtomGroup : public AtomPointerVector {
	public:
		AtomGroup();
		AtomGroup(Residue * _pParentResidue);
		~AtomGroup();

		/*
		 * should we add copy constructors?
		 * where would the pointer inside the
		 * atom point? the old AtomGroup or the
		 * new one?
		AtomGroup(const AtomPointerVector & _AV);
		AtomGroup(const AtomGroup & _AG);
		void operator=(const AtomGroup & _AG); // assignment
		 */
		
		// if the stamp is identical to the current stamp, a precalculated
		// center is given, otherwise it is calculated and the result is cached
	//	CartesianPoint getGeometricCenter(unsigned int _stamp=0);

		void setResidueName(std::string _resname);
		std::string getResidueName() const;

		void setResidueNumber(int _resnum);
		int getResidueNumber() const;

		void setResidueIcode(std::string _icode);
		std::string getResidueIcode() const;

		void setChainId(std::string _chainId);
		std::string getChainId() const;

		void setNameSpace(std::string _nameSpace);
		std::string getNameSpace() const;

		void setGroupNumber(unsigned int _groupNum);
		unsigned int getGroupNumber() const;

		unsigned int getIdentityIndex(); // return the index of its parent identity in the position

		void push_back(Atom * _atom);

		void setParentResidue(Residue * _parent);
		Residue * getParentResidue() const;
		Position * getParentPosition() const;
		Chain * getParentChain() const;
		System * getParentSystem() const;

		void updateResidueMap(Atom * _atom);

		/***************************************************
		 * As atoms have alternate conformations, residues can
		 * have alternate identities.
		 *
		 * Atoms can become not active when they belonging to 
		 * and identity that is not currently active (i.e. a
		 * VAL atom when LEU is active at the position)
		 ***************************************************/
		bool getActive() const;
		bool getHidden() const;

		/***************************************************
		 *  Ask the chain if this position is N- or C-terminal
		 ***************************************************/
		bool isPositionNterminal() const;
		bool isPositionCterminal() const;

	private:

		Residue * pParentResidue;
		
		// allows to save the geometric center and return the value in memory (cachedCenter) without 
		// recalculating it if the same stamp is given
		unsigned int stamp; 
		CartesianPoint cachedCenter;

		std::string nameSpace;  // pdb, charmm19, etc., mainly for name converting upon writing a pdb or crd
	
		unsigned int groupNumber;

		std::string residueName;
		int residueNumber;
		std::string residueIcode;
		std::string chainId;

		
};

inline Residue * AtomGroup::getParentResidue() const {return pParentResidue;}
inline void AtomGroup::push_back(Atom * _atom) {
	_atom->setParentGroup(this);
	AtomPointerVector::push_back(_atom);
}
inline void AtomGroup::setGroupNumber(unsigned int _groupNum) {groupNumber = _groupNum;}
inline void AtomGroup::setParentResidue(Residue * _parent) {pParentResidue = _parent;}

}

#endif
