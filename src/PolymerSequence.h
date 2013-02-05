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


#ifndef POLYMERSEQUENCE_H
#define POLYMERSEQUENCE_H



// MSL INCLUDES
#include "AtomPointerVector.h"
#include "System.h"
#include "MslTools.h"
#include "PDBReader.h"

// STL INCLUDES
#include <iostream>
#include <string>
#include <vector>


namespace MSL { 
	class PolymerSequence {
		public:
			PolymerSequence();
			PolymerSequence(std::string _sequence);
			PolymerSequence(System &_sys);
			PolymerSequence(const AtomPointerVector &_atoms);
			PolymerSequence(const PolymerSequence & _seq);
			PolymerSequence(System &_sys, std::vector<std::pair<std::string,std::string> > &_addTerminalResidues); // add additional residues to polymer sequence
			PolymerSequence(System &_sys, std::map<std::string,std::map<int,int> > &variablePositionMap, std::vector<std::vector<std::string> > &_identitesAtVariablePositions);
			PolymerSequence(System &_sys, std::map<std::string,int> &variablePositionMap, std::vector<std::vector<std::string> > &_identitesAtVariablePositions);
			~PolymerSequence();

			void operator=(const PolymerSequence & _seq);
			friend std::ostream & operator<<(std::ostream &_os, const PolymerSequence & _seq)  {_os << _seq.toString(); return _os;};

			void setSequence(std::string _sequence);
			void setSequence(System &_sys);
			void setSequence(const AtomPointerVector &_atoms);
			bool readSequenceFromPDB(std::string _pdbfile);

			unsigned int size() const;
			unsigned int chainSize(unsigned int _chainIndex) const;
			unsigned int positionSize(unsigned int _chainIndex, unsigned int _positionIndex) const;
			void setChainId(unsigned int _chainIndex, std::string _chainId);
			void setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, int _resnum);
			void setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, std::string _resnum);
			void setPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex, std::string _resname);
			void addPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex , std::string _resname);
			std::string getChainId(unsigned int _chainIndex) const;
			std::string getResidueNumber(unsigned int _chainIndex, unsigned int _positionIndex) const;
			std::string getPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex) const;

			void setName(std::string _name);
			std::string getName();

			std::vector<std::vector<std::vector<std::string> > > getSequence() const;

			void setPDBNamesFlag(bool _pdbNamesFlag);
			bool getPDBNamesFlag();

			void setReferenceSequence(std::string _refSeq, std::string _refName, int _startRefResidueNumber, int _equivalentRefRes, int _equivalentPolyRes);

			std::string getReferenceHeader();
			std::string toString() const;
			std::string getFormattedString() { return sequenceFormatted;}

			static std::string toThreeLetterCode(AtomPointerVector &_av,std::string _residueDefiningAtomType="CA");
			static std::string toOneLetterCode(AtomPointerVector &_av,std::string _residueDefiningAtomType="CA");

			

		private:
			void setup(std::string _sequence);
			void copy(const PolymerSequence & _seq);

			void parseString(std::string _sequence);

			std::vector<std::vector<std::vector<std::string> > > sequence; // std::vector of std::vector of std::vector of sting for the chain/position/identity dimensions
			std::vector<std::vector<std::string> > residueNumbers;
			std::vector<std::string> chainIds;

			std::string sequenceName;
			std::string sequenceFormatted;

			// Store a reference sequence . It would be nice to have this a PolymerSequence object..
			bool refSeqFlag;
			std::string refSequence;
			std::string refName;

			int refStartResNum;
			int refEquilResNum;
			int seqEquilResNum;

			bool pdbNamesFlag;

	};


	inline unsigned int PolymerSequence::size() const {return sequence.size();}
	inline unsigned int PolymerSequence::chainSize(unsigned int _chainIndex) const {return sequence[_chainIndex].size();}
	inline unsigned int PolymerSequence::positionSize(unsigned int _chainIndex, unsigned int _positionIndex) const {return sequence[_chainIndex][_positionIndex].size();}
	inline void PolymerSequence::setChainId(unsigned int _chainIndex, std::string _chainId) {chainIds[_chainIndex] = _chainId;}
	inline void PolymerSequence::setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, int _resnum) {residueNumbers[_chainIndex][_positionIndex] = MslTools::intToString(_resnum);}
	inline void PolymerSequence::setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, std::string _resnum) {residueNumbers[_chainIndex][_positionIndex] = _resnum;}
	inline void PolymerSequence::setPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex, std::string _resname) {sequence[_chainIndex][_positionIndex][_identityIndex] = _resname;}
	inline void PolymerSequence::addPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex , std::string _resname) {sequence[_chainIndex][_positionIndex].push_back(_resname); }
	inline std::string PolymerSequence::getChainId(unsigned int _chainIndex) const {return chainIds[_chainIndex];}
	inline std::string PolymerSequence::getResidueNumber(unsigned int _chainIndex, unsigned int _positionIndex) const {return residueNumbers[_chainIndex][_positionIndex];}
	inline std::string PolymerSequence::getPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex) const {return sequence[_chainIndex][_positionIndex][_identityIndex];}
	inline std::vector<std::vector<std::vector<std::string> > > PolymerSequence::getSequence() const {return sequence;}
	inline void PolymerSequence::setName(std::string _name) { sequenceName = _name; }
	inline std::string PolymerSequence::getName()           { return sequenceName; }
	inline bool PolymerSequence::readSequenceFromPDB(std::string _pdbfile) {
		PDBReader reader(_pdbfile);
		if (!reader.read()) {
			return false;
		}
		setSequence(reader.getAtomPointers());
		return true;
	}


	inline void PolymerSequence::setPDBNamesFlag(bool _flag) { pdbNamesFlag = _flag;}
	inline bool PolymerSequence::getPDBNamesFlag() { return pdbNamesFlag;}
}


#endif

