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

#ifndef POLYMERSEQUENCE_H
#define POLYMERSEQUENCE_H



// MSL INCLUDES
#include "AtomVector.h"
#include "System.h"
#include "MslTools.h"


// STL INCLUDES
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class PolymerSequence {
	public:
		PolymerSequence();
		PolymerSequence(string _sequence);
		PolymerSequence(System &_sys);
		PolymerSequence(const AtomVector &_atoms);
		PolymerSequence(const PolymerSequence & _seq);
		PolymerSequence(System &_sys, vector<pair<string,string> > &_addTerminalResidues); // add additional residues to polymer sequence
		PolymerSequence(System &_sys, map<string,map<int,int> > &variablePositionMap, vector<vector<string> > &_identitesAtVariablePositions);
		~PolymerSequence();

		void operator=(const PolymerSequence & _seq);
		friend ostream & operator<<(ostream &_os, const PolymerSequence & _seq)  {_os << _seq.toString(); return _os;};

		void setSequence(string _sequence);
		void setSequence(System &_sys);
		void setSequence(const AtomVector &_atoms);

		unsigned int size() const;
		unsigned int chainSize(unsigned int _chainIndex) const;
		unsigned int positionSize(unsigned int _chainIndex, unsigned int _positionIndex) const;
		void setChainId(unsigned int _chainIndex, string _chainId);
		void setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, int _resnum);
		void setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, string _resnum);
		void setPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex, string _resname);
		string getChainId(unsigned int _chainIndex) const;
		string getResidueNumber(unsigned int _chainIndex, unsigned int _positionIndex) const;
		string getPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex) const;

		void setName(string _name);
		string getName();

		vector<vector<vector<string> > > getSequence() const;


		void setReferenceSequence(string _refSeq, string _refName, int _startRefResidueNumber, int _equivalentRefRes, int _equivalentPolyRes);

		string getReferenceHeader();
		string toString() const;

		static string toThreeLetterCode(AtomVector &_av,string _residueDefiningAtomType="CA");
		static string toOneLetterCode(AtomVector &_av,string _residueDefiningAtomType="CA");

		

	private:
		void setup(string _sequence);
		void copy(const PolymerSequence & _seq);

		void parseString(string _sequence);

		vector<vector<vector<string> > > sequence; // vector of vector of vector of sting for the chain/position/identity dimensions
		vector<vector<string> > residueNumbers;
		vector<string> chainIds;

		string sequenceName;

		// Store a reference sequence . It would be nice to have this a PolymerSequence object..
		bool refSeqFlag;
		string refSequence;
		string refName;

		int refStartResNum;
		int refEquilResNum;
		int seqEquilResNum;


};


inline unsigned int PolymerSequence::size() const {return sequence.size();}
inline unsigned int PolymerSequence::chainSize(unsigned int _chainIndex) const {return sequence[_chainIndex].size();}
inline unsigned int PolymerSequence::positionSize(unsigned int _chainIndex, unsigned int _positionIndex) const {return sequence[_chainIndex][_positionIndex].size();}
inline void PolymerSequence::setChainId(unsigned int _chainIndex, string _chainId) {chainIds[_chainIndex] = _chainId;}
inline void PolymerSequence::setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, int _resnum) {residueNumbers[_chainIndex][_positionIndex] = MslTools::intToString(_resnum);}
inline void PolymerSequence::setPositionNumber(unsigned int _chainIndex, unsigned int _positionIndex, string _resnum) {residueNumbers[_chainIndex][_positionIndex] = _resnum;}
inline void PolymerSequence::setPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex, string _resname) {sequence[_chainIndex][_positionIndex][_identityIndex] = _resname;}
inline string PolymerSequence::getChainId(unsigned int _chainIndex) const {return chainIds[_chainIndex];}
inline string PolymerSequence::getResidueNumber(unsigned int _chainIndex, unsigned int _positionIndex) const {return residueNumbers[_chainIndex][_positionIndex];}
inline string PolymerSequence::getPositionIdentity(unsigned int _chainIndex, unsigned int _positionIndex, unsigned int _identityIndex) const {return sequence[_chainIndex][_positionIndex][_identityIndex];}
inline vector<vector<vector<string> > > PolymerSequence::getSequence() const {return sequence;}
inline void PolymerSequence::setName(string _name) { sequenceName = _name; }
inline string PolymerSequence::getName()           { return sequenceName; }
#endif

