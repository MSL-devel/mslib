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

#ifndef TWOBODYDISTANCEDEPENDENTPOTENTIALTABLE_H
#define TWOBODYDISTANCEDEPENDENTPOTENTIALTABLE_H

#include "PotentialTable.h"

#include "TBDReader.h"
#include "AtomPointerVector.h"
#include "System.h"

// Forward declaration
//class System;

namespace MSL { 
class TwoBodyDistanceDependentPotentialTable : public PotentialTable {

	public:
		TwoBodyDistanceDependentPotentialTable();

		void readPotentialTable(std::string _fileName, std::string _potentialName);
		void readPotentialTable(std::string _fileName);

               	int getResidueSkippingNumber();
               	void setResidueSkippingNumber(int _skipNum);

		double getMinDistCutoff();
		double getValueBelowCutoff();
		void setMinDistCutoffAndValue(double _minDist, double _value);

		double getMaxDistCutoff();
		double getValueAboveCutoff();
		void setMaxDistCutoffAndValue(double _maxDist, double _value=0);

		void addBin(double startDistance, double endDistance); // No checking for overlapping bins!
		int getBin(double _distance);

		double getPotential(std::string _body1, std::string _body2, int _dist);
		void addPotential(std::string _name1, std::string _name2, int _distBin, double _value);

		double getEnergyBetweenResidues();

		double calculateSelfEnergy(System &_sys, int _position, int _rotamer);
		double calculateTemplateEnergy(System &_sys, int _position, int _rotamer, bool _calcAllForFixed=false, bool _countLocalSCBB=false);
		double calculatePairEnergy(System &_sys, int _position1, int _rotamer1, int _position2, int _rotamer2, bool _countLocalSCBB=false);

		double calculateBackgroundEnergy(System &_sys, int _position, int _rotamer, bool _countLocalSCBB=false);
		double calculateSurroundingEnergy(System &_sys, int _position, int _rotamer, std::vector< std::vector< std::vector< std::vector<double> > > > & rotamerInteractions, std::vector<uint> & currentAllRotamers, bool _countLocalSCBB=false);

	private:

		double calculatePairwiseNonBondedEnergy(System &_sys, AtomPointerVector &_a, AtomPointerVector &_b, bool _sameSet=false, bool _countLocalSCBB=false);

		bool isBackbone(std::string atomname);

		void setup(int _resSkipNum, double _minDistCutoff, double _valueBelowCutoff, double _maxDistCutoff);
		void copy(const TwoBodyDistanceDependentPotentialTable &_pairDisPot);

		int residueSkippingNumber;

		double minDistCutoff;
		double valueBelowCutoff;

		double maxDistCutoff;  
		double valueAboveCutoff; // Default value is 0 for two atoms beyond this distance.
		

		struct distanceBin {
			int index;
			double startDistance;
			double endDistance;
		};

		std::vector<distanceBin> distBins; // Read in from potential file

		TBDReader reader;

};

// INLINES
inline TwoBodyDistanceDependentPotentialTable::TwoBodyDistanceDependentPotentialTable() {}

inline int TwoBodyDistanceDependentPotentialTable::getResidueSkippingNumber() { return residueSkippingNumber; }
inline void TwoBodyDistanceDependentPotentialTable::setResidueSkippingNumber(int _resSkip) { residueSkippingNumber = _resSkip; }

inline double TwoBodyDistanceDependentPotentialTable::getMinDistCutoff() { return minDistCutoff; }
inline double TwoBodyDistanceDependentPotentialTable::getValueBelowCutoff() { return valueBelowCutoff; }

inline void TwoBodyDistanceDependentPotentialTable::setMinDistCutoffAndValue(double _minDist, double _value) { minDistCutoff = _minDist; valueBelowCutoff = _value; }

inline double TwoBodyDistanceDependentPotentialTable::getMaxDistCutoff() { return maxDistCutoff; }
inline double TwoBodyDistanceDependentPotentialTable::getValueAboveCutoff() { return valueAboveCutoff; }
inline void TwoBodyDistanceDependentPotentialTable::setMaxDistCutoffAndValue(double _maxDist, double _value) { maxDistCutoff = _maxDist; valueAboveCutoff = _value; }

inline void TwoBodyDistanceDependentPotentialTable::addBin(double _start, double _end){
	distanceBin db;
	db.index = distBins.size();
	db.startDistance = _start;
	db.endDistance = _end;

	distBins.push_back(db);
}

inline int TwoBodyDistanceDependentPotentialTable::getBin(double _dist){

	int bin = -1;
	for (uint i = 0; i < distBins.size();i++){
		if (_dist >= distBins[i].startDistance && _dist <= distBins[i].endDistance){
			bin = i;
			break;
		}
	}
	return bin;
}

}

#endif
