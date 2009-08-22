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

#ifndef TWOBODYDISTANCEDEPENDENTPOTENTIALTABLE_H
#define TWOBODYDISTANCEDEPENDENTPOTENTIALTABLE_H

#include "PotentialTable.h"

#include "TBDReader.h"

class TwoBodyDistanceDependentPotentialTable : public PotentialTable {

	public:
		TwoBodyDistanceDependentPotentialTable();

		void readPotentialTable(string _fileName, string _potentialName);
		void readPotentialTable(string _fileName);

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

		double getPotential(string _body1, string _body2, int _dist);
	

		void addPotential(string _name1, string _name2, int _distBin, double _value);
		
		

	private:

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

		vector<distanceBin> distBins; // Read in from potential file

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

#endif
