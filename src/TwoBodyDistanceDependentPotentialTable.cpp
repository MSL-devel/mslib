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

#include "TwoBodyDistanceDependentPotentialTable.h"
#include "MslExceptions.h"


void TwoBodyDistanceDependentPotentialTable::readPotentialTable(string _fileName, string _potentialName){
	fileName      = _fileName;
	potentialName = _potentialName;

	
	reader.open(_fileName);
	reader.read(this);
	reader.close();


}

void TwoBodyDistanceDependentPotentialTable::readPotentialTable(string _fileName){
	fileName      = _fileName;
	
	reader.open(_fileName);
	reader.read(this);
	reader.close();


	cout << "Potential: "<<potentialName<<endl;
	cout << "Residue Skipping Number: "<<residueSkippingNumber<<endl;
	cout << "Bins: "<<endl;
	for (uint i = 0 ; i < distBins.size();i++){
		cout << "\t"<<distBins[i].index<<" [ "<<distBins[i].startDistance<<" , "<<distBins[i].endDistance<<" ] "<<endl;
	}

}


void TwoBodyDistanceDependentPotentialTable::addPotential(string _name1, string _name2, int _distBin, double _value){

	if (_distBin+1 >=  distBins.size() ){
		stringstream ss;
		ss << "TwoBodyDistanceDependentPotentialTable::addPotential() , asking for a distance bin that is not specified: "<<_distBin<<" size of vector : "<<distBins.size();
		throw MslSizeException(ss.str());
	}

	
	stringstream key;
	key << _name1 <<":"<<_name2<<":"<<_distBin;

	PotentialTable::addPotential(key.str(),_value);
}

double TwoBodyDistanceDependentPotentialTable::getPotential(string _name1, string _name2, int _distBin){

	if (distBins.size() < _distBin){
		stringstream ss;
		ss << "TwoBodyDistanceDependentPotentialTable::getPotential() , asking for a distance bin that is not specified: "<<_distBin<<" size of vector : "<<distBins.size();
		throw MslSizeException(ss.str());
	}

	
	stringstream key;
	key << _name1 <<":"<<_name2<<":"<<_distBin;

	return PotentialTable::getPotential(key.str());
}
