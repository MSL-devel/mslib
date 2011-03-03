/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2011 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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

#include "ChiStatistics.h"

using namespace MSL;
using namespace std;


ChiStatistics::ChiStatistics(){

}


ChiStatistics::~ChiStatistics(){
}

ChiStatistics::ChiStatistics(const ChiStatistics &_chiStat){
	copy(_chiStat);
}

void ChiStatistics::operator=(const ChiStatistics &_chiStat){
	copy(_chiStat);
}

void ChiStatistics::copy(const ChiStatistics &_chiStat){
	// Emtpy for now..
}

int ChiStatistics::getNumberChis(Residue &_n){

	map<string, vector< vector<string> > >::iterator it;
	it = pTop.getChis().find(_n.getResidueName());
	if  (it == pTop.getChis().end()){
		//cerr << "ERROR 4232 residue "<<_n.getResidueName()<<" is not found in chi table in ChiStatistics."<<endl;
		return -1;
	}

	return it->second.size();
}

bool ChiStatistics::atomsExist(Residue &_n, int _chiNumber){

	// Change from 1 based chi numbering to zero based..
	_chiNumber -= 1;

	map<string, vector< vector<string> > >::iterator it;
	it = pTop.getChis().find(_n.getResidueName());
	if  (it == pTop.getChis().end()){
		cerr << "ERROR 4232 residue "<<_n.getResidueName()<<" is not found in chi table in ChiStatistics."<<endl;
		return false;
	}
	if (it->second.size() < _chiNumber){
		cerr << "ERROR 4233 residue "<<_n.getResidueName()<<" does not have chi angle ("<<_chiNumber<< ") in chi table in ChiStatistics."<<endl;
		return false;
	}


	if (!  (  _n.atomExists( (it->second)[_chiNumber][0] ) &&
		  _n.atomExists( (it->second)[_chiNumber][1] ) &&
		  _n.atomExists( (it->second)[_chiNumber][2] ) &&
		  _n.atomExists( (it->second)[_chiNumber][3] ))){
		cerr << "ERROR 4234 residue "<<_n.getResidueName()<<" does not have an atom it needs"<<endl;
		return false;
	}
	return true;
}
double ChiStatistics::getChi(Residue &_n, int _chiNumber,bool _angleInRadians){

	// Change from 1 based chi numbering to zero based..
	_chiNumber -= 1;

	map<string, vector< vector<string> > >::iterator it;
	it = pTop.getChis().find(_n.getResidueName());
	if  (it == pTop.getChis().end()){
		cerr << "ERROR 4232 residue "<<_n.getResidueName()<<" is not found in chi table in ChiStatistics."<<endl;
		return MslTools::doubleMax;
	}
	if (it->second.size() < _chiNumber){
		cerr << "ERROR 4233 residue "<<_n.getResidueName()<<" does not have chi angle ("<<_chiNumber<< ") in chi table in ChiStatistics."<<endl;
		return MslTools::doubleMax;
	}


	if (!  (  _n.atomExists( (it->second)[_chiNumber][0] ) &&
		  _n.atomExists( (it->second)[_chiNumber][1] ) &&
		  _n.atomExists( (it->second)[_chiNumber][2] ) &&
		  _n.atomExists( (it->second)[_chiNumber][3] ))){
		cerr << "ERROR 4234 One of the atoms ("<<(it->second)[_chiNumber][0]<<","<<(it->second)[_chiNumber][1]<<","<<(it->second)[_chiNumber][2]<<","<<(it->second)[_chiNumber][3]<<") doesn't exist in residue "<<_n.getResidueName()<<" chi angle "<<_chiNumber<<endl;
		return MslTools::doubleMax;
	}
	
	if (_angleInRadians) {
		return _n((it->second)[_chiNumber][0]).dihedralRadians(_n((it->second)[_chiNumber][1]), _n((it->second)[_chiNumber][2]),_n((it->second)[_chiNumber][3]));
	}

	return _n((it->second)[_chiNumber][0]).dihedral(_n((it->second)[_chiNumber][1]), _n((it->second)[_chiNumber][2]),_n((it->second)[_chiNumber][3]));

}




