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

#include "ChiStatistics.h"

using namespace MSL;
using namespace std;


ChiStatistics::ChiStatistics(){
	SysEnv env;
	dofReader.read(env.getEnv("MSL_PDB_2.3_DOF"));
}

bool ChiStatistics::read(string _dofFile){
	return dofReader.read(_dofFile);
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
	dofReader = _chiStat.dofReader;
}

int ChiStatistics::getNumberChis(Residue &_n){
	vector< vector<string> > chis = dofReader.getChiAtoms(_n.getResidueName());
	return chis.size();
}

bool ChiStatistics::atomsExist(Residue &_n, int _chiNumber){

	// Change from 1 based chi numbering to zero based..
	_chiNumber -= 1;

	vector< vector<string> > chis = dofReader.getChiAtoms(_n.getResidueName());
	if  (chis.size() < 1){
		cerr << "ERROR 4232 residue "<<_n.getResidueName()<<" is not found in chi table in ChiStatistics."<<endl;
		return false;
	}
	if (chis.size() < _chiNumber){
		cerr << "ERROR 4233 residue "<<_n.getResidueName()<<" does not have chi"<<_chiNumber << " in chi table in ChiStatistics."<<endl;
		return false;
	}


	if (!  (  _n.atomExists( chis[_chiNumber][0] ) &&
		  _n.atomExists( chis[_chiNumber][1] ) &&
		  _n.atomExists( chis[_chiNumber][2] ) &&
		  _n.atomExists( chis[_chiNumber][3] ))){
		cerr << "ERROR 4234 residue "<<_n.getResidueName()<<" does not have an atom it needs"<<endl;
		return false;
	}
	return true;
}
double ChiStatistics::getChi(Residue &_n, int _chiNumber,bool _angleInRadians){

	if(!atomsExist(_n,_chiNumber)) {
		return MslTools::doubleMax;
	}

	// Change from 1 based chi numbering to zero based..
	_chiNumber -= 1;

	vector< vector<string> > chis = dofReader.getChiAtoms(_n.getResidueName());
		
	if (_angleInRadians) {
		return _n(chis[_chiNumber][0]).dihedralRadians(_n(chis[_chiNumber][1]), _n(chis[_chiNumber][2]),_n(chis[_chiNumber][3]));
	}

	return _n(chis[_chiNumber][0]).dihedral(_n(chis[_chiNumber][1]), _n(chis[_chiNumber][2]),_n(chis[_chiNumber][3]));

}

vector<double> ChiStatistics::getChis(Residue &_n,bool _angleInRadians) {
	int nChis = getNumberChis(_n);
	vector<double> chis;
	for(unsigned i = 1 ; i <= nChis; i++) {
		chis.push_back(getChi(_n,i,_angleInRadians));
	}
	return chis;
}



