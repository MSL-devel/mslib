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


#include "BaselineEnergyBuilder.h"

using namespace MSL;
using namespace std;


BaselineEnergyBuilder::BaselineEnergyBuilder() {
	setup();
}

BaselineEnergyBuilder::BaselineEnergyBuilder(System & _system, string _baselineFile) {
	setup();
	pSystem = &_system;
	readParameters(_baselineFile);
}
BaselineEnergyBuilder::BaselineEnergyBuilder(System & _system, map<string,map<string,double> >& _pBaselineEnergies) {
	setup();
	pSystem = &_system;
	pBaselineEnergies = _pBaselineEnergies;
}
BaselineEnergyBuilder::BaselineEnergyBuilder( BaselineEnergyBuilder & _builder) {
	setup();
	copy(_builder);
}

BaselineEnergyBuilder::~BaselineEnergyBuilder() {
}

void BaselineEnergyBuilder::operator=( BaselineEnergyBuilder & _builder) {
	copy(_builder);
}


void BaselineEnergyBuilder::setup() {
	pSystem = NULL;
	pBaselineEnergies.clear();
}

void BaselineEnergyBuilder::reset() {
	setup();
}

void BaselineEnergyBuilder::copy( BaselineEnergyBuilder & _builder) {
	reset();
	pSystem = _builder.pSystem;
	pBaselineEnergies = _builder.pBaselineEnergies;
}

bool BaselineEnergyBuilder::readParameters(std::string _baselineFile) {
	Reader pBaselineReader(_baselineFile);
	pBaselineReader.open();
	if (!(pBaselineReader.is_open())) {
		cerr << "WARNING: Unable to open " << _baselineFile << endl;
		return false;
	} 

	pBaselineEnergies.clear();
		
	vector<string> lines = pBaselineReader.getAllLines();

	for(int i = 0; i < lines.size(); i++) {
	//	cout << "UUU Line " << lines[i] << endl;

		lines[i] = MslTools::uncomment(lines[i],"!");
		vector <string> tokens = MslTools::tokenize(lines[i]," \t");
		if(tokens.size() < 1) {
			continue;
		}
		if(tokens.size() != 3 ) { 
			cerr << "WARNING 14353 : Line \"" << lines[i] << "\" is not in FORMAT: <ResName(string) or IdentityId> <AtomName(string)> <Energy(double)> " << endl; 
			continue;
		}
		pBaselineEnergies[MslTools::toUpper(tokens[0])][MslTools::toUpper(tokens[1])] = MslTools::toDouble(tokens[2]);
	}

	pBaselineReader.close();
	return true;
}

void BaselineEnergyBuilder::printParameters () {
	for(map<string,map<string,double> >::iterator it = pBaselineEnergies.begin(); it != pBaselineEnergies.end(); it++) {
		for(map<string,double>::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
			cout << it->first << " " << it1->first << " " << it1->second << endl;
		}
	}
}

bool BaselineEnergyBuilder::buildInteractions() {
	vector<Position*>& positions = pSystem->getPositions();
	EnergySet* ESet = pSystem->getEnergySet();
	//cout << "UUUU Size: " << atoms.size() << endl;

	for(vector<Position*>::iterator p = positions.begin(); p != positions.end(); p++) {
		for(unsigned int i = 0; i < (*p)->identitySize(); i++) {
			Residue& res = (*p)->getIdentity(i);
			string identityId = res.getIdentityId();
			// the baseId could be a resName or a identityId
			string baseId = res.getResidueName();
			if(pBaselineEnergies.find(identityId) != pBaselineEnergies.end()) {
				baseId = identityId;
			} else if(pBaselineEnergies.find(baseId) == pBaselineEnergies.end()) {
				continue;
			}
			for(map<string,double>::iterator it = pBaselineEnergies[baseId].begin(); it != pBaselineEnergies[baseId].end(); it++) {
				Atom *a = NULL;
				if(res.atomExists(it->first)) {
					a = &res.getLastFoundAtom();
				} else {
					continue;
				}
				ESet->addInteraction(new BaselineInteraction(*a,it->second));
			}
		}
	}
	return true;
}

