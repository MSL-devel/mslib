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


#include "HydrogenBondBuilder.h"

using namespace MSL;
using namespace std;


HydrogenBondBuilder::HydrogenBondBuilder() {
	setup();
}

HydrogenBondBuilder::HydrogenBondBuilder(System & _system, string _scwrl4ParameterFile) {
	setup();
	pSystem = &_system;
	readParameters(_scwrl4ParameterFile);
}

HydrogenBondBuilder::HydrogenBondBuilder( HydrogenBondBuilder & _sysBuild) {
	setup();
	copy(_sysBuild);
}

HydrogenBondBuilder::~HydrogenBondBuilder() {
	deletePointers();
}

void HydrogenBondBuilder::operator=( HydrogenBondBuilder & _sysBuild) {
	copy(_sysBuild);
}


void HydrogenBondBuilder::setup() {
	pSystem = NULL;
}

void HydrogenBondBuilder::reset() {
	acceptors.clear();
	donors.clear();
	acceptorData.clear();
	donorData.clear();
	donorAtom.clear();
	lonePairAngleAtom.clear();
	lonePairDihedralAtom.clear();

	deletePointers();
	setup();
}

void HydrogenBondBuilder::copy( HydrogenBondBuilder & _sysBuild) {
	reset();
	pSystem = _sysBuild.pSystem;
	acceptors = _sysBuild.acceptors;
	donors = _sysBuild.donors;
	acceptorData = _sysBuild.acceptorData;
	donorData = _sysBuild.donorData;
	donorAtom = _sysBuild.donorAtom;
	lonePairAngleAtom = _sysBuild.lonePairAngleAtom;
	lonePairDihedralAtom = _sysBuild.lonePairDihedralAtom;
}

void HydrogenBondBuilder::deletePointers() {
}

bool HydrogenBondBuilder::readParameters(std::string _parameterFile) {
	Reader pParReader(_parameterFile);
	pParReader.open();
	if (!(pParReader.is_open())) {
		cerr << "WARNING: Unable to open " << _parameterFile << endl;
		return false;
	} 
		
	const vector<string> lines = pParReader.getAllLines();
	bool foundResidue = false;
	string resName = "";

	for(int i = 0; i < lines.size(); i++) {
	//	cout << "UUU Line " << lines[i] << endl;
		MslTools::uncomment(lines[i],"!");
		
		vector <string> tokens = MslTools::tokenize(lines[i]);

		if(tokens.size() < 2) {
			continue;
		}
		if(MslTools::toUpper(tokens[0]).substr(0,4) == "RESI") {
			foundResidue = true;
			resName = MslTools::toUpper(tokens[1]);
		} else if (MslTools::toUpper(tokens[0]) == "DONOR" && foundResidue) {
			donorAtom[resName + " " + tokens[1]] = tokens[2];
			for(int j = 3; j < tokens.size(); j++) {
				donorData[resName + " " + tokens[1] ].push_back(MslTools::toDouble(tokens[j]));
			}
		} else if (MslTools::toUpper(tokens[0]) == "ACCEPTOR" && foundResidue) {
			lonePairAngleAtom[resName + " " + tokens[1]] = tokens[2];
			lonePairDihedralAtom[resName + " " + tokens[1]] = tokens[3];
			for(int j = 4; j < tokens.size(); j++) {
				acceptorData[resName + " " + tokens[1] ].push_back(MslTools::toDouble(tokens[j]));
			}

		} else  {
			if(foundResidue) {
				cerr << "WARNING: 134464 Unexpected Line Type " <<  tokens[0] << " in " << _parameterFile << " line " << i <<  endl;
			}
		}
	}

	pParReader.close();
	return true;
}

void HydrogenBondBuilder::printParameters () {
	cout << "Donors " << endl;
	
	for(map<string,vector<double> >::iterator it = donorData.begin(); it != donorData.end(); it++) {
			cout << it->first << " ";
			for(vector<double>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++) {
				cout << *it1 << " " ;
			}
			cout << endl;
	}
	cout << "Acceptors " << endl;
	for(map<string,vector<double> >::iterator it = acceptorData.begin(); it != acceptorData.end(); it++) {
			cout << it->first << " ";
			for(vector<double>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++) {
				cout << *it1 << " " ;
			}
			cout << endl;
	}

}

void HydrogenBondBuilder::collectDonorsAndAcceptors() {
	acceptors.clear();
	donors.clear();
	
	
	AtomPointerVector& atoms = pSystem->getAllAtomPointers();
	for(int i = 0; i < atoms.size(); i++) {
		string atomId = atoms[i]->getResidueName() + " " + atoms[i]->getName();
		if(acceptorData.find(atomId) != acceptorData.end() ) {
			Residue* res = atoms[i]->getParentResidue();
			Atom* angleAtom = NULL;
			if(res->atomExists(lonePairAngleAtom[atomId])) {
				angleAtom = &(res->getLastFoundAtom());
			}
			Atom* dihedralAtom = NULL;
			if(res->atomExists(lonePairDihedralAtom[atomId])) {
				dihedralAtom = &(res->getLastFoundAtom());
			}
			if(atoms[i] && angleAtom && dihedralAtom) {
				acceptors.push_back(vector<Atom*>());
				acceptors.back().push_back(atoms[i]);
				acceptors.back().push_back(angleAtom);
				acceptors.back().push_back(dihedralAtom);
			}

			continue;
		}
		if(donorData.find(atomId) != donorData.end()) {
			Residue* res = atoms[i]->getParentResidue();
			Atom* dAtom = NULL;
			if(res->atomExists(donorAtom[atomId])) {
				dAtom = &(res->getLastFoundAtom());
			}
			if(atoms[i] && dAtom) {
				donors.push_back(vector<Atom*>());
				donors.back().push_back(atoms[i]);
				donors.back().push_back(dAtom);
			}
			continue;
		}
	}
}

bool HydrogenBondBuilder::buildInteractions(double _cutoff) {
	collectDonorsAndAcceptors();
	return update(_cutoff);
}
bool HydrogenBondBuilder::update(double _cutoff) {
//	cout << "Num Acceptors " << acceptors.size() << endl;
//	cout << "Num Donors " << donors.size() << endl;
	EnergySet* ESet = pSystem->getEnergySet();
	// delete all existing scwrl4HBondinteractions
	ESet->eraseTerm("SCWRL4_HBOND");
	//cout << "UUUU Size: " << atoms.size() << endl;

	// create a map of acceptors and 

	for(int i = 0; i < acceptors.size(); i++) {
		Atom* acceptor = acceptors[i][0];
		string accPosId = acceptor->getPositionId();
		string accResName =  acceptor->getResidueName();
		string acceptorId = accResName + " " + acceptor->getName();
		vector<double>& accData = acceptorData[acceptorId];

		for(int j = 0; j < donors.size(); j++) {
			Atom* hydrogen = donors[j][0];
			string donPosId = hydrogen->getPositionId();
			string donResName =  hydrogen->getResidueName();
			if(accPosId == donPosId) { 
				// If acceptor and donor belong to the same position ,
				// make sure they are indeed distinct 
				// make sure that we build interactions only if the identities are same
				if(accResName != donResName) {
					// atoms from different identity - skip
					continue;	
				}
				if(donors[j][1]->getName() == acceptor->getName()) {
					// both atoms are from the same identity
					// make sure they are indeed different
					continue;
				}
			}
			string hydrogenId = hydrogen->getResidueName() + " " + hydrogen->getName(); 
			
			if(_cutoff > 0) {
				// a cutoff has been specified
				if(hydrogen->distance(*acceptor) > _cutoff) {
					continue;
				}
			}

			vector<double>& donData = donorData[hydrogenId];

			Scwrl4HBondInteraction* interaction = new Scwrl4HBondInteraction(*hydrogen,*donors[j][1],*acceptor,*acceptors[i][1],*acceptors[i][2], // atoms involved
				accData[0],(accData[1] * M_PI)/180.0,(accData[2] * M_PI)/180.0,(accData[3] * M_PI)/180.0, // acceptor data - angles in radians
				donData[0],donData[1],donData[2],(donData[3] * M_PI)/180.0,(donData[4] * M_PI)/180.0); // donor data - angles in radians
			ESet->addInteraction(interaction);
		}
	}
	return true;
}

