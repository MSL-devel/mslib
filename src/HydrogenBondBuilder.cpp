/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
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
	deletePointers();
	setup();
}

void HydrogenBondBuilder::copy( HydrogenBondBuilder & _sysBuild) {
	reset();
	pSystem = _sysBuild.pSystem;
	// copy donors
	for(map<string,map<string,string> >::iterator it = _sysBuild.donorData.begin(); it != _sysBuild.donorData.end(); it++) {
		for(map<string,string>::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++ ) {
			donorData[it->first][it1->first] = it1->second; 
		}
	}
	// copy acceptors
	for(map<string,map<string,vector<string> > >::iterator it = _sysBuild.acceptorData.begin(); it != _sysBuild.acceptorData.end(); it++) {
		for(map<string,vector<string> >::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
			for(int j = 0; j < _sysBuild.acceptorData[it->first][it1->first].size(); j++) {
				acceptorData[it->first][it1->first].push_back(_sysBuild.acceptorData[it->first][it1->first][j]);
			}	
		}
		
	}
}

void HydrogenBondBuilder::deletePointers() {
	for(map<string,map<string,string> >::iterator it = donorData.begin(); it != donorData.end(); it++) {
		donorData[it->first].clear();
	}
	donorData.clear();
	for(map<string,map<string,vector<string> > >::iterator it = acceptorData.begin(); it != acceptorData.end(); it++) {
		for(map<string,vector<string> >::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
			acceptorData[it->first][it1->first].clear();
		}
		acceptorData[it->first].clear();
	}
	
	acceptorData.clear();
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

		vector <string> tokens = MslTools::tokenize(lines[i]);
		if(tokens.size() < 2 || lines[i].substr(0,1) == "!") {
			continue;
		}
		if(MslTools::toUpper(tokens[0]) == "RESI") {
			foundResidue = true;
			resName = MslTools::toUpper(tokens[1]);
		} else if (MslTools::toUpper(tokens[0]) == "DONOR" && foundResidue) {
			donorData[resName][tokens[1]] = tokens[2];
		} else if (MslTools::toUpper(tokens[0]) == "ACCEPTOR" && foundResidue) {
			for(int j = 2; j < tokens.size(); j++) {
				acceptorData[resName][tokens[1]].push_back(tokens[j]);
			}

		} else {
			if(foundResidue) {
				cerr << "WARNING: 134464 Unexpected Line Type " <<  tokens[0] << " in " << _parameterFile << " line " << i <<  endl;
			} else {
				cerr << "WARNING: 145834 No RESI line in SCWRL4 Parameter File " << endl;
			}
		}
	}

	pParReader.close();
	return true;
}

void HydrogenBondBuilder::printParameters () {
	cout << "Donors " << endl;
	
	for(map<string,map<string,string> >::iterator it = donorData.begin(); it != donorData.end(); it++) {
		for(map<string,string>::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
			cout << it->first << " " << it1->first << " " << it1->second << endl;
		}
	}
	cout << "Acceptors " << endl;
	for(map<string,map<string,vector<string> > >::iterator it = acceptorData.begin(); it != acceptorData.end(); it++) {
		for(map<string,vector<string> >::iterator it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
			cout << it->first << " " << it1->first << " " ;
			for(vector<string>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				cout << *it2 << " " ;
			}
			cout << endl;
		}
	}
	

}

bool HydrogenBondBuilder::buildInteractions(double _cutoff) {
	AtomPointerVector& atoms = pSystem->getAllAtomPointers();
	EnergySet* ESet = pSystem->getEnergySet();
	//cout << "UUUU Size: " << atoms.size() << endl;

	for(int i = 0; i < atoms.size() -1; i++) {
		string res1Name = atoms[i]->getResidueName();
		bool acceptor = false;
		string atom1Name = atoms[i]->getName();
		if(acceptorData.find(res1Name) != acceptorData.end() && acceptorData[res1Name].find(atom1Name) != acceptorData[res1Name].end()) {
				acceptor = true;
		} else if (donorData.find(res1Name) != donorData.end() && donorData[res1Name].find(atom1Name) == donorData[res1Name].end()) {
				continue;
		} else {
			continue;
		}
		
		for(int j = i + 1; j < atoms.size(); j++) {
			// atoms[i] is either a donor or acceptor if control reaches here.
			// if atom[i] is an acceptor, atom[j] needs to be a donor or vice-versa
			string res2Name = atoms[j]->getResidueName(); 
			string atom2Name = atoms[j]->getName(); 
			if(acceptor) {
				if (donorData.find(res2Name) != donorData.end() && donorData[res2Name].find(atom2Name) == donorData[res2Name].end()) {
						continue;
				}	
			} else {
				if(acceptorData.find(res2Name) != acceptorData.end() && acceptorData[res2Name].find(atom2Name) == acceptorData[res2Name].end()) {
					continue;	
				}
			} 
			// we have one acceptor and one donor atom now
			Atom* donor_1 = NULL; // actual donor
			Atom* donor_2 = NULL; // bonded to donor
			Atom* acceptor_1 = NULL; // actual acceptor
			Atom* acceptor_2 = NULL; // bonded to acceptor
			Atom* acceptor_3 = NULL; // bonded to acceptor_2

			Residue *res1 = atoms[i]->getParentResidue();
			Residue *res2 = atoms[j]->getParentResidue();
			vector<double> data;

			if(acceptor) {
				// atom[i] is an acceptor and atom[j] is a donor
				donor_1 = atoms[j];
				acceptor_1 = atoms[i];
				if(_cutoff > 0) {
					if(donor_1->distance(*acceptor_1) > _cutoff) {
						continue;
					}
				}

				if(res1->atomExists(acceptorData[res1Name][atom1Name][0])) {
					acceptor_2 = &(res1->getLastFoundAtom());
				} else {
					continue;
				}

				if(res1->atomExists(acceptorData[res1Name][atom1Name][1])) {
					acceptor_3 = &(res1->getLastFoundAtom()); 
				} else {
					continue;
				}

				if(res2->atomExists(donorData[res2Name][atom2Name])) {
					donor_2 = &(res2->getLastFoundAtom());
				} else {
					continue;
				}

				vector<string>& accData = acceptorData[res1Name][atom1Name];
				data.push_back(MslTools::toDouble(accData[2]));
				for(int l = 3; l < accData.size(); l++) {
					data.push_back(MslTools::toDouble(accData[l]) * M_PI /180.0);
				}
			} else {
				donor_1 = atoms[i];
				acceptor_1 = atoms[j];
				if(_cutoff > 0) {
					if(donor_1->distance(*acceptor_1) > _cutoff) {
						continue;
					}
				}

				if(res2->atomExists(acceptorData[res2Name][atom2Name][0])) {
					acceptor_2 = &(res2->getLastFoundAtom());
				} else {
					continue;
				}

				if(res2->atomExists(acceptorData[res2Name][atom2Name][1])) {
					acceptor_3 = &(res2->getLastFoundAtom());
				} else {
					continue;
				}

				if(res1->atomExists(donorData[res1Name][atom1Name])) {
					donor_2 = &(res1->getLastFoundAtom());
				} else {
					continue;
				}

				vector<string>& accData = acceptorData[res1Name][atom1Name];
				data.push_back(MslTools::toDouble(accData[2]));
				for(int l = 3; l < accData.size(); l++) {
					data.push_back(MslTools::toDouble(accData[l]) * M_PI /180.0);
				}

			}

			ESet->addInteraction(new Scwrl4HBondInteraction(*donor_1,*donor_2,*acceptor_1,*acceptor_2,*acceptor_3,data[0],data[1],data[2],data[3]));

		}
	}
	return true;
}

