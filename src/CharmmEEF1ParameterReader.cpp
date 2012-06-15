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


#include "CharmmEEF1ParameterReader.h"

using namespace MSL;
using namespace std;

const string CharmmEEF1ParameterReader::defaultSolvent = "WATER";


CharmmEEF1ParameterReader::CharmmEEF1ParameterReader() {
}

CharmmEEF1ParameterReader::CharmmEEF1ParameterReader(const string & _filename) {
	open(_filename);
}

CharmmEEF1ParameterReader::CharmmEEF1ParameterReader(const CharmmEEF1ParameterReader & _par) {
	copy(_par);
}

CharmmEEF1ParameterReader::~CharmmEEF1ParameterReader() {
}

void CharmmEEF1ParameterReader::operator=(const CharmmEEF1ParameterReader & _par) {
	copy(_par);
}

void CharmmEEF1ParameterReader::copy(const CharmmEEF1ParameterReader & _par) {
	reset();
	EEF1Map = _par.EEF1Map;
	//EEF1PairMap = _par.EEF1PairMap;
}

void CharmmEEF1ParameterReader::reset() {
	close();
	EEF1Map.clear();
	//EEF1PairMap.clear();
}

void CharmmEEF1ParameterReader::addEEF1(string _solvent, string _atomType, double _V, double _Gref, double _Gfree, double _Href, double _CPref, double _Sigw) {
	EEF1Map[_solvent][_atomType].push_back(_V    );
	EEF1Map[_solvent][_atomType].push_back(_Gref );
	EEF1Map[_solvent][_atomType].push_back(_Gfree);
	EEF1Map[_solvent][_atomType].push_back(_Href );
	EEF1Map[_solvent][_atomType].push_back(_CPref);
	EEF1Map[_solvent][_atomType].push_back(_Sigw );
}


bool CharmmEEF1ParameterReader::read() {
	if (!is_open()) {
		// file wasn't open
		return false;
	}
	try { 

		vector<string> lines;
		vector<vector<string> > splitFile;
		while (!endOfFileTest()){

			string line = Reader::getLine();
			if (line[0] == '*') {
				// skip title
				continue;
			}
			line = MslTools::toUpper(line);
			lines.push_back(line);
		}

		lines = MslTools::joinConnectedLines(lines, "-");
		lines = MslTools::uncomment(lines, "!");
		for (unsigned int i=0; i<lines.size(); i++) {
			vector<string> tokens = MslTools::tokenize(lines[i]," \t");  
			if (tokens.size() > 0) {
				splitFile.push_back(tokens);
			}
		}

		bool solventFound = false;
		string currentSolvent = "";
		for (vector<vector<string> >::iterator k=splitFile.begin(); k!=splitFile.end(); k++) {

			if (k->size() == 0) {
				continue;
			}
			//For now we process only Bonds and since Bonds is the first block break here
			if (!solventFound) {
				// section begins
				currentSolvent = (*k)[0];
				solventFound = true;
				continue;
			} else {
				if ((*k)[0].substr(0, 3) == "END") {
					// section ends
					currentSolvent = "";
					solventFound = false;
					continue;
				} else {
					if (k->size() < 7) {
						cerr << "WARNING 4761: invalid format of solvation input file in bool CharmmEEF1ParameterReader::read()" << endl;
					} else {
						double V     = MslTools::toDouble((*k)[1]);
						double Gref  = MslTools::toDouble((*k)[2]);
						double Gfree = MslTools::toDouble((*k)[3]);
						double Href  = MslTools::toDouble((*k)[4]);
						double CPref = MslTools::toDouble((*k)[5]);
						double Sigw  = MslTools::toDouble((*k)[6]);
						addEEF1(currentSolvent, (*k)[0], V, Gref, Gfree, Href, CPref, Sigw);
					}
				}
			}
		}

	} catch(...){
		cerr << "ERROR 8123 in bool CharmmEEF1ParameterReader::read()\n";
		exit(8123);
	}

	return true;

}

bool CharmmEEF1ParameterReader::EEF1Param(vector<double> & _params, string _type, string _solvent) const{

	if (_solvent == "") {
		_solvent = defaultSolvent;
	}
	map<string, map<string, vector<double> > >::const_iterator found1;
	map<string, vector<double> > ::const_iterator found2;

	if ((found1 = EEF1Map.find(_solvent)) != EEF1Map.end() && (found2 = (found1->second).find(_type)) != (found1->second).end()) {
		_params = found2->second;
		return true;
	} else {
		_params = vector<double>(8,0.0);
		//cerr << "WARNING 48230: EEF1 paramenter not found for types type (" << _solvent << "," << _type << ") in vector<double> CharmmEEF1ParameterReader::EEF1Param(string _type, string _solvent) const" << endl;
		//return(out);
		return false;
	}

}
/*
void CharmmEEF1ParameterReader::createEEF1ParamPairs(){


	map<string,vector<double> >::iterator it;
	map<string,vector<double> >::iterator it2;
	map<string,map<string,vector<double> > >::iterator pairIt;
	map<string,vector<double> >::iterator pairIt2;
	for (it = vdwParamMap.begin();it != vdwParamMap.end();it++){

		vector<double> p1 = it->second;
		vector<double> tmp(4,0.0);
		for (it2 = it;it2 != vdwParamMap.end();++it2){
			vector<double> p2 = it2->second;

			tmp[0] = sqrt( p1[0] * p2[0]);
			tmp[2] = sqrt( p1[2] * p2[2]);
			tmp[1] = p1[1] + p2[1];

			tmp.push_back(p1[3] + p2[3]);


			//cout << "Adding Pair: "<<it->first<<" "<<it2->first<<" "<<tmp[0]<< " vs [ "<<p1[0]<<", "<<p2[0]<<" ] ----> "<<sqrt( p1[0] * p2[0])<<endl;
			vdwParamPairMap[it->first][it2->first] = tmp;
		}

	}


}
vector<double> CharmmEEF1ParameterReader::EEF1ParamPair(string _type1, string _type2) const {

	map<string,map<string,vector<double> > >::const_iterator found1;
	map<string,vector<double> >::const_iterator found2;
	
	if ( (found1 = vdwParamPairMap.find(_type1)) != vdwParamPairMap.end()) {
		

		if ( (found2=found1->second.find(_type2)) != found1->second.end()){
			return found2->second;
		}
	}
		
	if ( (found1 = vdwParamPairMap.find(_type2)) != vdwParamPairMap.end()){
			
		if ( (found2 = found1->second.find(_type1)) != found1->second.end()){
			return found2->second;
		}

	}


	vector<double> out(4,0.0);
	cerr << "vdwPairParams not found for types " << _type1 << " "<<_type2<<endl;
	return(out);	
}
*/
