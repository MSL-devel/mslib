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

#include "CharmmParameterReader.h"

using namespace MSL;
using namespace std;


CharmmParameterReader::CharmmParameterReader() {
}

CharmmParameterReader::CharmmParameterReader(const string & _filename) {
	open(_filename);
}

CharmmParameterReader::CharmmParameterReader(const CharmmParameterReader & _par) {
	copy(_par);
}

CharmmParameterReader::~CharmmParameterReader() {
}

void CharmmParameterReader::operator=(const CharmmParameterReader & _par) {
	copy(_par);
}

void CharmmParameterReader::copy(const CharmmParameterReader & _par) {
	reset();
	bondParamMap = _par.bondParamMap;
	angleParamMap = _par.angleParamMap;
	dihedralParamMap = _par.dihedralParamMap;
	improperParamMap = _par.improperParamMap;
	vdwParamMap = _par.vdwParamMap;
	vdwParamPairMap = _par.vdwParamPairMap;
}

void CharmmParameterReader::reset() {
//	deletePointers();
	bondParamMap.clear();
	angleParamMap.clear();
	dihedralParamMap.clear();
	improperParamMap.clear();
	vdwParamMap.clear();
	vdwParamPairMap.clear();
}

void CharmmParameterReader::addBond(string _type1, string _type2, double _Kb, double _B0) {
	bondParamMap[_type1][_type2].clear();
	bondParamMap[_type1][_type2].push_back(_Kb);
	bondParamMap[_type1][_type2].push_back(_B0);
	if (_type1 != _type2) {
		bondParamMap[_type2][_type1].clear();
		bondParamMap[_type2][_type1].push_back(_Kb);
		bondParamMap[_type2][_type1].push_back(_B0);
	//	cout << "bondParamMap["<< _type2 << "][" << _type1 << "]" << bondParamMap[_type2][_type1][0] << bondParamMap[_type2][_type1][1] << endl;
	}
//	cout << "bondParamMap["<< _type1 << "][" << _type2 << "]" << bondParamMap[_type1][_type2][0] << bondParamMap[_type1][_type2][1] << endl;
}

void CharmmParameterReader::addAngle(string _type1, string _type2, string _type3, double _Ktheta, double _Theta0, double _Kub, double _S0) {
	angleParamMap[_type1][_type2][_type3].clear();
	angleParamMap[_type1][_type2][_type3].push_back(_Ktheta);
	angleParamMap[_type1][_type2][_type3].push_back(_Theta0);
	angleParamMap[_type1][_type2][_type3].push_back(_Kub);
	angleParamMap[_type1][_type2][_type3].push_back(_S0);

	if(_type1 != _type3) {
		angleParamMap[_type3][_type2][_type1].clear();
		angleParamMap[_type3][_type2][_type1].push_back(_Ktheta);
		angleParamMap[_type3][_type2][_type1].push_back(_Theta0);
		angleParamMap[_type3][_type2][_type1].push_back(_Kub);
		angleParamMap[_type3][_type2][_type1].push_back(_S0);
	}
}

void CharmmParameterReader::addDihedral(string _type1, string _type2, string _type3, string _type4, double _Kchi, double _N, double _Delta) {
	
	if (_type2 == "X" || _type3 == "X") {
		cerr << "type2 or type3 or both are wildcards in addDihedral (" << _type1 << "," << _type2 << "," << _type3 << "," << _type4 << ")" << endl;
	}
	
	vector<double> temp;
	
	temp.push_back(_Kchi);
	temp.push_back(_N);
	temp.push_back(_Delta);

	dihedralParamMap[_type4][_type3][_type2][_type1].push_back(temp);
	if(!(_type1 == _type4 && _type2 == _type3 )) {
	//	cout << "Adding Dihedral Map (" << _type1 << "," << _type2 << "," << _type3 << "," << _type4 << ")" << endl;
		dihedralParamMap[_type1][_type2][_type3][_type4].push_back(temp);
	}
}

void CharmmParameterReader::addImproper(string _type1, string _type2, string _type3, string _type4, double _Kpsi, double _Psi0) {

	improperParamMap[_type4][_type3][_type2][_type1].clear();
	improperParamMap[_type4][_type3][_type2][_type1].push_back(_Kpsi);
	improperParamMap[_type4][_type3][_type2][_type1].push_back(_Psi0);
	if(!(_type1 == _type4 && _type2 == _type3)) {
	//	cout << "Adding Improper Map (" << _type1 << "," << _type2 << "," << _type3 << "," << _type4 << ")" << endl;
		improperParamMap[_type1][_type2][_type3][_type4].clear();
		improperParamMap[_type1][_type2][_type3][_type4].push_back(_Kpsi);
		improperParamMap[_type1][_type2][_type3][_type4].push_back(_Psi0);
	}
}

void CharmmParameterReader::addVdw(string _type1, double _Eps, double _Rmin, double _Eps14, double _Rmin14) {
	vdwParamMap[_type1].clear();	
	vdwParamMap[_type1].push_back(_Eps);	
	vdwParamMap[_type1].push_back(_Rmin);	
	vdwParamMap[_type1].push_back(_Eps14);	
	vdwParamMap[_type1].push_back(_Rmin14);	
}


bool CharmmParameterReader::read() {
	vector<vector<string> > splitFile;
	enum BlockTypes { Bonds = 0, Angles = 1, Dihedrals = 2, Improper = 3, NonBonded = 4, HBond = 5} block;

	try { 
		vector<string> lines;
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

		for (vector<vector<string> >::iterator k=splitFile.begin(); k!=splitFile.end(); k++) {
			
			//For now we process only Bonds and since Bonds is the first block break here
						
			if ((*k)[0].substr(0, 3) == "END") {
				//cout << "End of file." << endl;
				break;
			}
			//If a block is beginning
			if ((*k)[0].substr(0, 4) == "BOND") {
				block = Bonds;
				//cout << "Reading Bonds" << endl;
				continue;
			}
			
			if ((*k)[0].substr(0, 4) == "ANGL" || (*k)[0].substr(0, 4) == "THET" ) {
				block = Angles;
				//cout << "Reading Angles" << endl;
				continue;
			}

			if ((*k)[0].substr(0, 4) == "DIHE"|| (*k)[0] == "PHI" ) {
				block = Dihedrals;
				//cout << "Reading Dihedrals" << endl;
				continue;
			}

			if ((*k)[0].substr(0, 4) == "IMPR"|| (*k)[0].substr(0, 4) == "IMPH" ) {
				block = Improper;
				//cout << "Improper" << endl;
				continue;
			}
			if ((*k)[0].substr(0, 4) == "NONB") {
				block = NonBonded;
				//cout << "NonBonded" << endl;
				continue;
			}
			if ((*k)[0].substr(0, 4) == "HBON") {
				block = HBond;
				//cout << "HBond" << endl;
				continue;
			}


			//Process the Blocks here
			if (block == Bonds) {
				if ((*k).size() == 4) {
					addBond((*k)[0],(*k)[1],MslTools::toDouble((*k)[2]),MslTools::toDouble((*k)[3]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << endl;
				} else {
					cerr << "Wrong number of params in the Bonds block. Should be 4 but it is " << (*k).size() << endl;
				}
				continue;
			} else if (block == Angles) {
				if ((*k).size() == 5) {
					addAngle((*k)[0],(*k)[1],(*k)[2],MslTools::toDouble((*k)[3]),MslTools::toDouble((*k)[4]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << " "<< (*k)[4] << endl;
				} else if ((*k).size() == 7) {
					addAngle((*k)[0],(*k)[1],(*k)[2],MslTools::toDouble((*k)[3]),MslTools::toDouble((*k)[4]),MslTools::toDouble((*k)[5]),MslTools::toDouble((*k)[6]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << " "<< (*k)[4] << " " << (*k)[5] << " "<< (*k)[6] << endl;
				} else {
					cerr << "Wrong number of params in the Angles block. Should be 5 or 7 but it is " << (*k).size() << endl;
				}
				continue;
			} else if (block == Dihedrals) {
				if ((*k).size() == 7) {
					addDihedral((*k)[0],(*k)[1],(*k)[2],(*k)[3],MslTools::toDouble((*k)[4]),MslTools::toDouble((*k)[5]),MslTools::toDouble((*k)[6]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << " "<< (*k)[4] << " " << (*k)[5] << " "<< (*k)[6] << endl;
				} else {
					cerr << "Wrong number of params in the Dihedrals block. Should be 7 but it is " << (*k).size() << endl;
				}
				continue;
			} else if (block == Improper) {
				if ((*k).size() == 7) {
					addImproper((*k)[0],(*k)[1],(*k)[2],(*k)[3],MslTools::toDouble((*k)[4]),MslTools::toDouble((*k)[6]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << " "<< (*k)[4] << " " << (*k)[5] << " "<< (*k)[6] << endl;
				} else {
					cerr << "Wrong number of params in the Improper block. Should be 7 but it is " << (*k).size() << endl;
				}
				continue;
			} else if (block == NonBonded) {
				if ((*k).size() == 4) {
					addVdw((*k)[0],MslTools::toDouble((*k)[2]),MslTools::toDouble((*k)[3]),MslTools::toDouble((*k)[2]),MslTools::toDouble((*k)[3]));
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << endl;

				} else if((*k).size() == 7) {
				//	cout << (*k)[0] << " " << (*k)[1] << " "<< (*k)[2] << " " << (*k)[3] << " "<< (*k)[4] << " " << (*k)[5] << " "<< (*k)[6] << endl;
					addVdw((*k)[0],MslTools::toDouble((*k)[2]),MslTools::toDouble((*k)[3]),MslTools::toDouble((*k)[5]),MslTools::toDouble((*k)[6]));
				} else {
					cerr << "Wrong number of params in the NonBonded block. Should be 4 or 7 but it is " << (*k).size() << (*k)[0] << endl;
				}	
				continue;
			} else if (block == HBond) {
				continue;
			} else {
				cerr << "what block is it???" << endl;
				continue;
			}
		}
			

	} catch(...){
		cerr << "ERROR 8123 in void CharmmParameterReader::read()\n";
		exit(8123);
	}

	return true;
}

vector<double> CharmmParameterReader::vdwParam(string _type) const{
	
	map<string,vector<double> >::const_iterator found;

	if ((found= vdwParamMap.find(_type)) != vdwParamMap.end()) {
		return(found->second);
	} else {
		vector<double> out(4,0.0);
		cerr << "vdwParams not found for type type " << _type << endl;
		return(out);
	}
}

vector<double> CharmmParameterReader::bondParam(string _type1, string _type2) const{

	map<string, map<string, vector<double> > >::const_iterator found1;
	map<string, vector<double> > ::const_iterator found2;

	if ((found1 = bondParamMap.find(_type1)) != bondParamMap.end() && (found2 = (found1->second).find(_type2)) != (found1->second).end()) {
		return(found2->second);
	} else {
		vector<double> out(2,0.0);
		cerr << "bondParams not found for types type (" << _type1 << "," << _type2 << ")" << endl;
		return(out);
	}

}


vector<double> CharmmParameterReader::angleParam(string _type1, string _type2, string _type3) const{
	vector<double> out(4, 0.0);
	//vector<double> out(2, 0.0);
	map<string , map<string, map<string, vector<double> >  >  >::const_iterator found1; 
	map<string, map<string, vector<double> > >::const_iterator found2;
	map<string, vector<double> >::const_iterator found3;

	if ((found1 = angleParamMap.find(_type1)) != angleParamMap.end() && (found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end()) {
		
		//Assuming the size of the vector is going to be 4. Dangerous??
		out[0] = found3->second[0];
		out[1] = found3->second[1];
	//	out.push_back(found3->second[0]);
	//	out.push_back(found3->second[1]);
	}

	return(out);

}


vector<double> CharmmParameterReader::ureyBradleyParam(string _type1, string _type2, string _type3) const{
	//vector<double> out(4, 0.0);
	vector<double> out(2, 0.0);
	map<string , map<string, map<string, vector<double> >  >  >::const_iterator found1; 
	map<string, map<string, vector<double> > >::const_iterator found2;
	map<string, vector<double> >::const_iterator found3;

	if ((found1 = angleParamMap.find(_type1)) != angleParamMap.end() && (found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end()) {
		
		//Assuming the size of the vector is going to be 4. Dangerous??
		out[0] = found3->second[2];
		out[1] = found3->second[3];
	//	out.push_back(found3->second[2]);
	//	out.push_back(found3->second[3]);

//	} else {
//		out.push_back(0.0);
//		out.push_back(0.0);
	}

	return(out);
}

vector<double> CharmmParameterReader::angleAndUreyBradleyParam(string _type1, string _type2, string _type3) const{
	vector<double> out(4, 0.0);
	map<string , map<string, map<string, vector<double> >  >  >::const_iterator found1; 
	map<string, map<string, vector<double> > >::const_iterator found2;
	map<string, vector<double> >::const_iterator found3;

	if ((found1 = angleParamMap.find(_type1)) != angleParamMap.end() && (found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end()) {
		
		//Assuming the size of the vector is going to be 4. Dangerous??
		out = found3->second;
		//out.push_back(found3->second[0]);
		//out.push_back(found3->second[1]);
		//out.push_back(found3->second[2]);
		//out.push_back(found3->second[3]);
	}

	return(out);
}


vector <vector<double> > CharmmParameterReader::dihedralParam(string _type1, string _type2, string _type3, string _type4) const {

	


	map< string, map<string , map<string, map<string, vector<vector<double> > >  >  >  >::const_iterator found1;
	map<string , map<string, map<string, vector<vector<double> > >  >  >::const_iterator found2; 
	map<string, map<string, vector<vector<double> > > >::const_iterator found3;
	map<string, vector<vector<double>  >  >::const_iterator found4;

	// if there are multiple entries charmm uses all of them and sum up the energies
	// but if there are full entries it ignore "wild-carded" (X) entries
		

	if ((found1 = dihedralParamMap.find(_type1)) != dihedralParamMap.end() && ( found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end() && (found4 = (found3->second).find(_type4)) != (found3->second).end()) {

		return(found4->second);

	} else {
	
		//wildcard  3 combinations 1, 2,3, x or x,2,3,4|x

		vector <vector<double> > out;
		vector<string> tempType1;
		tempType1.push_back(_type1);
		tempType1.push_back("X");
		tempType1.push_back("X");
		
		vector<string> tempType4;
		tempType4.push_back("X");
		tempType4.push_back(_type4);
		tempType4.push_back("X");
		
		for ( int j= 0; j < 3 ; j++) {   // Loop over each possible wildcard combination
			string type1 = tempType1[j];
			string type4 = tempType4[j]; 

		//	cout << "Checking Dihedral Map for (" << type1 << "," << _type2 << "," << _type3 << "," << type4 << ")" << endl;
			if ((found1 = dihedralParamMap.find(type1)) != dihedralParamMap.end() && ( found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end() && (found4 = (found3->second).find(type4)) != (found3->second).end()) {

				return(found4->second);
			}

		}
			vector<double> temp(3,0.0);
			out.push_back(temp);
			
			cerr << "Dihedral Params Not found for type types (" << _type1 << "," << _type2 << "," << _type3 << "," << _type4 << ")" << endl;
			return(out);
	}
}


vector<double> CharmmParameterReader::improperParam(string _type1, string _type2, string _type3, string _type4) const{

	map< string, map<string , map<string, map<string, vector<double> >  >  >  >::const_iterator found1;
	map<string , map<string, map<string, vector<double> >  >  >::const_iterator found2; 
	map<string, map<string, vector<double> > >::const_iterator found3;
	map<string, vector<double> >::const_iterator found4;

	// if there are multiple entries charmm uses all of them and sum up the energies
	// but if there are full entries it ignore "wild-carded" (X) entries
		

	if ((found1 = improperParamMap.find(_type1)) != improperParamMap.end() && ( found2 = (found1->second).find(_type2)) != (found1->second).end() && (found3 = (found2->second).find(_type3)) != (found2->second).end() && (found4 = (found3->second).find(_type4)) != (found3->second).end()) {

		return(found4->second);

	} else {
	
		//wildcard  3 combinations 1, 2|x ,3 | x, 4 

		vector<double>  out(2,0.0);
		vector<string> tempType2;
		tempType2.push_back(_type2);
		tempType2.push_back("X");
		tempType2.push_back("X");
		
		vector<string> tempType3;
		tempType3.push_back("X");
		tempType3.push_back(_type3);
		tempType3.push_back("X");
		
		for ( int j= 0; j < 3 ; j++) {   // Loop over each possible wildcard combination
			string type2 = tempType2[j];
			string type3 = tempType3[j]; 

		//	cout << "Checking Improper Map for (" << _type1 << "," << type2 << "," << type3 << "," << _type4 << ")" << endl;
			if ((found1 = improperParamMap.find(_type1)) != improperParamMap.end() && ( found2 = (found1->second).find(type2)) != (found1->second).end() && (found3 = (found2->second).find(type3)) != (found2->second).end() && (found4 = (found3->second).find(_type4)) != (found3->second).end()) {

				return(found4->second);
			}

		}
			
			cerr << "Improper Params Not found for type types (" << _type1 << "," << _type2 << "," << _type3 << "," << _type4 << ")" << endl;
			return(out);
	}
}
		
		
	/*	
		vector<double> EEF1Param(string _type); // default solvent WATER
		vector<double> EEF1Param(string _type, string _solvent);
		vector<double> IMM1Param(string _type);
		vector<double> IMM1Param(string __type, string _solvent1, string _solvent2);
	*/


void CharmmParameterReader::createVdwParamPairs(){


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
vector<double> CharmmParameterReader::vdwParamPair(string _type1, string _type2) const {

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
