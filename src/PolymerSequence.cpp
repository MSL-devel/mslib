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

#include "PolymerSequence.h"
#include <map>

PolymerSequence::PolymerSequence() {
	setup("");
}

PolymerSequence::PolymerSequence(string _sequence) {
	setup(_sequence);
}

PolymerSequence::PolymerSequence(System &_sys) {

	stringstream seq;
	for (uint c = 0; c< _sys.size();c++){

		Chain ch = _sys.getChain(c);



		seq << ch.getChainId()<<" " << ch.getPositionByIndex(0).getResidueNumber()<<": ";


		for (uint p = 0 ; p < ch.size();p++){
			Position pos = ch.getPositionByIndex(p);

			seq << "{"<<pos.getResidueNumber();
			if (pos.getResidueIcode() != ""){
				seq << pos.getResidueIcode();
			}
			seq << "}";
			if (pos.size() > 1){
				seq << " [";
			}

			for (uint i = 0; i < pos.size();i++){

				string residueName = pos.getIdentity(i).getResidueName();
				// Replace HIS with HSD, by default
				if (residueName == "HIS") {
					cerr << "WARNING 78234: Residue HIS does not exist in topology, replacing with HSD by default!" << endl;
					residueName = "HSD";
				}
				seq << " "<<residueName;

			}

			if (pos.size() > 1){
				seq << "]";
			}

		}
		seq << "\n";
	}
	//cout << "SEQ FROM SYS "<<seq.str()<<endl;
	setup(seq.str());
}

PolymerSequence::PolymerSequence(System &_sys, vector<pair<string,string> > &_addTerminalResidues){

	stringstream seq;
	for (uint c = 0; c< _sys.size();c++){

		Chain ch = _sys.getChain(c);


		seq << ch.getChainId()<<": ";

		if (_addTerminalResidues.size() > c && _addTerminalResidues[c].first != ""){ 
			cout << "ADDINGN "<<_addTerminalResidues[c].first<<"."<<endl;
			seq << _addTerminalResidues[c].first<<" ";
		}


		for (uint p = 0 ; p < ch.size();p++){
			Position pos = ch.getPositionByIndex(p);

			if (pos.size() > 1){
				seq << " [";
			}

			for (uint i = 0; i < pos.size();i++){
				seq << " "<<pos.getIdentity(i).getResidueName();
			}

			if (pos.size() > 1){
				seq << "]";
			}

		}

		if (_addTerminalResidues.size() > c && _addTerminalResidues[c].second != ""){  
			cout << "ADDINGC "<<_addTerminalResidues[c].second<<"."<<endl;
			seq << " "<<_addTerminalResidues[c].second;
		}
	}

	setup(seq.str());
}

PolymerSequence::PolymerSequence(System &_sys, map<string,map<int,int> > _variablePositionMap, vector<vector<string> > _identitesAtVariablePositions){

	stringstream seq;

	map<string,map<int,int> >::iterator it1;
	map<int,int>::iterator it2;
	for (uint c = 0; c< _sys.size();c++){
		Chain &ch = _sys.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<" "<<ch.getPositionByIndex(0).getResidueNumber()<<": ";

		for (uint p = 0 ; p < ch.size();p++){
			Position &pos = ch.getPositionByIndex(p);

			cout << "Position: "<<pos.getResidueNumber()<<endl;
			
			
			int index = -1;
			it1 = _variablePositionMap.find(pos.getChainId());
			if (it1 != _variablePositionMap.end()){
				it2 = it1->second.find(pos.getResidueNumber());
				if (it2 != it1->second.end()){

					index = it2->second;
				}
			}
			if (index != -1){

				seq << " [";

				for (uint i = 0; i < _identitesAtVariablePositions[index].size();i++){
					seq << " "<<_identitesAtVariablePositions[index][i];
				}

				seq << "]";				
			} else {

				seq << " "<<pos.getCurrentIdentity().getResidueName();
			}

		}
		seq << "\n";
		
	}

	setup(seq.str());
}
PolymerSequence::PolymerSequence(const PolymerSequence & _seq) {
	copy(_seq);
}

PolymerSequence::~PolymerSequence() {
}

void PolymerSequence::setup(string _sequence) {
	parseString(_sequence);
	sequenceName = "SEQ";
	refSeqFlag  = false;
	refSequence = "";
	refName     = "";

	refStartResNum = 1;
	refEquilResNum = 1;
	seqEquilResNum = 1;

}

void PolymerSequence::copy(const PolymerSequence & _seq) {
	sequence = _seq.sequence;
	residueNumbers = _seq.residueNumbers;
	chainIds = _seq.chainIds;

	refSeqFlag  = _seq.refSeqFlag;
	refSequence = _seq.refSequence;
	refName     = _seq.refName;

	refStartResNum = _seq.refStartResNum;
	refEquilResNum = _seq.refEquilResNum;
	seqEquilResNum = _seq.seqEquilResNum;

}

void PolymerSequence::setSequence(string _sequence) {
	parseString(_sequence);
}


void PolymerSequence::parseString(string _sequence) {
	sequence.clear();
	residueNumbers.clear();
	chainIds.clear();
	vector<string> lines = MslTools::tokenize(_sequence, "\n", true);
	//cout << "====================" << endl;
	//cout << _sequence;
	//cout << "====================" << endl;
	for (vector<string>::iterator k=lines.begin(); k!=lines.end(); k++) {
		//cout << "$ " << *k << endl;
	}
	//lines = MslTools::joinBackslashedLines(lines, " ");
	lines = MslTools::joinConnectedLines(lines);
	bool createNew = true;
	string chainId = "A";
	string resnum = "1";
	string iCode = "";
	vector<int> resStart;
	for (vector<string>::iterator k=lines.begin(); k!=lines.end(); k++) {
		//cout << "% " << *k << endl;
		vector<string> col = MslTools::tokenize(*k, ":");
		if (col.size() > 1) {
			// the line is in the format
			//    A: ALA LEU VAL...
			// or
			//    A 12: ALA LEU VAL...
			createNew = true;
			chainId = col[0];
			col.erase(col.begin());
			*k = MslTools::joinLines(col);
			vector<string> split = MslTools::tokenize(chainId);
			
			if (split.size() > 1) {
				chainId = split[0];
				resnum = split[1];
				cerr << "WARNING : Depricated API. Use '"<<chainId<<": {"<<split[1]<<"} XXX YYY ZZZ'"<<endl;
			} else {
				chainId = split[0];
				resnum = "1";
			}
		}

		


		/******************************************************
		 * Resolve multiple identities if they were given
		 * 
		 * The dualLevelTokenizer creates a 2D vector based on
		 * square brackets
		 *
		 *    0,0  1,0 1,1  2,0  3,0 3,1 3,2 3,3
		 *    ILE [ALA VAL] LEU [PHE SER ILE TRP]
		 ******************************************************/
		vector<vector<string> > res2 = MslTools::dualLevelTokenizer(*k);
	//	cout << "UUU Line " << k-lines.begin() << ", " << res2.size() << " tokens" << endl;
		if (res2.size() == 0) {
			createNew = true;
			//cout << "UUU emtpy turn createNew on" << endl;
			continue;
		} else {
			if (col.size() > 0 || createNew) {
			//	cout << "UUU createNew on, new chain" << endl;
				sequence.push_back(vector<vector<string> >());
				chainIds.push_back(chainId);
				residueNumbers.push_back(vector<string>());
				resStart.push_back(MslTools::toInt(resnum));
				MslTools::splitIntAndString(resnum, resStart.back(), iCode);	
			}


			// Check for {}, residue numbers + insertion codes.
			int lastResNum = resStart.back()-1;
			for (uint s = 0; s < res2.size();s++){

				//cout << "RES2: "<<res2[s][0]<<endl;
				// This will give back before,inside,after the {}, so HI{FOO}BYE = [ HI, FOO, BYE ].
				vector<string> residueSequenceInfo = MslTools::extractBraketed(res2[s][0],"{","}");

				if (residueSequenceInfo[1] != ""){
					residueNumbers.back().push_back(residueSequenceInfo[1]);
					int resn;
					string icode;
					MslTools::splitIntAndString(residueSequenceInfo[1],resn,icode);
					lastResNum = resn;


					// Remove {..} from res2[s][0]...
					res2[s][0] = residueSequenceInfo[2];

				} else {
					
					lastResNum++;
					char tmp[5];
					sprintf(tmp,"%04d",(int)(lastResNum));
					residueNumbers.back().push_back(tmp);
				}


			}


			sequence.back().insert(sequence.back().end(), res2.begin(), res2.end());
			//cout << "UUU size of chain " << sequence.back().size() << endl;
			//cout << endl;
		}
	}
	for (vector<vector<vector<string> > >::iterator chain=sequence.begin(); chain!=sequence.end(); chain++) {
		for (vector<vector<string>  >::iterator residue=chain->begin(); residue!=chain->end(); residue++) {

			// DON'T CHANGE THIS FORMATING UNTIL YOU LOOK AT the toString function. or else printing out aligned sequences will be broken.
			//char tmp[5];
			//sprintf(tmp,"%04d",(int)(resStart[chain - sequence.begin()] + residue - chain->begin()));
			//residueNumbers[chain - sequence.begin()].push_back(tmp);
		}
	}
}


string PolymerSequence::getReferenceHeader(){
		
	stringstream ss;
	stringstream numberLine;
	stringstream dashLine;
	for (uint i =0; i < refSequence.length();i++){
			
		if (i % 10 == 0 ){

			char tmp[5];
			sprintf(tmp, "%-4d",(refStartResNum+i));
			numberLine << tmp;
			dashLine   << "|";
		} else {
			// The res number is %-4d, so skip the next for characters..
			if (!(i % 10 == 1 || i % 10 == 2 || i % 10 == 3)) { 
				numberLine << " ";
			}

			if (i % 5 == 0){
				dashLine   << "+";
			} else {
				dashLine   << "-";
			}
		}


	}


	char buffer[80];
	sprintf(buffer, "                                                         ");
	ss << buffer << numberLine.str()<<endl;
	ss << buffer << dashLine.str()<<endl;
	char name[80];
	sprintf(name, "%-45s ( 0000)    ",refName.substr(0,45).c_str());
	ss << name << refSequence<<endl;

	return ss.str();

}
string PolymerSequence::toString() const {

	stringstream ss;
	if (refSeqFlag){
		/*
		  FORMAT



                                                        70        80        90       100       110       120       130     
                                                   +----|----+----|----+----|----+----|----+----|----+----|----+----|----+
REFNAME(STARTRES)                                  ABCDEFGHIGYABCDEFGHIGYABCDEFGHIGYABCDEFGHIGYABCDEFGHIGYABCDEFGHIGYABCDE
2e74-000_001-0032_0056_A-0080_0105_A.pdb(STARTRES)                       XFXXLGXITXXCFXIQXXTGXAMXX                         

		*/



		int diffPolySeq = seqEquilResNum - MslTools::toInt(residueNumbers[0][0]);
		int resNum      = refEquilResNum-diffPolySeq; 
		int offset      = resNum - refStartResNum;
		
		//		int diffInResNumbers = refStartPolySeqResidueNum - refStartRefResidueNum; 
		//		int resNum = MslTools::toInt(residueNumbers[0][0])+diffInResNumbers;
		

		char sname[80];
		sprintf(sname, "%-45s (%1s%04d)    ",sequenceName.substr(0,45).c_str(),getChainId(0).c_str(),seqEquilResNum);
		ss << sname;

		//cout << "NUMBERS: "<<residueNumbers[0][0]<<","<<seqEquilResNum<<" and "<<refStartResNum<<","<<refEquilResNum<<" ---> "<<diffPolySeq<<" "<<resNum<<" "<<offset<<endl;

		string seqBuffer;
		if (offset > -1){
			for (uint i = 0; i < offset;i++){
				seqBuffer += " ";

			}
			ss <<seqBuffer;
		}

		for (vector<vector<vector<string> > >::const_iterator chain=sequence.begin(); chain!=sequence.end(); chain++) {
			
			for (vector<vector<string>  >::const_iterator residue=chain->begin(); residue!=chain->end(); residue++) {
				if (offset >= 0){
					string res = MslTools::getOneLetterCode((*residue)[0]);
					ss << res;
				} else{
					offset++;
				}
			}
		}
		
		ss <<endl;
		


	} else{

		for (vector<vector<vector<string> > >::const_iterator chain=sequence.begin(); chain!=sequence.end(); chain++) {
			ss << "Chain " << chainIds[chain - sequence.begin()] << endl;
			for (vector<vector<string>  >::const_iterator residue=chain->begin(); residue!=chain->end(); residue++) {
				if (residue == chain->begin()) {
					//ss << resStart[chain - sequence.begin()] << " ";
					ss << *(residueNumbers[chain - sequence.begin()].begin()) << " ";
				}
				for (vector<string>::const_iterator id=residue->begin(); id!=residue->end(); id++) {
					if (residue->size() > 1 && id==residue->begin()) {
						ss << "[";
					}
					ss << *id;
					if (residue->size() > 1 && id==residue->end()-1) {
						ss << "] ";
					} else {
						ss << " ";
					}
				}
				if (residue == chain->end()-1) {
					//ss << " " << resStart[chain - sequence.begin()] + chain->end() - chain->begin() - 1;
					//ss << " " << *(residueNumbers[chain - sequence.begin()].end()-1) << " ";
					ss << (residueNumbers[chain - sequence.begin()].back())<<endl;
				}
			}
			ss << endl;
			ss << endl;
		}
	}
	return ss.str();
}








void PolymerSequence::setReferenceSequence(string _refSeq, string _refName, int _startRefResidueNumber, int _equivalentRefRes, int _equivalentPolyRes){

	refSeqFlag  = true;
	refSequence = _refSeq;
	refName     = _refName;

	refStartResNum = _startRefResidueNumber;
	refEquilResNum = _equivalentRefRes;
	seqEquilResNum = _equivalentPolyRes;

}



string PolymerSequence::toThreeLetterCode(AtomVector &_av,string _residueDefiningAtomType){

	stringstream ss;
	int j = 0;
	for (uint i = 0; i < _av.size();i++){if (_av(i).getName() != _residueDefiningAtomType) continue;

		if (j > 0){
			ss << " ";
		}
		ss<< _av(i).getResidueName();
		j++;
	}


	return ss.str();
}

string PolymerSequence::toOneLetterCode(AtomVector &_av,string _residueDefiningAtomType){


	stringstream ss;
	for (uint i = 0; i < _av.size();i++){


		if (_av(i).getName() != _residueDefiningAtomType) continue;

		ss<< MslTools::getOneLetterCode(_av(i).getResidueName());
	}


	return ss.str();
}
