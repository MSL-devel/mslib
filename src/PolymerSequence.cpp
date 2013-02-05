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


#include "PolymerSequence.h"
#include <map>

using namespace MSL;
using namespace std;


PolymerSequence::PolymerSequence() {
	setup("");
}

PolymerSequence::PolymerSequence(string _sequence) {
	setup(_sequence);
}

PolymerSequence::PolymerSequence(System &_sys) {
	setup("");
	setSequence(_sys);
	/*
	stringstream seq;
	for (uint c = 0; c< _sys.size();c++){

		Chain ch = _sys.getChain(c);


		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.size();p++){
			Position pos = ch.getPosition(p);

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

				seq << residueName << " ";

			}

			if (pos.size() > 1){
                                seq << "]";
                        }
		}
		seq << "\n";
	}
	//cout << "SEQ FROM SYS "<<seq.str()<<endl;
	setup(seq.str());
	*/
}

PolymerSequence::PolymerSequence(const AtomPointerVector &_atoms) {
	setup("");
	setSequence(_atoms);
}

PolymerSequence::PolymerSequence(System &_sys, vector<pair<string,string> > &_addTerminalResidues){

	stringstream seq;
	for (uint c = 0; c< _sys.chainSize();c++){

		Chain ch = _sys.getChain(c);


		seq << ch.getChainId()<<": ";

		if (_addTerminalResidues.size() > c && _addTerminalResidues[c].first != ""){ 
			cout << "ADDINGN "<<_addTerminalResidues[c].first<<"."<<endl;
			seq << _addTerminalResidues[c].first<<" ";
		}


		for (uint p = 0 ; p < ch.positionSize();p++){
			Position pos = ch.getPosition(p);

			if (pos.identitySize() > 1){
				seq << " [";
			}

			for (uint i = 0; i < pos.identitySize();i++){
				seq << " "<<pos.getIdentity(i).getResidueName();
			}

			if (pos.identitySize() > 1){
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

/*
PolymerSequence::PolymerSequence(System &_sys, map<string,map<int,int> > &_variablePositionMap, vector<vector<string> > &_identitesAtVariablePositions){
	
	stringstream seq;

	map<string,map<int,int> >::iterator it1;
	map<int,int>::iterator it2;
	for (uint c = 0; c< _sys.size();c++){
		Chain &ch = _sys.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.size();p++){
			Position &pos = ch.getPosition(p);

			cout << "Position: "<<pos.getResidueNumber()<<endl;
			
			seq << "{"<<pos.getResidueNumber();
			if (pos.getResidueIcode() != ""){
				seq << pos.getResidueIcode();
			}
			seq << "}";

			int index = -1;
			it1 = _variablePositionMap.find(pos.getChainId());
			if (it1 != _variablePositionMap.end()){
				it2 = it1->second.find(pos.getResidueNumber());
				if (it2 != it1->second.end()){

					index = it2->second;
				}
			}
			if (index != -1){


				if (_identitesAtVariablePositions[index].size() > 1){
					seq << "[";
				}

				for (uint i = 0; i < _identitesAtVariablePositions[index].size();i++){

					seq << _identitesAtVariablePositions[index][i]<<" ";
				}

				if (_identitesAtVariablePositions[index].size() > 1){
					seq << "] ";				
				}
			} else {

				string residueName = pos.getCurrentIdentity().getResidueName();
				// Replace HIS with HSD, by default
				if (residueName == "HIS") {
					cerr << "WARNING 78234: Residue HIS does not exist in topology, replacing with HSD by default!" << endl;
					residueName = "HSD";
				}


				seq << residueName<<" ";
			}

		}
		seq << "\n";
		
	}

	//cout << "Full STRING: "<<seq.str()<<endl;
	setup(seq.str());
}

*/

PolymerSequence::PolymerSequence(System &_sys, map<string,int> &_variablePositionMap, vector<vector<string> > &_identitesAtVariablePositions){
        setup("");

	stringstream seq;
	for (uint c = 0; c< _sys.chainSize();c++){
		Chain & ch = _sys.getChain(c);
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.positionSize();p++){
			Position *pos = &ch.getPosition(p);

			// If this is marked as a slave, we should use the master position
			if (pos->getLinkedPositionType() == Position::SLAVE){
			  vector<Position *> linked = pos->getLinkedPositions();
			  pos = linked[0];
			}

			seq << "{"<<pos->getResidueNumber();
			if (pos->getResidueIcode() != ""){
				seq << pos->getResidueIcode();
			}
			seq << "}";

			map<string,int>::iterator it = _variablePositionMap.find(pos->getPositionId());
			if (it != _variablePositionMap.end()){
			    seq << "[ ";
			    for (uint m = 0; m < _identitesAtVariablePositions[it->second].size();m++){
			      seq << _identitesAtVariablePositions[it->second][m];
			      if (m < _identitesAtVariablePositions[it->second].size()-1){
				seq << " ";
			      }
			      
			    }
			    seq << "] ";
			} else {
			    
				string residueName = pos->getIdentity(0).getResidueName();
				seq << residueName << " ";
			}

		}
	}

	sequenceFormatted = seq.str();
	parseString(seq.str());
}

/*
PolymerSequence::PolymerSequence(System &_sys, map<string,vector<string> > &_variablePositionMap){

	stringstream seq;

	map<string,map<int,int> >::iterator it1;
	map<int,int>::iterator it2;
	for (uint c = 0; c< _sys.size();c++){
		Chain &ch = _sys.getChain(c);
		cout << "Chain: "<<ch.getChainId()<<endl;

		
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.size();p++){
			Position &pos = ch.getPosition(p);

			cout << "Position: "<<pos.getResidueNumber()<<endl;
			
			seq << "{"<<pos.getResidueNumber();
			if (pos.getResidueIcode() != ""){
				seq << pos.getResidueIcode();
			}
			seq << "}";

			int index = -1;
			it1 = _variablePositionMap.find(pos.getChainId());
			if (it1 != _variablePositionMap.end()){
				it2 = it1->second.find(pos.getResidueNumber());
				if (it2 != it1->second.end()){

					index = it2->second;
				}
			}
			if (index != -1){


				if (_identitesAtVariablePositions[index].size() > 1){
					seq << "[";
				}

				for (uint i = 0; i < _identitesAtVariablePositions[index].size();i++){

					seq << _identitesAtVariablePositions[index][i]<<" ";
				}

				if (_identitesAtVariablePositions[index].size() > 1){
					seq << "] ";				
				}
			} else {

				string residueName = pos.getCurrentIdentity().getResidueName();
				// Replace HIS with HSD, by default
				if (residueName == "HIS") {
					cerr << "WARNING 78234: Residue HIS does not exist in topology, replacing with HSD by default!" << endl;
					residueName = "HSD";
				}


				seq << residueName<<" ";
			}

		}
		seq << "\n";
		
	}

	//cout << "Full STRING: "<<seq.str()<<endl;
	setup(seq.str());
}
*/
PolymerSequence::PolymerSequence(const PolymerSequence & _seq) {
	copy(_seq);
}

PolymerSequence::~PolymerSequence() {
}

void PolymerSequence::setup(string _sequence) {
	if (_sequence != "") {
		parseString(_sequence);
	}
	sequenceName = "SEQ";
	refSeqFlag  = false;
	refSequence = "";
	refName     = "";
	sequenceFormatted = _sequence;

	refStartResNum = 1;
	refEquilResNum = 1;
	seqEquilResNum = 1;

	pdbNamesFlag = false;

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

void PolymerSequence::setSequence(System &_sys) {
	stringstream seq;
	for (uint c = 0; c< _sys.chainSize();c++){
		Chain & ch = _sys.getChain(c);
		seq << ch.getChainId()<<": ";

		for (uint p = 0 ; p < ch.positionSize();p++){
			Position & pos = ch.getPosition(p);

			seq << "{"<<pos.getResidueNumber();
			if (pos.getResidueIcode() != ""){
				seq << pos.getResidueIcode();
			}
			seq << "}";
			if (pos.identitySize() > 1){
				seq << "[";
				for (uint i = 0; i < pos.identitySize();i++){
					string residueName = pos.getIdentity(i).getResidueName();

					// Replace HIS with HSD, by default
					if (pdbNamesFlag && residueName == "HIS") {

					  residueName = "HSD"; //default to HSD
					  if (pos.atomExists("HD1") && !pos.atomExists("HE2")){
						residueName = "HSD";
					  } else if (!pos.atomExists("HD1") && pos.atomExists("HE2")){
						residueName = "HSE";
					  } else if (pos.atomExists("HD1") && pos.atomExists("HE2")){
						residueName = "HSP";
					  }

					}

					seq << residueName;
					if (i<pos.identitySize()-1) {
						 seq << " ";
					}

				}
				seq << "] ";
			} else {
				string residueName = pos.getIdentity(0).getResidueName();

				// Replace HIS with HSD, by default
				if (pdbNamesFlag && residueName == "HIS") {

				  residueName = "HSD"; //default to HSD
				  if (pos.atomExists("HD1") && !pos.atomExists("HE2")){
				    residueName = "HSD";
				  } else if (!pos.atomExists("HD1") && pos.atomExists("HE2")){
				    residueName = "HSE";
				  } else if (pos.atomExists("HD1") && pos.atomExists("HE2")){
				    residueName = "HSP";
				  }
				}


				seq << residueName << " ";
			}

		}
	}
	parseString(seq.str());
}

void PolymerSequence::setSequence(const AtomPointerVector& _atoms) {
	System sys(_atoms);
	setSequence(sys);	
}

void PolymerSequence::parseString(string _sequence) {
	/****************************************************************************
	 *  INPUT STRING FORMAT
	 *
	 *  - Each chain is separated by a chainid plus semicolon (A:)
	 *  - If the chain-id is not given, it is defaulted in alphabetical
	 *    order (i.e. if it is the 3rd chain ":" -> "C:", even if D: was given
	 *    before)
	 *  - Lines do not matter, multi-line or single line is parsed the same
	 *  - Residue numbers are given in {34} (with insertion code {34A}
	 *  - Multiple identities at the same position are given within [ILE VAL THR]
	 *
	 *  EXAMPLES:
	 *   A: ALA VAL ILE B: {4}TYR VAL LEU C: PRO TYR {7}VAL SER D: ALA LEU {2A}ILE SER VAL
	 *   
	 *   corresponds to
	 *   chain
	 *   | resnum
	 *   | |icode
	 *   | || resname
	 *   | || |
	 *   - -- ---
	 *   A 1  ALA
	 *   A 2  VAL
	 *   A 3  ILE
	 *   B 4  TYR
	 *   B 5  VAL
	 *   B 6  LEU
	 *   C 1  PRO
	 *   C 2  TYR
	 *   C 7  VAL
	 *   C 8  SER
	 *   D 1  ALA
	 *   D 2  LEU
	 *   D 2A ILE
	 *   D 3  SER
	 *   D 4  VAL
	 * 
	 *   Multiple identites 
	 *   A: ALA [VAL LEU TRP] ILE B: {4}[TYR SER ASP] VAL LEU
	 *
	 *   A 1  ALA
	 *   A 2  VAL LEU TRP
	 *   A 3  ILE
	 *   B 4  TYR SER ASP
	 *   B 5  VAL
	 *   B 6  LEU
	 * 
	 ****************************************************************************/
	sequence.clear();
	residueNumbers.clear();
	chainIds.clear();

	vector<string> defaultChains;
	defaultChains.push_back("A");
	defaultChains.push_back("B");
	defaultChains.push_back("C");
	defaultChains.push_back("D");
	defaultChains.push_back("E");
	defaultChains.push_back("F");
	defaultChains.push_back("G");
	defaultChains.push_back("H");
	defaultChains.push_back("I");
	defaultChains.push_back("J");
	defaultChains.push_back("K");
	defaultChains.push_back("L");
	defaultChains.push_back("M");
	defaultChains.push_back("N");
	defaultChains.push_back("O");
	defaultChains.push_back("P");
	defaultChains.push_back("Q");
	defaultChains.push_back("R");
	defaultChains.push_back("S");
	defaultChains.push_back("T");
	defaultChains.push_back("U");
	defaultChains.push_back("V");
	defaultChains.push_back("W");
	defaultChains.push_back("X");
	defaultChains.push_back("Y");
	defaultChains.push_back("Z");
	defaultChains.push_back("a");
	defaultChains.push_back("b");
	defaultChains.push_back("c");
	defaultChains.push_back("d");
	defaultChains.push_back("e");
	defaultChains.push_back("f");
	defaultChains.push_back("g");
	defaultChains.push_back("h");
	defaultChains.push_back("i");
	defaultChains.push_back("j");
	defaultChains.push_back("k");
	defaultChains.push_back("l");
	defaultChains.push_back("m");
	defaultChains.push_back("n");
	defaultChains.push_back("o");
	defaultChains.push_back("p");
	defaultChains.push_back("q");
	defaultChains.push_back("r");
	defaultChains.push_back("s");
	defaultChains.push_back("t");
	defaultChains.push_back("u");
	defaultChains.push_back("v");
	defaultChains.push_back("w");
	defaultChains.push_back("x");
	defaultChains.push_back("y");
	defaultChains.push_back("z");
	defaultChains.push_back("0");
	defaultChains.push_back("1");
	defaultChains.push_back("2");
	defaultChains.push_back("3");
	defaultChains.push_back("4");
	defaultChains.push_back("5");
	defaultChains.push_back("6");
	defaultChains.push_back("7");
	defaultChains.push_back("8");
	defaultChains.push_back("9");
	// note, segfault is we pass a string with more than 62 :, error is not catched

	// split by end-of-line and re-join into a single line
	vector<string> lines = MslTools::tokenize(_sequence, "\n", true);
	lines = MslTools::joinConnectedLines(lines);
	string rejoined = MslTools::joinLines(lines, " ");
	// identify the beginning of chains (A:)
	size_t pos  = rejoined.find(":");
	vector<string> newLines;
	if (pos == std::string::npos) {
		// no ":" given, add one at the beginning are we are done, one chain
		rejoined = defaultChains[0] + ": " + rejoined;
		newLines.push_back(rejoined);
	} else {
		unsigned int chainCounter = 0;
		unsigned int prevPos = 0;
		if (pos == 0) {
			// starts with a ": but not letter, add a letter
			rejoined = defaultChains[0] + rejoined;
			pos = 1;
		} else if (pos > 1 && MslTools::isWhiteSpaces(rejoined.substr(0, pos))) {
			// starts with "    :", add a letter
			rejoined = defaultChains[0] + rejoined.substr(pos, rejoined.size()-pos);
			pos = 1;
		} else if (pos > 1 && !MslTools::isWhiteSpaces(rejoined.substr(0, pos-1))) {
			// something like "ALA ILE B: ALA VAL"
			// add a letter at the beginning"A:ALA ILE B: ALA VAL"
			rejoined = defaultChains[0] + ": " + rejoined;
			pos = 1;
		}
		while (pos != std::string::npos) {
			// look for all the ":"
			string previous = rejoined.substr(pos-1, 1);
			if (pos == 0) {
				rejoined = defaultChains[chainCounter] + rejoined;
				pos++;
			} else {
				if (!MslTools::isAlphaNumericChars(rejoined.substr(pos-1, 1))) {
					// a semicolon without a letter was given, let's default the chain ID
					rejoined = rejoined.substr(0, pos-1) + " " + defaultChains[chainCounter] + rejoined.substr(pos, rejoined.size() -pos);
					pos++;
				}
			}
			if (prevPos> 0) {
				newLines.push_back(rejoined.substr(prevPos-1, pos-prevPos));
			}
			prevPos = pos;
			pos  = rejoined.find(":", pos+1);
			chainCounter++;
		}
		newLines.push_back(rejoined.substr(prevPos-1, rejoined.size()-prevPos+1));
	}
	lines = newLines;


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
			//    A: {12A}ALA {13}[VAL LEU] ...
			createNew = true;
			chainId = col[0];
			col.erase(col.begin());
			*k = MslTools::joinLines(col);
			vector<string> split = MslTools::tokenize(chainId);
			
			if (split.size() > 1) {
				chainId = split[0];
				resnum = split[1];
				cerr << "ERROR 43234: Depricated API. Use '"<<chainId<<": {"<<split[1]<<"} XXX YYY ZZZ'"<<endl;
				exit(43234);
			} else {
				chainId = split[0];
				resnum = "1";
			}
		}

		
		// remove all the spaces after the curlies "{34} LEU" -> "{34}LEU"
		for (unsigned int i=0; i<k->size(); i++) {
			if (k->substr(i,1) == (string)"}") {
				while (i<k->size() - 1 && k->substr(i+1,1) == " ") {
					k->erase(i+1, 1);
				}
			}
		}

		// put the curlies inside the squares "{34}[ILE LEU VAL]" -> "[{34}ILE LEU VAL]"
		unsigned int start = 0;
		for (unsigned int i=0; i<k->size(); i++) {
			if (k->substr(i,1) == (string)"{") {
				start = i;
			}
			if (k->substr(i,1) == (string)"}") {
				if (i<k->size() - 1 && k->substr(i+1,1) == "[") {
					string curlyNum = k->substr(start, i-start+1);
					k->erase(start, i-start+1);
					k->insert(start+1, curlyNum);
					i++;
					while (i<k->size() - 1 && k->substr(i+1,1) == " ") {
						k->erase(i+1, 1);
					}
				}
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
					sprintf(tmp,"%d",(int)(lastResNum));
					residueNumbers.back().push_back(tmp);
				}


			}


			sequence.back().insert(sequence.back().end(), res2.begin(), res2.end());
		}
	}
	/*
	for (vector<vector<vector<string> > >::iterator chain=sequence.begin(); chain!=sequence.end(); chain++) {
		for (vector<vector<string>  >::iterator residue=chain->begin(); residue!=chain->end(); residue++) {

			// DON'T CHANGE THIS FORMATING UNTIL YOU LOOK AT the toString function. or else printing out aligned sequences will be broken.
			//char tmp[5];
			//sprintf(tmp,"%04d",(int)(resStart[chain - sequence.begin()] + residue - chain->begin()));
			//residueNumbers[chain - sequence.begin()].push_back(tmp);
		}
	}
	*/
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
			//ss << "Chain " << chainIds[chain - sequence.begin()] << endl;
			ss << chainIds[chain - sequence.begin()] << ": ";
			int prevResnum = 0;
			string prevICode = "";
			for (vector<vector<string>  >::const_iterator residue=chain->begin(); residue!=chain->end(); residue++) {
				//ss << ">" << residue-chain->begin() << "<";
				int resNumPart = 0;
				string iCodePart = "";
				MslTools::splitIntAndString(residueNumbers[chain - sequence.begin()][residue - chain->begin()], resNumPart, iCodePart);
				bool printResnum = false;
				/* 
				  Print resnum and icode for the residue if
				   - it is the first of the chain
				   - it is the last of the chain
				   - it has a non-blank insertion code
				   - the resnum is not incremented from the previous one by 1
				*/
				if (residue == chain->begin() || residue == chain->end()-1 || iCodePart != "" || resNumPart != prevResnum + 1) {
					printResnum = true;
				} else {
				}
				if (printResnum) {
					ss << "{" << residueNumbers[chain - sequence.begin()][residue - chain->begin()] << "}";
				}
				prevResnum = resNumPart;
				prevICode = iCodePart;
				for (vector<string>::const_iterator id=residue->begin(); id!=residue->end(); id++) {
					if (residue->size() > 1 && id==residue->begin()) {
						ss << "[";
					}
					ss << *id;
					if (residue->size() > 1) {
						if (id==residue->end()-1) {
							ss << "]";
						} else {
							ss << " ";
						}
					}
				}
				if (residue!=chain->end()-1) {
					ss << " ";
				}

				//if (residue == chain->end()-1) {
				//	//ss << " " << resStart[chain - sequence.begin()] + chain->end() - chain->begin() - 1;
				//	//ss << " " << *(residueNumbers[chain - sequence.begin()].end()-1) << " ";
				//	ss << (residueNumbers[chain - sequence.begin()].back())<<endl;
				//}
			}
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



string PolymerSequence::toThreeLetterCode(AtomPointerVector &_av,string _residueDefiningAtomType){

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

string PolymerSequence::toOneLetterCode(AtomPointerVector &_av,string _residueDefiningAtomType){


	stringstream ss;
	for (uint i = 0; i < _av.size();i++){


		if (_av(i).getName() != _residueDefiningAtomType) continue;
		ss<< MslTools::getOneLetterCode(_av(i).getResidueName());
	}


	return ss.str();
}
