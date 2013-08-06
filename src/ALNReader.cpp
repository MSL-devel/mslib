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
#include "ALNReader.h"

using namespace MSL;
using namespace std;



ALNReader::ALNReader() : Reader(){
  version = "";
  remark  = "";
  
}

ALNReader::ALNReader(const string &_filename) : Reader(_filename){
}

ALNReader::ALNReader(const ALNReader &_alnreader){
}

ALNReader::~ALNReader(){}


bool ALNReader::read(){

	if (!is_open()) {
		return false;
	}

	try { 
		string headerExp = " \\(*([0-9\\.])\\)*\\s+([\\S\\s]+)$";
		string dataExp   = "^(\\S+)\\s+(\\S+)\\s([0-9]+)$";
		string equivExp  = "^([\\s\\S]+)$";
		int seqIndex = 0;
		int lineIndex = 0;
		while (!endOfFileTest()){
			string line = Reader::getLine();

			// Skip blank lines.
			if (line.size() != 0) {
				// Storage for regular expression matches for this line
				vector<string> reTokens;
				
				// Need to parse header line (CLUSTAL W 2.1 ...)
				if ((lineIndex == 0) && MslTools::regex(line,headerExp,reTokens)){
					version = reTokens[0];
					remark  = reTokens[1];
					lineIndex++;
				} else if (MslTools::regex(line,dataExp,reTokens)){
					// Parse data lines
					map<string,string>::iterator it;
					it = sequences.find(reTokens[0]);
					if (it == sequences.end()){
					  sequences[reTokens[0]] = reTokens[1];
					} else {
					  sequences[reTokens[0]] += reTokens[1];
					}

					// Find index of first character not key+spaces
					seqIndex = line.find(reTokens[1]);
				} else if (MslTools::regex(line.substr(seqIndex),equivExp,reTokens)){
					// Parse Residue equivalency lines
					map<string,string>::iterator it;
					it = sequences.find("equivRes");
					if (it == sequences.end()){
					  sequences["equivRes"] = reTokens[0];
					} else {
					  sequences["equivRes"] += reTokens[0];
					}
				}
			}
		}
	} catch(...){
	  cerr << "ERROR 9090 in ALNReader::read()\n";
	  return false;
	}

	return true;
}
