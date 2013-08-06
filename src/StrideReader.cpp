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
#include "StrideReader.h"

using namespace MSL;
using namespace std;

StrideReader::StrideReader() : Reader(){
}

StrideReader::StrideReader(const string &_filename) : Reader(_filename){
	open(_filename);
}

StrideReader::StrideReader(const StrideReader &_alnreader){
}

StrideReader::~StrideReader(){}


bool StrideReader::read(){

	if (!is_open()) {
		cout << "Not open!" << endl;
		return false;
	}

	try { 
		string dataExp   = "^ASG\\s+\\S+\\s+(\\S+)\\s+(\\S+)\\s+\\S+\\s+(\\S+).*$";
		while (!endOfFileTest()){
			string line = Reader::getLine();

			// Skip blank lines.
			if (line.size() != 0) {
				// Storage for regular expression matches for this line
				vector<string> reTokens;
				
				if (MslTools::regex(line,dataExp,reTokens)){
					// Parse data lines
					secondaryStructures[reTokens[0] + "," + reTokens[1]] = reTokens[2];
				}
			}
		}
	} catch(...){
	  cerr << "ERROR 9090 in StrideReader::read()\n";
	  return false;
	}

	return true;
}
