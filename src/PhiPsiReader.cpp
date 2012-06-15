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


#include "PhiPsiReader.h"

using namespace MSL;
using namespace std;


bool PhiPsiReader::read() {

	if (!is_open()) {
		return false;
	}

	
	try { 
		string residue = "";
		while (!endOfFileTest()){

			string line = Reader::getLine();
			if (line[0] == '#') continue;
			if (line.length() <= 1) continue;

			vector<string> toks = MslTools::tokenize(line," ",false);
			if (toks[0] == "RESN"){
				if (toks.size() != 2){
					cerr << "EXIT 1343 PhiPsiReader::read() , RESN line not enough tokens: "<<line<<endl;
					exit(1343);
				}
				residue = toks[1];
			}

			if (toks[0] == "STAT"){
				if (toks.size() != 4){
					cerr << "EXIT 1344 PhiPsiReader::read() , STAT line not enough tokens("<<toks.size()<<"): "<<line<<endl;
					exit(1344);
				}
				phiPsiStat.addStatisitics(residue,toks[1],toks[2],MslTools::toInt(toks[3]));
			}
			
		}
	} catch(...){
		cerr << "ERROR 5723 in PhiPsiReader::read()\n";
		exit(5723);
	}


	phiPsiStat.computeTotalCounts();

    return true;
}

