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

#include "PhiPsiWriter.h"

using namespace MSL;
using namespace std;




bool PhiPsiWriter::write(PhiPsiStatistics &_stat) {


	map<string,int> counts = _stat.getPhiPsiCounts();
	map<string,int>::iterator it;

	map<string,bool> aaFound;
	map<string,bool>::iterator it2;
	for (it = counts.begin(); it != counts.end();it++){
		vector<string> tokens = MslTools::tokenize(it->first, ":");
		if (tokens.size() != 3) continue;
		if (tokens[0] == "ALL") continue;
		cout << "\tWriting "<<tokens[0]<<" "<<tokens[1]<<" "<<tokens[2]<<endl;
		it2 = aaFound.find(tokens[0]);
		if (it2 == aaFound.end()){
			aaFound[it2->first] = true;
			string line = MslTools::stringf("RESN %3s", tokens[0].c_str());
			writeln(line);
		}
		string line = MslTools::stringf("STAT %s %s %10d\n", tokens[1].c_str(),tokens[2].c_str(),it->second);
		cout << "LINE: "<<line<<endl;
		writeln(line);
	}

	return true;
}


