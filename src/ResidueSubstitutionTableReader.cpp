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

#include "ResidueSubstitutionTableReader.h"

using namespace std;

bool ResidueSubstitutionTableReader::read(){
		
	try { 
		int lineCount = 0;
		vector<string> labels;
		while (!endOfFileTest()){

			string line = Reader::getLine();
			if (line[0] == '#') continue;
			if (line.length() <= 1) continue;
			lineCount++;

			vector<string> toks = MslTools::tokenize(line," ",false);
                        SubstitutionTable.addSubstitutionEntry(toks);
		}
	} catch(...){
		cerr << "ERROR 6723 in ResidueSubstitutionTableReader::read()\n";
		exit(6723);
	}

    return true;
}
