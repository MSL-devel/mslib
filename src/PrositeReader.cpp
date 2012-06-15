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


#include "PrositeReader.h"

using namespace MSL;
using namespace std;





bool PrositeReader::read() {
	if (!is_open()) {
		cerr << "WARNING: 14567 PrositeReader::read () File not open" << endl;
		return false;
	}

	try { 

	       string id = "";
	       string de = "";
	       string pa = "";
	       while (!endOfFileTest()){
			string line = Reader::getLine();
			// Skip comments
			if (line.substr(0,2) == "CC") continue;

			// New section...
			if (line.substr(0,2) == "//") {

			  // Add last entry
			  if (id != "" && de != "" && pa != ""){
			    prosite_data[id][de] = pa;
			    id = "";
			    de = "";
			    pa = "";
			  }
			  continue;

			}

			// do line parsing
			if (line.substr(0,2) == "ID"){
			  vector<string> tokens = MslTools::tokenize(line, ";");
			  if (tokens.size() < 2) {
			    cerr << "ERROR PrositeReader::read() line: "<<line<<endl;
			    exit(13423);
			  }
			  if (tokens[1] != " PATTERN.") continue;
			  
			  
			  id = tokens[0].substr(2,tokens[0].length()-1);
			  id = MslTools::trim(id);
			  
			}
			if (line.substr(0,2) == "DE"){
			  de = line.substr(2,line.length()-1);
			  de = MslTools::trim(de);
			}
			if (line.substr(0,2) == "PA"){

			  string pa_line = line.substr(2,line.length()-1);
			  pa_line = MslTools::trim(pa_line);

			  // Convert prosite expression to regular expression
			  MslTools::replace(pa_line," ","",true);
			  MslTools::replace(pa_line,".","",true);
			  MslTools::replace(pa_line,"-","",true);
			  MslTools::replace(pa_line,"x",".",true);
			  MslTools::replace(pa_line,"{","[^",true);
			  MslTools::replace(pa_line,"}","]",true);
			  MslTools::replace(pa_line,"(","{",true);
			  MslTools::replace(pa_line,")","}",true);
			  MslTools::replace(pa_line,">","$",true);

			  pa += pa_line;
			}

		
		}
		return true;
		

	} catch(...){
		cerr << "ERROR 5623 in PrositeReader::read()\n";
		return false;
	}

	return true;
}
