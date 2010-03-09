/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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


#include "Reader.h"

using namespace MSL;
using namespace std;





string Reader::getLine() {
	string line;
	line = "END"; // Is this the best default value?

	switch (fileHandler) {
		case cstyle:
			break;
		case cppstyle:
			if (!fileStream.fail() && (!fileStream.eof())){
				getline(fileStream, line);
			}
			break;
		case stringstyle:
			if (!stringStreamPtr->fail()) {
				getline(*stringStreamPtr,line);
			}
			break;
	}
	return line;
}


void Reader::copy(const Reader &_anotherReader){
	
}

bool Reader::read(string &_inputString){
	fileName = "string";
	stringStreamPtr = new stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	bool out = read();
	delete stringStreamPtr;
	return out;
}
