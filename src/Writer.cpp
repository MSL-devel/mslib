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

#include "Writer.h"

// STL Includes
#include <iostream>

using namespace MSL;
using namespace std;


/*
  write(string)
  purpose to write string and check for errors.
  this is important when writing over a network.
 */
bool Writer::write(string &_formatedString){

	switch (fileHandler) {
		case cstyle:	
			/*

			char *buf = &_formatedString[0];
			int len   = _formatedString.size();
			while (true){

				size_t numCharWritten = fwrite(buf,1,len,fileStream);

				// If numCharWritten = length of string, then complete string written
				if (numCharWritten == len) return true;

				if (ferror(fileStream)){

					// Check for recoverable errors
					if (errno == EINTR || errno == EAGAIN) {
						cerr << "WARNING: Write::write(string) issue EINTR or EAGAIN error symbol found after write on string '"<<_formatedString<<"'\n";
						continue;
					}

			
					// Un-recoverable error
					cerr << "ERROR: Write::write(string) failed, number of characters written: "<<numCharWritten<< " from string '"<<_formatedString<<"'\n";
					return false;

				} else if ( numCharWritten < len ) {

					cerr << "WARNING: Write::write(string) issue, number of characters written: "<<numCharWritten<<" from string '"<<_formatedString<<"'\n";

					// If partial write try to re-write unwritten portion of string
					buf += numCharWritten;
					len -= numCharWritten;
					continue;
				} else {

					// Catch-all ; not sure why this would happen
					cerr << "ERROR: Write::write(string) failed, not sure why number of characters written is less than length of string?: "<<numCharWritten<<" '"<<_formatedString<<"'\n";
					return false;
				}
		
			}
			// Just in case... Most likely only get here.
			cerr << "ERROR: Write::write(string) failed, tried to write string; after infinte loop. '"<<_formatedString<<"'"<<endl;

			*/

			break;
		case cppstyle:
			fileStream << _formatedString;
			return true;
			break;
		case stringstyle:
			(*stringStreamPtr) << _formatedString;
			return true;
			break;
	}

	return false;
}

/*
  write(string)
  purpose to write string and check for errors.
  this is important when writing over a network.
 */
bool Writer::writeln(string &_formatedString){
	
	_formatedString += "\n";
	
	return write(_formatedString);

}
