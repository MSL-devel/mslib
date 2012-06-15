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


#include "BBQTableReader.h"

using namespace std;

using namespace MSL;


bool BBQTableReader::read(BBQTable &_bbqTable) {
  	if (!is_open()) {
		return false;
	}

  _bbqTable.deleteTableEntries();
    
    try {
        int lineCount = 0;
        vector<string> labels;

        // Loop over all of the lines in the file.
        while (!endOfFileTest()) {
            string line = Reader::getLine();
            // Skip this line if it's a comment.
            if (line[0] == '#') continue;
            // If there's less than 1 character on the line,
            // then skip it.  Nothing to see here.
            if (line.length() <= 1) continue;

            // Split the line by spaces.
            vector<string> toks = MslTools::tokenize(line, " ", false);

            // Special instructions for the first line.
            // The first line (not counting comments) should
            // give dimensions.
            if (lineCount == 0) {
                // The first three numbers are floats which
                // give the bin size (in Angstroms).
                Real binSizes[3];
                for (int i = 0; i < 3; ++i)
                    binSizes[i] = MslTools::toReal(toks[i]);
                _bbqTable.setBinSizes(binSizes[0], binSizes[1], binSizes[2]);
            } else {    // This is not the first line.
                try {
                    CartesianPoint key;
                    AtomPointerVector *av = new AtomPointerVector();
                    // The first three words on the line are the table indices.
                    key.setX( MslTools::toReal(toks[0]) );
                    key.setY( MslTools::toReal(toks[1]) );
                    key.setZ( MslTools::toReal(toks[2]) );

                    // Following the indices is a list of atom names, and x, y, and z coords.
                    for (int i = 3; i < toks.size(); i += 4) {
                        Atom *a = new Atom(toks[i], MslTools::toReal(toks[i + 1]), MslTools::toReal(toks[i + 2]), MslTools::toReal(toks[i + 3]));
                        av->push_back(a);
                    }
                    
                    _bbqTable[key] = av;
                } catch (...) {
                    cerr << "ERROR 6723 in BBQTableReader::read()\n";
                    exit(6723);
                }
            }
            lineCount++;
        }
    } catch (...) {
        cerr << "ERROR 6723 in BBQTableReader::read()\n";
        exit(6723);
    }

    return true;
}
