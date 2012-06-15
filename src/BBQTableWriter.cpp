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

#include "BBQTableWriter.h"
#include "Real.h"
#include <sstream>

using namespace std;

using namespace MSL;


/**
 * This function will write a BBQTable to file.
 *
 * @param _bbqTable  The BBQ Table to write.
 */
bool BBQTableWriter::write(BBQTable &_bbqTable) {
    stringstream ss;
    string lineOfText;
    Real binSizes[3];
    _bbqTable.getBinSizes(binSizes[0], binSizes[1], binSizes[2]);

    // The first line gives the sizes of the bins.
    ss << binSizes[0] << " " << binSizes[1] << " " << binSizes[2] << "\n";
    lineOfText = ss.str();
    writeln(lineOfText);
    ss.str("");

    // Now loop over every entry in the BBQ Table and write it to file.
    for (BBQTable::iterator currIter = _bbqTable.begin(); currIter != _bbqTable.end(); ++currIter) {
        CartesianPoint coord = currIter->first;
        AtomPointerVector *av = currIter->second;

        // Write out the key (the r coords) for this entry.
        ss << coord.getX() << " " << coord.getY() << " " << coord.getZ();

        // Now write out the atoms and their coords for this entry.
        for (AtomPointerVector::iterator currAvIter = av->begin(); currAvIter != av->end(); ++currAvIter) {
            coord = (*currAvIter)->getCoor();
            ss << " " << (*currAvIter)->getName() << " " << coord.getX() << " " << coord.getY() << " " << coord.getZ();
        }
        // New line for the next entry.
        lineOfText = ss.str();
        // Actually write this line.
        writeln(lineOfText);
        ss.str("");
    }

    return true;
}

