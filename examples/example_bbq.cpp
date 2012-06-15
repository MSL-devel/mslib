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

#include "System.h"
#include "BBQTable.h"

using namespace MSL;

int main() {
    System sys;
    // Read in a pdb file that only includes C-alpha atoms.
    sys.readPdb("../exampleFiles/example0009_caOnly.pdb");
    BBQTable bbq("../tables/PiscesBBQTable.txt");

    // Now fill in the missing bacbone atoms for each chain
    for(int chainNum = 0; chainNum < sys.chainSize(); ++ chainNum) {
        bbq.fillInMissingBBAtoms(sys.getChain(chainNum));
    }

    // Now output a pdb with all of the backbone atoms.
    // Note: Due to the way the BBQ algorithm works, no backbone
    // atoms will be generated for the first and last resiude in a chain.
    sys.writePdb("output.pdb");
}
