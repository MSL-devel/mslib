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

#include "ResidueSubstitutionTable.h"

using namespace std;

using namespace MSL;


void ResidueSubstitutionTable::addSubstitutionEntry(std::vector<string> residueEntry) {
    string badResidueName;
    ReplacementInfo ri;
    // If there are no entries, then just return.
    if(residueEntry.size() < 1)
        return;

    // The first string is the name of the residue we want to replace.
    badResidueName = residueEntry[0];

    // If there are more strings, the 2nd is the residue name we should replace
    // the bad residue with, and the remaining items are the atoms of the bad
    // residue which we should keep.
    if(residueEntry.size() > 1) {
        AtomHash ah;
        string goodResidueName = residueEntry[1];

        // Any additional information is a list of atoms that we should keep.
        for(int i=2; i < residueEntry.size(); ++i) {
            ah.insert( pair<string, bool>(residueEntry[i], true) );
        }

        // Create our replacement info, a pair of the residue name we should
        // switch to, and a hash of the atoms we should keep.
        ri.first = goodResidueName;
        ri.second = ah;
    }

    // Insert the replacement info into our map.
    insert(pair<string, ReplacementInfo>(badResidueName, ri) );
}

bool ResidueSubstitutionTable::isResidueInSubstitutionTable(Residue &_res) {
    return( find(_res.getResidueName()) != end() );
}

Residue ResidueSubstitutionTable::replaceResidue(Residue &_res) {
    // If this residue isn't in our list of residues to be replaced,
    // just return it!
    if( !isResidueInSubstitutionTable(_res) )
        return _res;

    // Clone the current residue, and find the replacement info.
    Residue newRes;
    ReplacementInfo ri = find( _res.getResidueName() )->second;
    string newName = ri.first;
    AtomHash ah = ri.second;

    // Loop over all of the atoms in the residue.
    for(int i=0; i < _res.size(); ++i ){
        Atom &currAtom = _res.getAtom(i);
        string currAtomName = currAtom.getName();
        // If the current atom is in the list of atoms we should keep,
        // add it.
        if( ah.find(currAtomName) != ah.end() )
            newRes.addAtom( currAtom );
    }

    // Set the new name of the residue.
    newRes.setResidueNumber( _res.getResidueNumber() );
    newRes.setResidueIcode( _res.getResidueIcode() );
    newRes.setChainId( _res.getChainId() );
    newRes.setResidueName( newName );

    return newRes;
}
