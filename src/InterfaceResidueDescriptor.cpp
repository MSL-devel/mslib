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

#include "InterfaceResidueDescriptor.h"
#include <cctype>
#include <algorithm>

using namespace MSL;
using namespace std;


#define std_tolower ((int(*)(int))std::tolower)

/**
 * This method will look for a given string in all of the descriptions held by this residue. 
 * It will iterate over all keys in the description table and see if it can find the given 
 * string in any of the descriptions. Currently this method is case sensitive. It will also 
 * check to see if the given string is found in the given PDB ID. This allows users to easily 
 * filter for information from a given PDB ID.
 *
 * @param     _query 	The string to look for in the description table and PDB ID.
 * @return Returns true if the _query string was found in any description and false if it was not. 
 */
bool InterfaceResidueDescriptor::isStringInResidueDescription(string &_query) {
    bool stringInResidue = false;
    string lowerCasePdb = pdbId;

    // Convert pdb to lower case.
    std::transform( pdbId.begin(), pdbId.end(), lowerCasePdb.begin(), std_tolower );
    // Convert _query to lower case.
    std::transform( _query.begin(), _query.end(), _query.begin(), std_tolower );

    if( lowerCasePdb.find(_query) != string::npos )
        stringInResidue = true;

    for( map<string, string>::iterator it=descriptionTable.begin(); it != descriptionTable.end(); ++it) {
        string currDescription = it->second;
        std::transform( currDescription.begin(), currDescription.end(), currDescription.begin(), std_tolower );
        if( currDescription.find(_query) != string::npos )
            stringInResidue = true;
    }

    return stringInResidue;
}

/**
 * This method will create a clone of a given InterfaceResidueDescriptor. 
 * All information will be duplicated, so it is safe to delete the original copy.
 *
 * @param _ird 	The InterfaceResidueDescriptor to be cloned. 
 */
void InterfaceResidueDescriptor::copy(InterfaceResidueDescriptor &_ird){
    pdbId                     = _ird.getPdbId();
    chainId                   = _ird.getChainId();
    resolution                = _ird.getResolution();
    resolutionIsValid         = _ird.isResolutionValid();
    residueName               = _ird.getResidueName();
    residueNumber             = _ird.getResidueNumber();
    residueIcode              = _ird.getResidueIcode();
    secondaryStructure        = _ird.getSecondaryStructure();
    singleChainDeltaSolventAccessibility = _ird.getSingleChainDeltaSolventAccessibility();
    // Copy over the otherChain dASA
    for(map<string, double>::iterator currIt = _ird.otherChainDeltaSolventAccessibility.begin();
        currIt != _ird.otherChainDeltaSolventAccessibility.end(); ++ currIt) {
        otherChainDeltaSolventAccessibility[currIt->first] = currIt->second;
    }

    // Copy over the aaProbs
    for(map<string, double>::iterator currIt = _ird.aaProbs.begin();
        currIt != _ird.aaProbs.end(); ++ currIt) {
        aaProbs[currIt->first] = currIt->second;
    }

    consScore = _ird.getConsScore();
    consScoreConfInterval = _ird.getConsScoreConfInterval();

    chainType                 = _ird.getChainType();
    numberBackboneContacts    = _ird.getNumberBackboneContacts();
    numberSidechainContacts   = _ird.getNumberSidechainContacts();
    numberMixedContacts       = _ird.getNumberMixedContacts();
    archiveType               = _ird.getArchiveType();
    descriptionTable.clear();
    map<string,string> desTable = _ird.getDescriptionTable();
    map<string,string>::iterator it;
    for (it = desTable.begin(); it != desTable.end();it++){
        descriptionTable[it->first] = it->second;
    }
}

