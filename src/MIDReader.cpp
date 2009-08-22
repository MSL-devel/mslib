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

#include "MIDReader.h"
#include <iostream>

/**
 * This method will read in the MoleculeInterfaceDatabase from 
 * the file, sstream, or string given by the constructor.
 */
bool MIDReader::read(string _archiveType){
    try {
        int lineCount = 1;
        string lastPDB = "";
        while (!endOfFileTest()){
            string line = MslTools::trim(Reader::getLine());
            vector<string> fields = MslTools::tokenizeAndTrim(line,";",true);

            // Skip over comment lines..
            if (line[0] == '#') continue;

            // Skip over lines that don't have enough fields.
            if (fields.size() < 17) continue;

            // Check for blank residue number.. should check all fields?
            if (fields[4] == "") continue;

            if (lastPDB != fields[0]){
                cout << "Working on PDB: "<<fields[0]<<endl;
                lastPDB = fields[0];
            }

            InterfaceResidueDescriptor *ird = new InterfaceResidueDescriptor();
            ird->setPdbId(fields[0]);
            ird->setResolution( fields[1] );
            ird->setChainId(fields[2]);
            ird->setResidueName( fields[3] );
            ird->setArchiveType(_archiveType);

            int residueNumber;
            string residueIcode;
            splitIntAndString(fields[4], residueNumber, residueIcode);
            ird->setResidueNumber(residueNumber);
            ird->setResidueIcode(residueIcode);
            ird->setSecondaryStructure(fields[5]);
            if(fields[6] != "")
                ird->setSingleChainDeltaSolventAccessibility(MslTools::toDouble(fields[6],"PIDReader::read() trying to parse 4th token for solvent accessibility on line "+MslTools::intToString(lineCount)+"\n"));
            else
                ird->setSingleChainDeltaSolventAccessibility(0.0f);
            // Now get otherChainDeltaSolventAccessibility, if it exists.
            if(fields[7] != "") {
                vector<string> openBraces = MslTools::tokenize(fields[7], "{", true);
                vector<string> closeBraces = MslTools::tokenize(openBraces[1], "}", true);
                vector<string> toks = MslTools::tokenize(closeBraces[0], ",", true);
                for (uint i = 0; i < toks.size(); i++){
                    vector<string> keyValuePair = MslTools::tokenize(toks[i], "=>", true);
                    string chain = trim(keyValuePair[0], " ");
                    string dASAString = trim(keyValuePair[1], " ");
                    double dASA = MslTools::toDouble(dASAString, "PIDReader::read() trying to parse other chain dASA on line " + MslTools::intToString(lineCount)+"\n");
                    ird->setOtherChainDeltaSolventAccessibility(chain, dASA);
                }
            }

            ChainType_t type = Other;
            if (fields[8].find("Protein",0) != string::npos){
                type = Protein;
            } else if (fields[8].find("Peptide") != string::npos){
                type = Peptide;
            }

            ird->setChainType(type);                
            ird->setNumberBackboneContacts(MslTools::toInt(fields[9], "PIDReader::read() trying to parse 9th token as int for backbone contacts on line "+MslTools::intToString(lineCount)+"\n"));
            ird->setNumberSidechainContacts(MslTools::toInt(fields[10], "PIDReader::read() trying to parse 10th token as int for sidechain contacts on line "+MslTools::intToString(lineCount)+"\n"));
            ird->setNumberMixedContacts(MslTools::toInt(fields[11], "PIDReader::read() trying to parse 11th token as int for mixed contacts on line "+MslTools::intToString(lineCount)+"\n"));
            
            // Now get aaProbs, if it exists.
            if( (fields[12] != "") && (fields[12] != "{}") ) {
                vector<string> openBraces = MslTools::tokenize(fields[12], "{", true);
                vector<string> closeBraces = MslTools::tokenize(openBraces[1], "}", true);
                vector<string> toks = MslTools::tokenize(closeBraces[0], ",", true);
                for (uint i = 0; i < toks.size(); i++){
                    vector<string> keyValuePair = MslTools::tokenize(toks[i], "=>", true);
                    string aaName = trim(keyValuePair[0], " ");
                    string aaProbString = trim(keyValuePair[1], " ");
                    double aaProb = MslTools::toDouble(aaProbString, "PIDReader::read() trying to parse aa prob " + MslTools::intToString(lineCount)+"\n");
                    ird->setAAProb(aaName, aaProb);
                }
            }

            if( fields[13] != "" ) { 
                ird->setConsScore( MslTools::toDouble(fields[13]) );
            }

            if( (fields[14] != "" ) && (fields[15] != "" ) ){
                ird->setConsScoreConfInterval( MslTools::toDouble(fields[14]), MslTools::toDouble(fields[15]) );
            }

            // Now get elements that are in {}.
            vector<string> openBraces = MslTools::tokenize(fields[16], "{", true);
            vector<string> closeBraces = MslTools::tokenize(openBraces[1], "}", true);
            vector<string> toks = MslTools::tokenize(closeBraces[0], ";", true);
            for (uint i = 0; i < toks.size();i++){
                vector<string> keyValuePairs = MslTools::tokenize(toks[i], "=>", true);
                for(uint j = 0; j < keyValuePairs.size(); j+=2) {
                    string key   = trim(keyValuePairs[j],"'");
                    string value = trim(keyValuePairs[j+1],"'");
                    //cout << "Key,Value: "<<key<<","<<value<<endl;
                    ird->addToDescriptionTable(key, value);
                }
            }
            
            // This will make a copy 
            mid.addInterafaceResidueDescriptor(ird);


            // Delete the temp pointer
            delete(ird);

            // Increment line counting
            lineCount++;


            
        }
    } catch(...){
        cerr << "ERROR 5624 in PIDReader::read();\n";
        exit(5623);
    }

    return true;
}
