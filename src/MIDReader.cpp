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


#include "MIDReader.h"
#include <iostream>

using namespace MSL;
using namespace std;


/**
 * This method will read in the MoleculeInterfaceDatabase from 
 * the file, sstream, or string given by the constructor.
 */
bool MIDReader::read(string _archiveType){
 	if (!is_open()) {
		return false;
	}

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
            MslTools::splitIntAndString(fields[4], residueNumber, residueIcode);
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
                    string chain = MslTools::trim(keyValuePair[0], " ");
                    string dASAString = MslTools::trim(keyValuePair[1], " ");
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
                    string aaName = MslTools::trim(keyValuePair[0], " ");
                    string aaProbString = MslTools::trim(keyValuePair[1], " ");
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
                    string key   = MslTools::trim(keyValuePairs[j],"'");
                    string value = MslTools::trim(keyValuePairs[j+1],"'");
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
