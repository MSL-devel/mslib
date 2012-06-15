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
#include "RandomSeqGenerator.h"
#include "Reader.h"
#include "MslTools.h"
#include <sstream>

using namespace std;

using namespace MSL;


RandomSeqGenerator::RandomSeqGenerator(const RandomSeqGenerator & _rsg) {
    seqTable = _rsg.seqTable;
    rng = _rsg.rng;
}

RandomSeqGenerator::RandomSeqGenerator(std::string &_seqTableFileName) {
    readSeqTable(_seqTableFileName);
    rng.seed(static_cast<unsigned int>(std::time(0)));
}

void RandomSeqGenerator::readSeqTable(std::string &_seqTableFileName) {
    Reader r(_seqTableFileName);
    r.open();
    double cumProb = 0.0;
    std::map<double, std::string> tempTable;
    
    string line = r.getLine();
    while(!r.endOfFileTest()) {
        // Skip over comments.
        if(line.find("#") != 0) {
            // Table should be comma separated pairs.
            vector<string> toks = MslTools::tokenize(line, ",", false);
            if(toks.size() == 2) {
                string seq = toks[0];
                double prob = MslTools::toDouble(toks[1]);
                cumProb += prob;

                // Skip entries with a probability of 0.
                // If we didn't explicity skip them, then
                // they would actually overwrite the previous entry!
                if(prob != 0.0)
                    tempTable[cumProb] = seq;
            }
        }
        
        line = r.getLine();
    }

    // Normalize the table so that the last entry has a cumulative probability of 1.
    for(std::map<double, std::string>::iterator it = tempTable.begin(); it != tempTable.end(); ++it) {
        double prob = it->first / cumProb;
        seqTable[prob] = it->second;
    }
    
    r.close();
}

void RandomSeqGenerator::generateSeq(std::string &_seq, int _seqLength) {
    stringstream ss;
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);

    for(int i = 0; i < _seqLength; ++i) {
        std::map<double, string>::iterator it = seqTable.begin();
        double randNumber = uni();
        while(it->first <= randNumber ) {
            ++it;
        }
        
        ss << it->second;
    }
    
    _seq = ss.str();
}
