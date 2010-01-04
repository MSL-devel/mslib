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
#include "RandomSeqGenerator.h"
#include "Reader.h"
#include "MslTools.h"
#include <sstream>
#include <ctime>
#include <boost/random.hpp>

using namespace std;

RandomSeqGenerator::RandomSeqGenerator(const RandomSeqGenerator & _rsg) {
    seqTable = _rsg.seqTable;
}

RandomSeqGenerator::RandomSeqGenerator(std::string &_seqTableFileName) {
    readSeqTable(_seqTableFileName);
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
    boost::mt19937 rng;
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng, uni_dist);
    rng.seed(static_cast<unsigned int>(std::time(0)));

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
