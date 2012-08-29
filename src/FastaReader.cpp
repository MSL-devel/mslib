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
#include "FastaReader.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"

static MslOut MSLOUT("FastaReader");

// BOOST Includes
#ifdef __BOOST__
#include <boost/regex.hpp>
#endif

FastaReader::FastaReader() : Reader(){
}

FastaReader::FastaReader(const string &_filename) : Reader(_filename){
}

FastaReader::FastaReader(const FastaReader &_Fastareader){
}

FastaReader::~FastaReader(){}


bool FastaReader::read(string _regex){

	if (!is_open()) {
		return false;
	}

	try { 
	       string headerExp = "^>([\\s\\S]+)$";
	       string name = "";
		while (!endOfFileTest()){
			string line = Reader::getLine();

			// Skip blank lines.
			if (line.size() == 0) {
			  continue;
			}

			// Storage for regular expression matches for this line
			vector<string> reTokens;
			
			// Need to parse name 
			if (MslTools::regex(line,headerExp,reTokens)){
			  name = MslTools::trim(reTokens[0]);
			  
			  continue;
			}
			
			if (_regex != "" && !MslTools::regex(name,_regex,reTokens)){
			  MSLOUT.stream() << "FastaReader::read() name("<<name<<") filtered due to regex("<<_regex<<")"<<endl;
			  continue;
			}

			// Any non-name or non-blank string should be a 
			map<string,string>::iterator it;
			it = sequences.find(name);
			if (it == sequences.end()){
			    sequences[name] = line;
			} else {
			  sequences[name] += line;
			  /*
			  int i = 1;
			  string tmpName = "";
			  while (i < 100) {
			    tmpName = MslTools::stringf("%s_%d",name.c_str(),i++);
			    it = sequences.find(tmpName);
			    if (it == sequences.end()){
			      sequences[tmpName] = line;
			      break;
			    }
			  } 
			  cerr << "Duplicate name in fasta: "<<name<<" adding as "<<tmpName<<endl;
			  */
			}
		}
	} catch(...){
	  cerr << "ERROR 9090 in FastaReader::read()\n";
	  return false;
	}

	return true;
}
string FastaReader::getPositionId(Chain &_ch, int _index, std::string _key){

  // get fasta seq.
  string fasta = getSequence(_key);

  if (fasta == ""){
    cerr << "Fasta is blank!"<<endl;
    return "";
  } 

  //cout << "Fasta is: "<<fasta<< ". from key: "<<_key<<endl;


  // Test if index is even a potential position
  if (fasta[_index] == '-'){
    return "";
  }

#ifdef __BOOST__
  // Find index of first non '-' character
  boost::regex re("^\\-*(\\S)\\S*$");


  boost::cmatch matches;
  if (!boost::regex_match(fasta.c_str(),matches,re)) {
    cerr << "ERROR 9234 with boost match"<<endl;
    exit(9234);
  }

  //cout << "Match: "<<matches.position()<< " fasta AA is: "<<fasta[matches.position()]<<endl;

  // Go from matches[0].first to posIndex OR end of string
  int i = matches.position();
  int chainIndex = 0;
  while (i != fasta.size()){

    if (i == _index){
      break;
    }
    
    if (fasta[i] != '-'){
      chainIndex++;
    }
    


    i++;
  }

  //cout << "ChainIndex is: "<<chainIndex<<" index="<<_index<< " in fasta: "<<fasta[_index]<<endl;
  if (i == _index){

    // Check it the position and this index have the same AA, if not then probably there is an issue.
    if (MslTools::getOneLetterCode(_ch.getPosition(chainIndex).getResidueName()) != MslTools::stringf("%c",fasta[_index])){
      cerr << "ERROR 9255 Position AA is: "<<MslTools::getOneLetterCode(_ch.getPosition(chainIndex).getResidueName())<< " and Fasta AA is: "<<fasta[_index]<<" index: "<<_index<<" chainIndex: "<<chainIndex<<endl;
      cerr << "Position.toString() = "<<_ch.getPosition(chainIndex).toString()<<" posId: "<<_ch.getPosition(chainIndex).getPositionId()<<endl;
      exit(9255);
    } 

    // Return the index
    return _ch.getPosition(chainIndex).getPositionId();
  }
#endif

  return "";
}

int FastaReader::getIndex(Chain &_ch, std::string _posId, std::string _key){

  Position &pos = _ch.getPosition(_posId);
  
  // get fasta seq.
  string fasta = getSequence(_key);

#ifdef __BOOST__
  int posIndex = _ch.getPositionIndex(&pos);

  // Find index of first non '-' character
  boost::regex re("^-*(!-)");

  boost::cmatch matches;
  if (!boost::regex_match(fasta.c_str(),matches,re)) {
    cerr << "ERROR 9234 with boost match"<<endl;
    exit(9234);
  }

  //cout << "Match: "<<matches.position()<< " fasta AA is: "<<fasta[matches.position()]<<endl;

  // Go from matches[0].first to posIndex OR end of string
  int index = matches.position();
  int posIndex_counter = posIndex;
  while (index != fasta.size()){
    
    if (fasta[index] != '-'){
      posIndex_counter--;
    }

    // Check if we are here!
    if (posIndex_counter == 0){
      break;
    }
    index++;
  }

  if (posIndex_counter == 0){

    // Check it the position and this index have the same AA, if not then probably there is an issue.
    if (MslTools::getOneLetterCode(pos.getResidueName()) != MslTools::stringf("%c",fasta[index])){
      cerr << "ERROR 9255 Position AA is: "<<MslTools::getOneLetterCode(pos.getResidueName())<< " and Fasta AA is: "<<fasta[index]<<endl;
      exit(9255);
    } 

    // Return the index
    return index;
  }
#endif

  return -1;
}
