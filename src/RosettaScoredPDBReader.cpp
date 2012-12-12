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


#include "RosettaScoredPDBReader.h"
#include "PDBFormat.h"

using namespace MSL;
using namespace std;

/**
 * Simple constructor.
 */
RosettaScoredPDBReader::RosettaScoredPDBReader() : Reader() {
}

RosettaScoredPDBReader::RosettaScoredPDBReader(const std::string &_filename) : Reader(_filename) {
}


RosettaScoredPDBReader::RosettaScoredPDBReader(const RosettaScoredPDBReader &_reader) {

  //residueScores = _reader.getResidueScores();
}


bool RosettaScoredPDBReader::read(){


	if (!is_open()) {
		return false;
	}


	try { 

	        bool begin_score_section = false;
		string currentResidue = "";

		while (!endOfFileTest()){
		  string line = Reader::getLine();
		  if (line.size() == 0) continue;

		  //cout << "LINE: "<<line<<endl;

		  string header = "";
		  if (line.size() >= PDBFormat::S_RECORD_NAME + PDBFormat::L_RECORD_NAME) {
		    header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);


		    if (header == "ATOM  "){
		      PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line);
		    /*
		      atom.D_ATOM_NAME;
		      atom.D_RES_NAME;
		      atom.D_I_CODE;
		      atom.D_RES_SEQ;
		      atom.D_CHAIN_ID;
		    */

		      if (strcmp(atom.D_ATOM_NAME, "CA") == 0){

			string posId = MslTools::getPositionId(atom.D_CHAIN_ID,atom.D_RES_SEQ,atom.D_I_CODE);
			aaSeq.push_back(pair<string,string>(atom.D_RES_NAME,posId));
			//cout << "aaSeq["<<aaSeq.size()-1<<"] is "<<atom.D_RES_NAME<<" and posId: "<<posId<<endl;
		      }
		      continue;
		    } // END if HEADER == "ATOM"
		  } // END HEADER
		

		  // Tokenize line
		  vector<string> toks = MslTools::tokenize(line);

		  if (toks[0] == "#BEGIN_POSE_ENERGIES_TABLE"){
		    //cout << "BEGIN SCORE SECTION!!!!!!!!!!"<<endl;
		    begin_score_section = true;
		    continue;
		  }

		  if (toks[0] == "#END_POSE_ENERGIES_TABLE"){
		    break;
		  }

		  if (toks[0] == "weights" || toks[0] == "pose") continue;

		  if (!begin_score_section) continue;

		  if (toks[0] == "label"){
		    score_labels = toks;
		    continue;
		  }

		  boost::regex expression("[A-Z]{3}_");
		  boost::sregex_iterator r1(toks[0].begin(),toks[0].end(),expression);
		  boost::sregex_iterator r2;
		  string resAA = "";
		  string resNum = "";
		  if (r1 == r2){
		    cout << "NO REGEX MATCHES!!!!!"<<endl;
		  } else {

		    resAA  = toks[0].substr((*r1).position(),(*r1).position()+(*r1).length()-1);
		    boost::regex expression("_[0-9]+$");
		    boost::sregex_iterator r1(toks[0].begin(),toks[0].end(),expression);
		    boost::sregex_iterator r2;
		    if (r1 == r2){
		      cout << "NO REGEX MATCHES!!!!!"<<endl;
		    } else {
		      resNum = toks[0].substr((*r1).position()+1,(*r1).position()+(*r1).length());
		    }
		  }



		  int index = MslTools::toInt(resNum);
		  if (index > aaSeq.size()){
		    cerr << "ERROR 112211 RosettaScoredPDBReader::read()...Index is "<<index<<" number of amino acids in PDB portion of file is: "<<aaSeq.size()<<endl;
		    exit(112211);
		  }
		  //fprintf(stdout,"toks[0] = %s, ResAA: %s, ResNum: %s.....",toks[0].c_str(),resAA.c_str(),resNum.c_str());
		  //fprintf(stdout,"%4d, AA = %s, PosId = %s\n",index-1,aaSeq[index-1].first.c_str(),aaSeq[index-1].second.c_str());
		  
		  // tokenize rest of line
		  map<string,double> scores;
		  for (uint i = 1; i < toks.size();i++){
		    double score = MslTools::toDouble(toks[i]);
		    //fprintf(stdout, "Score: %s = %8.3f\n",score_labels[i].c_str(),score);
		    scores[score_labels[i]] = score;

		  }
		  residueScores[aaSeq[index-1].second] = scores;


		} // WHILE NOT END OF FILE

	} // END try block 
	catch(...){
		cerr << "ERROR 3265 in RosettaScoredPDBReader::read()\n";
		return false;
	}

	return true;
}


