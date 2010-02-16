#include "ALNReader.h"


ALNReader::ALNReader() : Reader(){
  version = "";
  remark  = "";
  
}

ALNReader::ALNReader(const string &_filename) : Reader(_filename){
}

ALNReader::ALNReader(const ALNReader &_alnreader){
}

ALNReader::~ALNReader(){}


bool ALNReader::read(){

	try { 
	       string headerExp = "^CLUSTAL W (\\S+)\\s+([\\S\\s]+)$";
	       string dataExp   = "^(\\S+)\\s+(\\S+)\\s([0-9]+)$";
	       string equivExp  = "^([\\s\\S]+)$";
	       int seqIndex = 0;
		while (!endOfFileTest()){
			string line = Reader::getLine();

			// Skip blank lines.
			if (line.size() == 0) {
			  continue;
			}

			// Storage for regular expression matches for this line
			vector<string> reTokens;
			
			// Need to parse header line (CLUSTAL W 2.1 ...)
			if (MslTools::regex(line,headerExp,reTokens)){
			  
			  version = reTokens[0];
			  remark  = reTokens[1];
			  continue;
			}
			
			// Parse data lines
			if (MslTools::regex(line,dataExp,reTokens)){

			  map<string,string>::iterator it;
			  it = sequences.find(reTokens[0]);
			  if (it == sequences.end()){
			    sequences[reTokens[0]] = reTokens[1];
			  } else {
			    sequences[reTokens[0]] += reTokens[1];
			  }

			  // Find index of first character not key+spaces
			  seqIndex = line.find(reTokens[1]);
			  
			  
			  continue;
			}

			// Parse Residue equivalency lines
			if (MslTools::regex(line.substr(seqIndex),equivExp,reTokens)){
			  
			  map<string,string>::iterator it;
			  it = sequences.find("equivRes");
			  if (it == sequences.end()){
			    sequences["equivRes"] = reTokens[0];
			  } else {
			    sequences["equivRes"] += reTokens[0];
			  }

			  continue;
			}

			



		}
	} catch(...){
	  cerr << "ERROR 9090 in ALNReader::read()\n";
	  return false;
	}


	return true;
}
