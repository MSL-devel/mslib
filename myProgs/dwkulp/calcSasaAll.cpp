#include <iostream>
#include <cstdlib>

#include "MslTools.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Timer.h"
#include "SasaCalculator.h"
#include "FastaReader.h"
#include "Position.h"
#include "calcSasaAll.h"


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("calcSasaAll");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnOn("calcSasaAll");

	Timer t;
	double start = t.getWallTime();

	// Read in a fasta-aligned file
	FastaReader fin;
	fin.open(opt.fasta);
	fin.read();
	fin.close();

	
	// Read in the pdb with a set of functional groups defined
	vector<string> list;
	MslTools::readTextFile(list, opt.list);

	System ref;
	ref.readPdb(list[0]);

	SasaCalculator sasaRef(ref.getAtomPointers(), 1.4, 200);
	sasaRef.calcSasa();

	map<string,string> seqs = fin.getSequences();
	
	string refSeq = seqs[opt.refSeqName];

	// first string is list[i] name, second is posId, double is tmpSasa - sasaRef.
	map<string, map<string, double> > sasaData;
	map<string, map<string, bool> > sasaResDiff;
	for (uint i = 1; i < list.size();i++){


	  System tmp;
	  tmp.readPdb(list[i]);
	  

	  SasaCalculator sasaCalc(tmp.getAtomPointers(), 1.4, 200);
	  sasaCalc.calcSasa();	  

	  // Get all sequences that have name in it...
	  vector<string> seqsFromThisPdb;
	  map<string,string>::iterator seqIt;
	  for (seqIt = seqs.begin();seqIt != seqs.end();seqIt++){
	    if (seqIt->first.substr(0,4) == MslTools::getFileName(list[i]).substr(0,4)){
	      seqsFromThisPdb.push_back(MslTools::trim(seqIt->first));
	    }
	  }
	  

	  // Get the reference sequence
	  // For each index in ref seq, get refPos and get tmpPos.
	  for (uint s = 0; s < refSeq.size();s++){
	    string refChain = opt.refSeqName.substr(opt.refSeqName.length()-1);
	    //MSLOUT.stream() << "Get position with chain: "<<refChain<<", index: "<<s<<" from seq: "<<opt.refSeqName<<endl;
	    string refPosId = fin.getPositionId(ref.getChain(refChain), s, opt.refSeqName);
	    Position &refPos = ref.getPosition(refPosId);
	    double refPosSasa = refPos.getSasa();



	    for (uint p = 0; p < seqsFromThisPdb.size();p++){

	      sasaData[refPosId][seqsFromThisPdb[p]] = 0.0; // default it to 0.0

	      string tmpChain = seqsFromThisPdb[p].substr(seqsFromThisPdb[p].length()-1);
	      //MSLOUT.stream() << "Get position with chain: "<<tmpChain<<", index: "<<s<<" from seq: "<<seqsFromThisPdb[p]<<endl;
	      string tmpPosId = fin.getPositionId(tmp.getChain(tmpChain), s, seqsFromThisPdb[p]);
	      if (tmpPosId == "") continue;

	      Position &tmpPos = tmp.getPosition(tmpPosId);
	      double tmpPosSasa = tmpPos.getSasa();	      


	      //MSLOUT.stream() << "Ref pos("<<refPos.getCurrentIdentity().getIdentityId()<<") = "<<seqsFromThisPdb[p]<<" ("<<tmpPos.getCurrentIdentity().getIdentityId()<<")";

	      double deltaSasa = refPosSasa - tmpPosSasa;
	      string tag = "";
	      if (refPos.getResidueName() != tmpPos.getResidueName()){
		sasaResDiff[refPosId][seqsFromThisPdb[p]] = true;
		tag = " **** ";

	      } else {
		sasaResDiff[refPosId][seqsFromThisPdb[p]] = false;
	      }
	      MSLOUT.fprintf(stdout, "%-20s: chain %s and resi %6d  === %-20s: chain %s and resi %6d : %8.3f %8.3f : %s\n",
			     opt.refSeqName.c_str(),refPos.getChainId().c_str(), refPos.getResidueNumber(),
			     seqsFromThisPdb[p].c_str(),tmpPos.getChainId().c_str(),tmpPos.getResidueNumber(),
			     refPosSasa, tmpPosSasa,tag.c_str());

	      sasaData[refPosId][seqsFromThisPdb[p]] = deltaSasa;
	     
	    }
	  }
	  
	} // FOR LIST[I]

	map<string,double>::iterator it;

	for (uint p = 0; p < ref.positionSize();p++){
	  Position &pos = ref.getPosition(p);

	  // Write labels...
	  if (p == 0){
	    fprintf(stdout, "          ");
	    for (it = sasaData[pos.getPositionId()].begin();it != sasaData[pos.getPositionId()].end();it++){
	      fprintf(stdout, " %8s  ", it->first.c_str());
	    }
	    fprintf(stdout, "\n");
	  }

	  // Write 
	  fprintf(stdout, "%8s ",pos.getPositionId().c_str());
	  for (it = sasaData[pos.getPositionId()].begin();it != sasaData[pos.getPositionId()].end();it++){
	    if (sasaResDiff[pos.getPositionId()][it->first]){
	      fprintf(stdout, " %8.3f* ",it->second);
	    }else{
	      fprintf(stdout, " %8.3f  ",it->second);
	    }
	  }
	  fprintf(stdout,"\n");
	}
	
	
	


}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "calcSasaAll --list LIST\n";
		exit(0);
	}
	opt.list= OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}
	opt.fasta = OP.getString("fasta");
	if (OP.fail()){
	  cerr << "ERROR 1111 fasta not specified.\n";
	  exit(1111);
	}
	opt.refSeqName = OP.getString("refSeqName");
	if (OP.fail()){
	  cerr << "ERROR 1111 refSeqName not specified.\n";
	  exit(1111);
	}
	return opt;
}
