#include <iostream>
#include <cstdlib>
#include <queue>

#include "System.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "PSSMCreator.h"
#include "PolymerSequence.h"
#include "AtomSelection.h"

#include "addSequenceConservation.h"

#include <queue>

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("addSequenceConservation");


priority_queue<pair<double,string> > orderMap(map<string,double> &_aMap);

int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);


  // Read multiple sequence alignment and create a PSSM
  PSSMCreator pssm;
  MSLOUT.stream() << "Add FASTA to PSSM"<<endl;
  pssm.addMultipleSequenceAlignment(opt.fasta,opt.regexForFasta);
  MSLOUT.stream() << "Add refAACounts to PSSM"<<endl;
  pssm.readReferenceCounts(opt.refAACounts);
  MSLOUT.stream() << "Create the PSSM"<<endl;
  if (opt.freq){
    pssm.createUsingMultipleSequenceAlignment(PSSMCreator::freq);
  } else {
    pssm.createUsingMultipleSequenceAlignment(PSSMCreator::logodds);
  }


  if (opt.pdb != ""){
    // Read in the reference pdb file
    System sys;
    sys.readPdb(opt.pdb);
  
    AtomSelection select(sys.getAtomPointers());
    AtomPointerVector ats = select.select(opt.sel);
    string seq = PolymerSequence::toOneLetterCode(ats);

    MSLOUT.stream() << "Getting scoring function for sequence: "<<seq<<endl;

    MSLOUT.stream() << "refSeqName is: "<<opt.refSeqName<<endl;
    vector<double> scores = pssm.getScoreFunction(seq,opt.refSeqName,opt.refSeqOffset);

    System selSys;
    selSys.addAtoms(ats);
    if (selSys.positionSize() != scores.size()){
      cerr << "ERROR 234 Selection position size is : "<<selSys.positionSize()<< " score vector is size: "<<scores.size()<<endl;
      exit(234);
    }

    // Clear all temp-factors
    for (uint i = 0; i < sys.positionSize();i++){
      Position &pos = sys.getPosition(i);
      for (uint a = 0; a < pos.atomSize();a++){
	pos.getAtom(a).setTempFactor(0.0);
      }      
    }    

    // Re-assign with score
    for (uint i = 0; i < selSys.positionSize();i++){

      Position &pos = sys.getPosition(selSys.getPosition(i).getPositionId());      
      for (uint c = 0; c < sys.chainSize();c++){
	
	if (!(opt.applyToAllChains || sys.getChain(c).getChainId() == pos.getChainId())) continue;

	Position &chainPos = sys.getPosition(MslTools::stringf("%s,%d%s",sys.getChain(c).getChainId().c_str(),pos.getResidueNumber(),pos.getResidueIcode().c_str()));
	for (uint a = 0; a < chainPos.atomSize();a++){
	  chainPos.getAtom(a).setTempFactor(scores[i]);
	}

      }
      
      
    }

    sys.writePdb(opt.outPdb);
  } // IF opt.pdb

  if (opt.seq != ""){
    vector<string> toks = MslTools::tokenize(opt.seq, ",");
    
    if (toks.size() != 3){
      cerr << "ERROR opt.seq needs to be in the format REFNAME,STARTINDEX,ENDINDEX and it is: "<<opt.seq<<endl;
    }
    string refName = toks[0];
    int start      = MslTools::toInt(toks[1]);
    int end        = MslTools::toInt(toks[2]);

    vector<map<string,double> > freqs = pssm.getFrequencies(refName,start,end);
    string seq = pssm.getSequence(refName,start,end);


    // ******   CRAZY FORMATTING CODE.... THIS SUCKS ***** //

    // For each position
    int largestVariation = 0;
    for (uint i = 0;i < freqs.size();i++){
      if (freqs[i].size() > largestVariation){
	largestVariation = freqs[i].size();
      }
    }

    int resNum = start;

    vector<vector<string> > output;
    for (uint i = 0;i < freqs.size();i++){

      vector<string> rows;
      priority_queue<pair<double,string> > orderedMap = orderMap(freqs[i]);
      int blankRows = largestVariation - orderedMap.size();
      for (uint v = 0; v < largestVariation;v++){
	if (v <= blankRows) {
	  rows.push_back(MslTools::stringf("      "));
	  continue;
	}
	
	while (orderedMap.size() > 0){
	  pair<double,string> data = orderedMap.top();
	  orderedMap.pop();
	  rows.push_back(MslTools::stringf("%1s %3d ", data.second.c_str(), data.first));
        }
      }
	
      // Add residue number at very end
      rows.push_back(MslTools::stringf(" %4d ",resNum++));
      output.push_back(rows);
    }

    for (uint r = 0; r < output.size();r++){
      string outrow = "";
      for (uint c = 0; c < output[r].size();c++){
	outrow += output[r][c];
      }
      fprintf(stdout,"%s\n",outrow.c_str());
    }

  }

  MSLOUT.stream() << "Done."<<endl;
}

Options setupOptions(int theArgc, char * theArgv[]){
  Options opt;

  OptionParser OP;


  OP.setRequired(opt.required);
  OP.setAllowed(opt.optional);
  OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
  OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
  OP.readArgv(theArgc, theArgv);

  if (OP.countOptions() == 0){
    cout << "Usage:" << endl;
    cout << endl;
    cout << "addSequenceConservation --pdb PDB --fasta FASTA_FILE --refCounts REF_COUNTS_FILE\n";
    exit(0);
  }


  opt.configfile = OP.getString("configfile");
  if (opt.configfile != "") {
    OP.readFile(opt.configfile);
    if (OP.fail()) {
      cerr << "ERROR couldn't read : "<<opt.configfile<<endl;
      exit(1);
    }
  }

  opt.fasta = OP.getString("fasta");
  if (OP.fail()){
    cerr << "ERROR 1111 fasta not specified.\n";
    exit(1111);
  }

  opt.refSeqName = OP.getString("refSeqName");
  if (OP.fail()){
    cerr << "WARNING 1111 refSeqName not specified.\n";
    opt.refSeqName = "";
  } 

  opt.seq = OP.getString("seq");
  if (OP.fail()){
    opt.seq = "";
  }

  if (opt.seq == ""){
    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
      cerr << "ERROR 1111 pdb not specified.\n";
      exit(1111);
    }


    opt.refAACounts = OP.getString("refAACounts");
    if (OP.fail()){
      cerr << "ERROR 1111 refAACounts not specified.\n";
      exit(1111);
    }

    opt.sel = OP.getString("sel");
    if (OP.fail()){
      cerr << "ERROR 1111 sel not specified.\n";
      exit(1111);
    }




    opt.refSeqOffset = OP.getInt("refSeqOffset");
    if (OP.fail()){
      cerr << "WARNING 1111 refSeqOffset not specified.\n";
      opt.refSeqOffset = 0;
    }
    opt.freq = OP.getBool("freq");
    if (OP.fail()){
      opt.freq = false;
    }

    opt.applyToAllChains = OP.getBool("applyToAllChains");
    if (OP.fail()){
      opt.applyToAllChains = false;
    }

    opt.regexForFasta = OP.getString("regexForFasta");
    if (OP.fail()){
      opt.regexForFasta = "";
    }

    opt.outPdb = OP.getString("outPdb");
    if (OP.fail()){
      opt.outPdb = MslTools::stringf("%s_seqcons.pdb",MslTools::getFileName(opt.pdb).c_str());
      cerr << "WARNING --outPdb is set to : "<<opt.outPdb<<endl;
    }

  }



  return opt;
}


priority_queue<pair<double,string> > orderMap(map<string,double> &_aMap){

  priority_queue<pair<double,string> > results;
  map<string,double>::iterator it;
  for (it = _aMap.begin();it != _aMap.end();it++){
    results.push(pair<double,string>(it->second,it->first));
  }

  return results;
}
