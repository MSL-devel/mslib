#include <iostream>
#include <cstdlib>
#include <fstream>
#include "OptionParser.h"
#include "FastaReader.h"
#include "System.h"
#include "MslOut.h"
#include "PolymerSequence.h"

#include "setupRosettaMSA.h"


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("setupRosettaMSA");

int main(int argc, char *argv[]) {
  Options opt = setupOptions(argc,argv);

  FastaReader fin(opt.msa);
  fin.open();
  if (!fin.read()){
    cerr << "ERROR reading "<<opt.msa<<endl;
    exit(143);
  }
  fin.close();

  string structureSeq = fin.getSequence(MslTools::getFileName(opt.pdb));
  if (structureSeq == ""){
    cerr << "ERROR couldn't find "<<opt.pdb<<" as a name in : "<<opt.msa<<endl;
    exit(2342);
  }

  // Get "real" length (no dashes)
  bool firstNonDash = false;
  int firstNonDashChar = 0;
  int lastNonDashChar  = 0;
  int structureSeqLength = 0;
  for (uint i = 0; i < structureSeq.length();i++){
    if (structureSeq.substr(i,1) != "-"){
      structureSeqLength++;
      lastNonDashChar = i;
    }

    if (!firstNonDash && structureSeq.substr(i,1) != "-"){
      firstNonDash = true;
      firstNonDashChar = i;
    }
    

  }

  System sys;
  sys.readPdb(opt.pdb);
  if (structureSeqLength != sys.positionSize()){
    cerr << "ERROR Sequence in MSA for "<<opt.pdb<<" has "<<structureSeqLength<< " residues, but the structure has "<<sys.positionSize()<<endl;
    cerr << "From MSA: "<<structureSeq<<endl;
    cerr << "From SYS: "<<PolymerSequence::toOneLetterCode(sys.getAtomPointers())<<endl;
    exit(23433);
  }


  // Create a resfile(blueprint)/pdb for each of sequences in the msa..
  map<string,string> seqs = fin.getSequences();
  map<string,string>::iterator it;
  for (it = seqs.begin(); it != seqs.end();it++){
    if (it->first == opt.pdb) continue;	
    //if (it->second == structureSeq) continue;	
    bool skip=true;
    if (opt.select_seqs.size() == 0) skip = false;
    for (uint i = 0; i < opt.select_seqs.size() ;i++){
      if (MslTools::trim(it->first) == MslTools::trim(opt.select_seqs[i])){
	skip=false;
      }
    }
    if (it->first == opt.pdb || it->first == MslTools::getFileName(opt.pdb)) skip = true;

    if (skip) continue;

    cout << "Setting up: "<<it->first<<endl;
    if (it->second.length() != structureSeq.length()){
      cerr << "ERROR sequence: "<<it->first<<" has "<<it->second.length()<<" residues, but structure has: "<<structureSeq.length()<<" , skipping"<<endl;
      continue;
    }

    // Create new files
    string basefile = MslTools::stringf("%s_%s_model", MslTools::getFileName(opt.pdb).c_str(),MslTools::trim(it->first).c_str());
    string pdbfile_fixbb = MslTools::stringf("%s.pdb", basefile.c_str());
    string resfile = MslTools::stringf("%s.resfile", basefile.c_str());
    string pdbfile_remodel = MslTools::stringf("%s_remodel.pdb", basefile.c_str());
    string bluefile = MslTools::stringf("%s.blueprint", basefile.c_str());

    // Store lines in map so we can revert some of them before printing.
    map<int,string> rfile;
    map<int,string> bfile;

    // Need a pdb that may have deletions
    System forRosettaFixbb;

    // Need a renumbered pdb with deletions
    System forRosettaRemodel;

    map<int,int> alreadyRemodeled;
    // Go through each residue of the structure
    int structureIndex = 0;
    int remodelIndex = 0;
    int lastStructuredResidue = -1;
    int lastStructuredResidueIndex = -1;
    int lastInsert = 0;
    bool lastDelete = false;
    for (uint i = 0; i < structureSeq.length();i++){
      cout << "I: "<<i<<endl;
      // Skip initial '-' and ending '-'.
      if (opt.skipBlankEnds && (i < firstNonDashChar || i > lastNonDashChar) ){
	  continue;
      }
      lastDelete=false;

      // Look up character in structureSeq alignment 
      string wtAA = structureSeq.substr(i,1);

      // Look up character in current sequence
      string seqAA = it->second.substr(i,1);

      cout << "1. seqAA = '"<<seqAA<<"' and wtAA = '"<<wtAA<<"'"<<endl;
      if (wtAA != "-" && seqAA != "-"){
	lastStructuredResidue = remodelIndex+1;
	lastStructuredResidueIndex = i;
      }


      // Skip in both
      if (wtAA == "-" && seqAA == "-"){
	continue;
      }

      // Handle insertion into WT structure
      if (wtAA == "-" && seqAA != "-"){

	// Blueprint .. build this AA please.. go to structureIndex-1 and reset to remodel it...
	if (i > opt.remodel_neighbors-1 && lastStructuredResidue > opt.remodel_neighbors-1){

	  if (structureSeq.substr(i-1,1) != "-") {
	    // We need to only use the last residues that HAVE structure as remodel-able neighbors..
	    for (uint r = opt.remodel_neighbors; r >= 1;r--){
	      if (!alreadyRemodeled[i-r]){
		bfile[i-r] = MslTools::stringf("%d %1s D PIKAA %s", remodelIndex-r+1,structureSeq.substr(i-r,1).c_str(),it->second.substr(i-r,1).c_str());
		alreadyRemodeled[i-r] = 1;
		cout << "RESET: "<<i<<" "<<lastStructuredResidue<<" "<<structureIndex-r<< " R: "<<r<<" "<<structureSeq.substr(lastStructuredResidue-r,1)<<" :: "<<bfile[i-r]<<endl;
	      } else {
		cout << "ALREADY REMODELED: "<<i<<" "<<lastStructuredResidue<<" "<<structureIndex-r<< " R: "<<r<<" "<<structureSeq.substr(lastStructuredResidue-r,1)<<" :: "<<bfile[i-r]<<endl;
	      }
	    }
	  }

	  bfile[lastStructuredResidueIndex] = MslTools::stringf("%d %1s D PIKAA %s", remodelIndex,structureSeq.substr(lastStructuredResidueIndex,1).c_str(),it->second.substr(lastStructuredResidueIndex,1).c_str());
	  cout << "INSERTION! "<<wtAA<<" "<<seqAA<<" "<<structureSeq.substr(i-1,1)<<" "<<it->second.substr(i-1,1)<<endl;
	  bfile[i] = MslTools::stringf("o x D PIKAA %s", seqAA.c_str());
	  lastInsert = true;
	}
	// Resfile  .. skip

	continue;
      }

      // Handle deletion in WT, so wtAA != '-' at this point.
      /*
	Remove residue 6:
	4 I .
	5 L H PIKAA L
	7 G H PIKAA G
	8 A .

	So we need to skip all the seqAA = '-' and wtAA != '-'
	
       */
      // First deleted residue
      if (seqAA == "-" && i > 0 && it->second.substr(i-1,1) != "-" && structureSeq.substr(i-1,1) != "-"){
	cout << "DELETION! First residue."<<endl;

	// Remodel the number of remodel neighbors..
	for (uint r = opt.remodel_neighbors; r >= 1;r--){
	  bfile[i-r] = MslTools::stringf("%d %1s D PIKAA %s", remodelIndex+1-r,it->second.substr(i-r,1).c_str(),it->second.substr(i-r,1).c_str());
	}

	lastInsert=false;
	lastDelete=true;
      } else if (seqAA != "-" && i > 0 && it->second.substr(i-1,1) == "-" && structureSeq.substr(i-1,1) != "-"){       // Last deleted residue
	cout << "DELETION! last residue."<<endl;

	// Remodel the number of remodel neighbors..
	for (uint r = 0; r < opt.remodel_neighbors;r++){
	  cout << "r: "<<r<<endl;
	    bfile[i+r] = MslTools::stringf("%d %1s D PIKAA %s", remodelIndex+1+r,structureSeq.substr(i+r,1).c_str(),structureSeq.substr(i+r,1).c_str());
	    alreadyRemodeled[i+r] = 1;
	}

	lastInsert=false;
	lastDelete=true;
      } else if (seqAA == "-"){     // In the middle of a deleted region..
	cout << "DELETION! middle residue."<<endl;
	//structureIndex++;
	lastInsert = false;
	lastDelete = true;
      }


      // Both wtAA and seqAA != "-" here
      string remodelTag = MslTools::stringf(". PIKAA %s", seqAA.c_str());
      if (wtAA == seqAA){
	remodelTag = ".";
      }

      if (lastInsert){
	for (uint r = 1; r <= opt.remodel_neighbors;r++){
	  if (i > opt.remodel_neighbors-1 && structureSeq.substr(i-r,1) == "-"){
	    remodelTag = MslTools::stringf("D PIKAA %s",seqAA.c_str());
	  }
	}
      }


      if (!lastDelete){
	if (wtAA == seqAA){

	  // Blueprint PIKAA
	  if (alreadyRemodeled[i] != 1){
	    bfile[i] = MslTools::stringf("%d %1s %s", remodelIndex+1,wtAA.c_str(),remodelTag.c_str());
	  }
	} else {
	 
	  // Blueprint PIKAA
	  if (alreadyRemodeled[i] != 1){
	    bfile[i] = MslTools::stringf("%d %1s %s", remodelIndex+1,wtAA.c_str(),remodelTag.c_str());
	  }
 
	  // Resfile PIKAA
	  rfile[i] = MslTools::stringf("%d%1s %1s PIKAA %s EX 1 LEVEL 4 EX 2 LEVEL 4 EX 3 EX 4",sys.getPosition(structureIndex).getResidueNumber(),sys.getPosition(structureIndex).getResidueIcode().c_str(),sys.getPosition(structureIndex).getChainId().c_str(),seqAA.c_str());
	}
      }

      lastInsert = false;

      // Add atoms..
      AtomContainer ats;
      ats.addAtoms(sys.getPosition(structureIndex).getAtomPointers());

      forRosettaFixbb.addAtoms(ats.getAtomPointers());

      cout << "Reset residue number: "<<sys.getPosition(structureIndex).toString()<<" to "<<remodelIndex+1<<endl;
      for (uint a = 0; a < ats.size();a++){
	ats(a).setResidueNumber(remodelIndex+1);
	ats(a).setChainId(" ");
	ats(a).setResidueIcode(" ");
      }
      
      forRosettaRemodel.addAtoms(ats.getAtomPointers());

      structureIndex++;
      remodelIndex++;

    }
   

    // Print out the files
    fstream fresfile;
    fresfile.open(resfile.c_str(),std::ios::out);
    fresfile << "NATRO" <<endl;
    fresfile << "start"<<endl;
    map<int,string>::iterator fit;
    for (fit = rfile.begin();fit != rfile.end();fit++){
      fresfile << fit->second <<endl;
    }
    fresfile.close();


    fstream fbluefile;
    fbluefile.open(bluefile.c_str(),std::ios::out);
    for (fit = bfile.begin();fit != bfile.end();fit++){
      fbluefile << fit->second <<endl;
    }
    fbluefile.close(); 


    // Write out PDBs
    forRosettaFixbb.writePdb(pdbfile_fixbb);
    forRosettaRemodel.writePdb(pdbfile_remodel);
    
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
		cout << "setupRosettaMSA --msa msa.fasta --pdb PDB\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.msa = OP.getString("msa");
	if (OP.fail()){
		cerr << "ERROR 1111 msa not specified.\n";
		exit(1111);
	}
	opt.select_seqs = OP.getMultiString("select_seqs");
	if (OP.fail()){
		cerr << "WARNING select_seqs not specified.\n";
	}
	opt.remodel_neighbors = OP.getInt("remodel_neighbors");
	if (OP.fail() || opt.remodel_neighbors < 1){
	  opt.remodel_neighbors = 1;
	}
	
	opt.skipBlankEnds = OP.getBool("skipBlankEnds");
	if (OP.fail()){
	  opt.skipBlankEnds = false;
	}
	
	return opt;
}
