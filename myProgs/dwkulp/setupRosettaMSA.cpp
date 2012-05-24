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

  string structureSeq = fin.getSequence(opt.pdb);
  if (structureSeq == ""){
    cerr << "ERROR couldn't find "<<opt.pdb<<" as a name in : "<<opt.msa<<endl;
    exit(2342);
  }

  // Get "real" length (no dashes)
  int structureSeqLength = 0;
  for (uint i = 0; i < structureSeq.length();i++){
    if (structureSeq.substr(i,1) != "-"){
      structureSeqLength++;
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
    if (it->second == structureSeq) continue;	

    if (it->second.length() != structureSeq.length()){
      cerr << "ERROR sequence: "<<it->first<<" has "<<it->second.length()<<" residues, but structure has: "<<structureSeq.length()<<" , skipping"<<endl;
      continue;
    }

    // Create new files
    string basefile = MslTools::stringf("%s_model", MslTools::trim(it->first).c_str());
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
    for (uint i = 0; i < structureSeq.length();i++){
      
      // Look up character in structureSeq alignment 
      string wtAA = structureSeq.substr(i,1);

      // Look up character in current sequence
      string seqAA = it->second.substr(i,1);

      // Skip in both
      if (wtAA == "-" && seqAA == "-"){
	continue;
      }

      // Handle insertion into WT structure
      if (wtAA == "-" && seqAA != "-"){

	// Blueprint .. build this AA please.. go to structureIndex-1 and reset to remodel it...
	if (i > opt.remodel_neighbors-1){

	  for (uint r = opt.remodel_neighbors; r >= 1;r--){
	    if (structureSeq.substr(i-r,1) != "-" && alreadyRemodeled[i-r] != 1){
	      bfile[i-r] = MslTools::stringf("%d %1s D PIKAA %s", structureIndex-r+1,structureSeq.substr(i-r,1).c_str(),it->second.substr(i-r,1).c_str());
	      alreadyRemodeled[i-r] = 1;
	      cout << "RESET: "<<i<<" "<<structureIndex-r<< " R: "<<r<<" "<<structureSeq.substr(i-r,1)<<" :: "<<bfile[i-r]<<endl;
	    }
	  }

	  bfile[i] = MslTools::stringf("o x D PIKAA %s", seqAA.c_str());
	}
	// Resfile  .. skip

	continue;
      }

      // Handle deletion in WT
      if (seqAA == "-"){

	// Blueprint ... remove or ignore?
	
	// Resfile .. skip

	structureIndex++;
	continue;
      }


      // Both wtAA and seqAA != "-" here
      string remodelTag = ".";
      for (uint r = 1; r <= opt.remodel_neighbors;r++){
	if (i > opt.remodel_neighbors-1 && structureSeq.substr(i-r,1) == "-"){
	  remodelTag = MslTools::stringf("D PIKAA %s",seqAA.c_str());
	}
      }

        if (wtAA == seqAA){
  	// Blueprint PIKAA
	bfile[i] = MslTools::stringf("%d %1s %s", structureIndex+1,wtAA.c_str(),remodelTag.c_str());
      } else {
	// Blueprint PIKAA
	bfile[i] = MslTools::stringf("%d %1s %s PIKAA %s", structureIndex+1,wtAA.c_str(),remodelTag.c_str(), seqAA.c_str());
 
	// Resfile PIKAA
	rfile[i] = MslTools::stringf("%d%1s %1s PIKAA %s",sys.getPosition(structureIndex).getResidueNumber(),sys.getPosition(structureIndex).getResidueIcode().c_str(),sys.getPosition(structureIndex).getChainId().c_str(),seqAA.c_str());
      }


      // Add atoms..
      AtomContainer ats;
      ats.addAtoms(sys.getPosition(structureIndex).getAtomPointers());

      forRosettaFixbb.addAtoms(ats.getAtomPointers());

      for (uint a = 0; a < ats.size();a++){
	ats(a).setResidueNumber(structureIndex+1);
	ats(a).setChainId(" ");
	ats(a).setResidueIcode(" ");
      }
      
      forRosettaRemodel.addAtoms(ats.getAtomPointers());

      structureIndex++;

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
	opt.remodel_neighbors = OP.getInt("remodel_neighbors");
	if (OP.fail() || opt.remodel_neighbors < 1){
	  opt.remodel_neighbors = 1;
	}

	
	return opt;
}
