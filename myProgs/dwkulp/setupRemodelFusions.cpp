#include <iostream>
#include <cstdlib>
#include "System.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include "MslOut.h"


#include "setupRemodelFusions.h"



using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("setupRemodelFusions");


// Utility functions
void sighandler(int sig);


int main(int argc, char *argv[]) {
  
  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // MslOut can suppress output, to make example output clean
  MSLOUT.turnOn("setupRemodelFusions");

  // Read in file
  System sys;
  sys.readPdb(opt.pdb);

  // Make sure insertion points exist in system.
  if (!sys.positionExists(opt.insertionPoint1domain1)){
    cerr << "ERROR 2222 insertionPoint1domain1("<<opt.insertionPoint1domain1<<") does not exist in "<<opt.pdb<<endl;
    exit(2222);
  }
  if (!sys.positionExists(opt.insertionPoint1domain2)){
    cerr << "ERROR 2222 insertionPoint1domain2("<<opt.insertionPoint1domain2<<") does not exist in "<<opt.pdb<<endl;
    exit(2222);
  }

  if (opt.insertionPoint2domain1 != "" && !sys.positionExists(opt.insertionPoint2domain1)){
    cerr << "ERROR 2222 insertionPoint2domain1("<<opt.insertionPoint2domain1<<") does not exist in "<<opt.pdb<<endl;
    exit(2222);
  }
  if (opt.insertionPoint2domain2 != "" && !sys.positionExists(opt.insertionPoint2domain2)){
    cerr << "ERROR 2222 insertionPoint2domain2("<<opt.insertionPoint2domain2<<") does not exist in "<<opt.pdb<<endl;
    exit(2222);
  }

  map<string,bool> extraDesignPositions;
  map<string,bool>::iterator extraIt;
  if (opt.extraDesignPositions != ""){
      vector<string> extraPos = MslTools::tokenize(opt.extraDesignPositions,":");
      for (uint i = 0; i < extraPos.size();i++){
	extraDesignPositions[extraPos[i]] = true;
	MSLOUT.stream() << "Extra design position is: "<<extraPos[i]<<endl;
      }
  }

  map<string,bool> extraRemodelPositions;
  map<string,bool>::iterator extraRemodelIt;
  if (opt.extraRemodelPositions != ""){
      vector<string> extraPos = MslTools::tokenize(opt.extraRemodelPositions,":");
      for (uint i = 0; i < extraPos.size();i++){
	extraRemodelPositions[extraPos[i]] = true;
	MSLOUT.stream() << "Extra remodel position is: "<<extraPos[i]<<endl;
      }
  }
  
  // Blueprint file (s)
  vector<string> linker1points = MslTools::tokenize(opt.linker1lengths,",");
  vector<string> linker2points = MslTools::tokenize(opt.linker2lengths,",");
  MSLOUT.stream() << "Linker sizes: "<<linker1points.size()<<" "<<linker2points.size()<<endl;
  map<int,map<int,ofstream *> > blueprint_files;
  map<int,map<int,ofstream *> >::iterator it1;
  map<int,ofstream *>::iterator it2;
  for (uint i = 0; i < linker1points.size(); i++){
    int length1 = MslTools::toInt(linker1points[i]);
    for (uint j = 0; j < linker2points.size(); j++){
      int length2 = MslTools::toInt(linker2points[j]);
      string filename = MslTools::stringf("blueprint_%s_%02d_%02d", MslTools::getFileName(opt.pdb).c_str(),length1,length2);
      ofstream *fs = new ofstream();
      fs->open(filename.c_str());
      blueprint_files[length1][length2] = fs;
    }
  }


  // Create different segments (part1domain1, part1domain2, part2domain1,part2domain2)  
  string domain1_chainid;
  int tmp_resnum;
  string tmp_icode;
  MslTools::parsePositionId(opt.insertionPoint1domain1, domain1_chainid, tmp_resnum, tmp_icode);

  string domain2_chainid;
  MslTools::parsePositionId(opt.insertionPoint1domain2, domain2_chainid, tmp_resnum, tmp_icode);

  MSLOUT.stream() << "Domain1 = "<<domain1_chainid<<" "<<opt.insertionPoint1domain1<<endl;
  MSLOUT.stream() << "Domain2 = "<<domain2_chainid<<" "<<opt.insertionPoint1domain2<<endl;


  // Design the interface
  map<string, string> designPosition;
  map<string,string>::iterator designIt;
  if (opt.designInterfaceDistance != 0.0){
    for (uint i = 0; i < sys.getChain(domain1_chainid).positionSize();i++){
      Position &pos1 = sys.getChain(domain1_chainid).getPosition(i);
      for (uint j = 0; j < sys.getChain(domain2_chainid).positionSize();j++){
	Position &pos2 = sys.getChain(domain2_chainid).getPosition(j);

	if (!pos1.atomExists("CA") || !pos2.atomExists("CA")) continue;

	
	double dist = pos1.getAtom("CA").distance(pos2.getAtom("CA"));
	if (dist <= opt.designInterfaceDistance){
	  designIt = designPosition.find(pos1.getPositionId());
	  if (designIt == designPosition.end()){
	    MSLOUT.stream() << MslTools::stringf("Designing position: %s\n",pos1.getPositionId().c_str());
	  }

	  designIt = designPosition.find(pos2.getPositionId());
	  if (designIt == designPosition.end()){
	    MSLOUT.stream() << MslTools::stringf("Designing position: %s\n",pos2.getPositionId().c_str());
	  }

	  designPosition[pos1.getPositionId()] = "ALLAA";
	  designPosition[pos2.getPositionId()] = "ALLAA";
	}

      }
    }
  }


  System newProtein;
  int resNum = 1;
  int startIndex = sys.getChain(domain1_chainid).getPosition(0).getIndexInSystem();
  int endIndex = sys.getPositionIndex(opt.insertionPoint1domain1);
  MSLOUT.stream() << MslTools::stringf("%s(%-4d) to %s(%-4d)\n",sys.getChain(domain1_chainid).getPosition(0).getPositionId().c_str(),startIndex, opt.insertionPoint1domain1.c_str(),endIndex);
  // part1_domain1 = resi 1 of domain1_chainid to insertionPoint1domain1
  for (uint i = startIndex ; i <= endIndex;i++){
    Position &pos = sys.getPosition(i);

    // Add this to part
    AtomContainer tmpAts;
    tmpAts.addAtoms(pos.getAtomPointers());
    
    for (uint a = 0; a < tmpAts.size();a++){
      tmpAts(a).setResidueNumber(resNum);
      tmpAts(a).setChainId(" ");
    }

    newProtein.addAtoms(tmpAts.getAtomPointers());
    
    string remodeltag = ".";
    designIt        = designPosition.find(pos.getPositionId());
    extraIt         = extraDesignPositions.find(pos.getPositionId());
    extraRemodelIt  = extraRemodelPositions.find(pos.getPositionId());
    if (i == endIndex) {
      remodeltag = "D ALLAA";
    }else if (designIt != designPosition.end()){
      remodeltag = ". ALLAA";
    }else if (extraIt != extraDesignPositions.end()){
      remodeltag = ". ALLAA";
    } else if (extraRemodelIt != extraRemodelPositions.end()){
      remodeltag = "D ALLAA";
    }



    
    for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
      for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
	*(it2->second) << MslTools::stringf("%-4d %1s %s\n", resNum, MslTools::getOneLetterCode(pos.getResidueName()).c_str(), remodeltag.c_str());
      }
    }

    resNum++;
  }


  // Insertion of residue to blueprint files

  for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
    for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
      for (uint i = 1; i <= it1->first;i++){
	(*it2->second) << MslTools::stringf("0 x L ALLAA\n");
      }
    }
  }  

  // domain2 = insertionPoint1domain2 to insertionPoint2domain2
  startIndex = sys.getPosition(opt.insertionPoint1domain2).getIndexInSystem();
  endIndex   = sys.getPosition(opt.insertionPoint2domain2).getIndexInSystem();
  MSLOUT.stream() << MslTools::stringf("%s(%-4d) to %s(%-4d)\n",opt.insertionPoint1domain2.c_str(),startIndex, opt.insertionPoint2domain2.c_str(),endIndex);
  for (uint i = startIndex ; i <= endIndex;i++){
    Position &pos = sys.getPosition(i);

    // Add this to part
    AtomContainer tmpAts;
    tmpAts.addAtoms(pos.getAtomPointers());
    
    for (uint a = 0; a < tmpAts.size();a++){
      tmpAts(a).setResidueNumber(resNum);
      tmpAts(a).setChainId(" ");
    }

    newProtein.addAtoms(tmpAts.getAtomPointers());

    string remodeltag = ".";
    designIt = designPosition.find(pos.getPositionId());
    extraIt  = extraDesignPositions.find(pos.getPositionId());
    extraRemodelIt  = extraRemodelPositions.find(pos.getPositionId());
    if ( i == startIndex || (opt.insertionPoint2domain1 != "" && i == endIndex)) {
      remodeltag = "D ALLAA";
    } else if (designIt != designPosition.end()){
      remodeltag = ". ALLAA";
    } else if (extraIt != extraDesignPositions.end()){
      remodeltag = ". ALLAA";
    } else if (extraRemodelIt != extraRemodelPositions.end()){
      remodeltag = "D ALLAA";
    } 

    for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
      for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
	(*it2->second) << MslTools::stringf("%-4d %1s %s\n", resNum, MslTools::getOneLetterCode(pos.getResidueName()).c_str(), remodeltag.c_str());
      }
    }

    resNum++;

  }

  if (opt.insertionPoint2domain1 != ""){

    // Insertion of residue to blueprint files
    for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
      for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
	for (uint i = 1; i <= it2->first;i++){
	  (*it2->second) << MslTools::stringf("0 x L ALLAA\n");
	}
      }
    }    


  
    // part2_domain1 = insertionPoint2domain1 to end_domain1
    startIndex = sys.getPosition(opt.insertionPoint2domain1).getIndexInSystem();
    endIndex = sys.getChain(domain1_chainid).getPosition(sys.getChain(domain1_chainid).positionSize()-1).getIndexInSystem();
    MSLOUT.stream() << MslTools::stringf("%s(%-4d) to %s(%-4d)\n",opt.insertionPoint2domain1.c_str(),startIndex, sys.getChain(domain1_chainid).getPosition(sys.getChain(domain1_chainid).positionSize()-1).getPositionId().c_str(),endIndex);
    for (uint i = startIndex ; i <= endIndex ;i++){
      Position &pos = sys.getPosition(i);

      // Add this to part
      AtomContainer tmpAts;
      tmpAts.addAtoms(pos.getAtomPointers());
    
      for (uint a = 0; a < tmpAts.size();a++){
	tmpAts(a).setResidueNumber(resNum);
	tmpAts(a).setChainId(" ");
      }

      newProtein.addAtoms(tmpAts.getAtomPointers());

      string remodeltag = ".";
      designIt = designPosition.find(pos.getPositionId());
      extraIt  = extraDesignPositions.find(pos.getPositionId());
      extraRemodelIt  = extraRemodelPositions.find(pos.getPositionId());
      if (i == startIndex) {
	remodeltag = "D ALLAA";
      } else if (designIt != designPosition.end()){
	remodeltag = ". ALLAA";
      }else if (extraIt != extraDesignPositions.end()){
	remodeltag = ". ALLAA";
      } else if (extraRemodelIt != extraRemodelPositions.end()){
	remodeltag = "D ALLAA";
      }

      for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
	for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
	  (*it2->second) << MslTools::stringf("%-4d %1s %s\n", resNum, MslTools::getOneLetterCode(pos.getResidueName()).c_str(),remodeltag.c_str());
	}
      }

      resNum++;

    }
  }


  // Add extra chains...
  for (uint c = 0;c < sys.chainSize();c++){
    if (sys.getChain(c).getChainId() == domain1_chainid ||
	sys.getChain(c).getChainId() == domain2_chainid) continue;

    MSLOUT.stream() << "Adding "<<sys.getChain(c).getChainId()<<" "<<sys.getChain(c).getPosition(0)<< " "<<sys.getChain(c).getPosition(sys.getChain(c).positionSize()-1)<<endl;
    for (uint i = 0; i < sys.getChain(c).positionSize();i++){

      Position &pos = sys.getChain(c).getPosition(i);

      // Add this to part
      AtomContainer tmpAts;
      tmpAts.addAtoms(pos.getAtomPointers());
    
      for (uint a = 0; a < tmpAts.size();a++){
	tmpAts(a).setResidueNumber(resNum);
	tmpAts(a).setChainId(" ");
      }

      newProtein.addAtoms(tmpAts.getAtomPointers());

      string remodeltag = ".";
      designIt = designPosition.find(pos.getPositionId());
      extraIt  = extraDesignPositions.find(pos.getPositionId());
      extraRemodelIt  = extraRemodelPositions.find(pos.getPositionId());
      if (designIt != designPosition.end()){
	remodeltag = ". ALLAA";
      } else if (extraIt != extraDesignPositions.end()){
	remodeltag = ". ALLAA";
      } else if (extraRemodelIt != extraRemodelPositions.end()){
	remodeltag = "D ALLAA";
      }
      for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
	for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
	  (*it2->second) << MslTools::stringf("%-4d %1s %s\n", resNum, MslTools::getOneLetterCode(pos.getResidueName()).c_str(),remodeltag.c_str());
	}
      }
      resNum++;
    }
  }

  // Close the files
  for (it1 = blueprint_files.begin(); it1 != blueprint_files.end();it1++){
    for (it2 = it1->second.begin(); it2 != it1->second.end();it2++){
      (*it2->second).close();
    }
  }


  newProtein.writePdb(MslTools::stringf("%s_remodel.pdb",MslTools::getFileName(opt.pdb).c_str()));
  MSLOUT.stream() << "Number of residues: "<<resNum<<"  . Done."<<endl;

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
    cout << "setupRemodelFusions --pdb PDB --insertionPoint1domain1 --insertionPoint1domain2 --insertionPoint2domain1 --insertionPoint2domain2\n";
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


  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 ref not specified.\n";
    exit(1111);
  }

  opt.insertionPoint1domain1 = OP.getString("insertionPoint1domain1");
  if (OP.fail()){
    cerr << "ERROR 1111 insertionPoint1domain1 not specified.\n";
    exit(1111);
  }

  opt.insertionPoint1domain2 = OP.getString("insertionPoint1domain2");
  if (OP.fail()){
    cerr << "ERROR 1111 insertionPoint1domain2 not specified.\n";
    exit(1111);
  }

  opt.insertionPoint2domain1 = OP.getString("insertionPoint2domain1");
  if (OP.fail()){
    opt.insertionPoint2domain1 = "";
  }

  opt.insertionPoint2domain2 = OP.getString("insertionPoint2domain2");
  if (OP.fail()){
    cerr << "ERROR 1111 insertionPoint2domain2 not specified.\n";
    exit(1111);
  }

  opt.linker1lengths = OP.getString("linker1lengths");
  if (OP.fail()){
    cerr << "ERROR 1111 linker1lengths not specified.\n";
    exit(1111);
  }
  opt.linker2lengths = OP.getString("linker2lengths");
  if (OP.fail()){
    cerr << "ERROR 1111 linker2lengths not specified.\n";
    exit(1111);
  }
  opt.designInterfaceDistance = OP.getDouble("designInterfaceDistance");
  if (OP.fail()){
    opt.designInterfaceDistance = 0.0;
  }

  opt.extraDesignPositions = OP.getString("extraDesignPositions");
  opt.extraRemodelPositions = OP.getString("extraRemodelPositions");

  return opt;
}
