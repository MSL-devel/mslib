#include <iostream>
#include <cstdlib>
#include <fstream>
#include <queue>
#include "OptionParser.h"
#include "FastaReader.h"
#include "System.h"
#include "MslOut.h"
#include "PolymerSequence.h"
#include "PyMolVisualization.h"
#include "compareRosettaModels.h"
#include "System.h"
#include "SysEnv.h"
#include "RosettaScoredPDBReader.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("compareRosettaModels");


int main(int argc, char *argv[]) {
  Options opt = setupOptions(argc,argv);
  SysEnv env;

  PyMolVisualization py;

  System sys1;
  sys1.readPdb(opt.pdb1);
  
  // Read score info
  RosettaScoredPDBReader rin1;
  rin1.open(opt.pdb1);
  rin1.read();
  rin1.close();  
  map<string,map<string,double> > rscores1 = rin1.getResidueScores();
  
  System sys2;
  sys2.readPdb(opt.pdb2);

  // Read score info
  RosettaScoredPDBReader rin2;
  rin2.open(opt.pdb2);
  rin2.read();
  rin2.close();  
  map<string,map<string,double> > rscores2 = rin2.getResidueScores();

  // Select interesting positions
  vector<string> interestingPositions;
  if (opt.resfile != ""){
  // Read resfile
    parseResfile(interestingPositions, opt.resfile);
  } else {

    // Add all positions from PDB # 1
    for (uint i = 0; i < sys1.positionSize();i++){
      interestingPositions.push_back(sys1.getPosition(i).getPositionId());
    }

  }


  stringstream pymolSel_rep;
  pymolSel_rep << "select repShiftedMutants, resi ";

  //cout << "Sys1: "<<sys1.toString()<<endl<<"Sys2: "<<sys2.toString()<<endl;
  ofstream fout;
  fout.open(MslTools::stringf("%s.scores", opt.out.c_str()).c_str());
  fout << MslTools::stringf("# File1   = %s\n",opt.pdb1.c_str());
  fout << MslTools::stringf("# File2   = %s\n",opt.pdb2.c_str());
  fout << MslTools::stringf("# Resfile = %s\n",opt.resfile.c_str());
  fout << MslTools::stringf("# Dir     = %s\n",env.getEnv("PWD").c_str());
  fout << MslTools::stringf("%10s %8s %8s %8s\n","Position","fa_rep1","fa_rep2","rep1-rep2");
  bool firstPyrep = true;
  for (uint i = 0; i < interestingPositions.size();i++){

    // Score fa_rep differences..
    double fa_rep1 = rscores1[interestingPositions[i]]["fa_rep"];
    double fa_rep2 = rscores2[interestingPositions[i]]["fa_rep"];

    if (!sys1.positionExists(interestingPositions[i])) {
      //cout << "Position: "<<interestingPositions[i]<<" does not exist in file: "<<opt.pdb1<<endl;
      continue;
    }
    if (!sys2.positionExists(interestingPositions[i])) {
      //cout << "Position: "<<interestingPositions[i]<<" does not exist in file: "<<opt.pdb2<<endl;
      continue;
    }

    if (fa_rep1-fa_rep2 < 3.0){

      if (!firstPyrep){
	  pymolSel_rep << "+";
      }
      firstPyrep = false;
      pymolSel_rep<<sys1.getPosition(interestingPositions[i]).getResidueNumber();
    }
    fout << MslTools::stringf("%25s %25s %10s %8.3f %8.3f %8.3f\n",MslTools::getFileName(opt.pdb1).c_str(),MslTools::getFileName(opt.pdb2).c_str(),interestingPositions[i].c_str(),fa_rep1,fa_rep2,fa_rep1-fa_rep2);
  }
  fout.close();

  map<string,int> positionChangeCount;
  fout.open(MslTools::stringf("%s.dist", opt.out.c_str()).c_str());
  fout << MslTools::stringf("# File1   = %s\n",opt.pdb1.c_str());
  fout << MslTools::stringf("# File2   = %s\n",opt.pdb2.c_str());
  fout << MslTools::stringf("# Resfile = %s\n",opt.resfile.c_str());
  fout << MslTools::stringf("# Dir     = %s\n",env.getEnv("PWD").c_str());
  fout << MslTools::stringf("%10s %10s %8s %8s %8s\n","Position1","Position2","Dist1","Dist2","D1-D2");
  for (uint i = 0; i < interestingPositions.size();i++){
    if (!sys1.positionExists(interestingPositions[i])) {
      //cout << "Position: "<<interestingPositions[i]<<" does not exist in file: "<<opt.pdb1<<endl;
      continue;
    }
    if (!sys2.positionExists(interestingPositions[i])) {
      //cout << "Position: "<<interestingPositions[i]<<" does not exist in file: "<<opt.pdb2<<endl;
      continue;
    }
    for (uint j = 0; j < interestingPositions.size();j++){

      if (!sys1.positionExists(interestingPositions[j])) {
	//cout << "Position: "<<interestingPositions[j]<<" does not exist in file: "<<opt.pdb1<<endl;
	continue;
      }

      if (!sys2.positionExists(interestingPositions[j])) {
	//cout << "Position: "<<interestingPositions[j]<<" does not exist in file: "<<opt.pdb2<<endl;
	continue;
      }

      double dist1 = sys1.getPosition(interestingPositions[i]).getAtom("CA").distance(sys1.getPosition(interestingPositions[j]).getAtom("CA"));
      double dist2 = sys2.getPosition(interestingPositions[i]).getAtom("CA").distance(sys2.getPosition(interestingPositions[j]).getAtom("CA"));
      
      if ( (dist1 < 12 || dist2 < 12) && 
           (abs(dist1-dist2) > 5 )  ){
	  fout << MslTools::stringf("%10s %10s %8.3f %8.3f %8.3f\n",interestingPositions[i].c_str(),interestingPositions[j].c_str(),dist1,dist2, dist1-dist2);
	  positionChangeCount[interestingPositions[i]]++;
      }
    }
	
  }
  fout.close();

  map<string,int>::iterator it;
  map<int,vector<string > > newSortedMap;
  for (it = positionChangeCount.begin();it != positionChangeCount.end();it++){
    newSortedMap[it->second].push_back(it->first);
  }
  stringstream pymolSel_env;
  pymolSel_env << "select envShiftedMutants, resi ";
  bool firstPysel = true;
  fout.open(MslTools::stringf("%s.contacts", opt.out.c_str()).c_str());
  fout << MslTools::stringf("# File1   = %s\n",opt.pdb1.c_str());
  fout << MslTools::stringf("# File2   = %s\n",opt.pdb2.c_str());
  fout << MslTools::stringf("# Resfile = %s\n",opt.resfile.c_str());
  fout << MslTools::stringf("# Dir     = %s\n",env.getEnv("PWD").c_str());
  fout << MslTools::stringf("%15s %4s\n","Position","ContactChanges");

  map<int,vector<string> >::iterator it2;
  for (it2 = newSortedMap.begin();it2 != newSortedMap.end();it2++){
    for (uint i = 0; i < it2->second.size();i++){
      fout << MslTools::stringf("%15s %4d\n",it2->second[i].c_str(),it2->first);
      if (it2->first >= 3){
	string chain;
	int resnum;
	string icode;
	MslTools::parsePositionId(it2->second[i],chain,resnum,icode);
	if (!firstPysel){
	  pymolSel_env << "+";
	}
	firstPysel = false;
	pymolSel_env<<resnum;
      }
    }
  }
  fout.close();
  fout.open(MslTools::stringf("%s.env.py", opt.out.c_str()).c_str());
  fout << MslTools::stringf("# File1   = %s\n",opt.pdb1.c_str());
  fout << MslTools::stringf("# File2   = %s\n",opt.pdb2.c_str());
  fout << MslTools::stringf("# Resfile = %s\n",opt.resfile.c_str());
  fout << MslTools::stringf("# Dir     = %s\n",env.getEnv("PWD").c_str());
  fout << "cmd.do(\""<<pymolSel_env.str()<<"\")"<<endl;
  fout << "cmd.do(\""<<pymolSel_rep.str()<<"\")"<<endl;
  fout.close();
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
		cout << "compareRosettaModels --pdb1 PDB --pdb2 PDB --resfile RESFILE\n";
		exit(0);
	}
	opt.pdb1 = OP.getString("pdb1");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb1 not specified.\n";
		exit(1111);
	}
	opt.pdb2 = OP.getString("pdb2");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb2 not specified.\n";
		exit(1111);
	}
	opt.out = OP.getString("out");
	if (OP.fail()){
	  opt.out = "compare";
	}
	opt.resfile = OP.getString("resfile");
	if (OP.fail()){
	  opt.resfile = "";
	}
	
	return opt;
}


void parseResfile(vector<string> &_interestingPositions, string _filename){

  vector<string> resfileLines;
  MslTools::readTextFile(resfileLines,_filename);
  bool start_parsing = false;
  for (uint i = 0; i < resfileLines.size();i++){
    if (resfileLines.size() == 0) continue;

    vector<string> toks = MslTools::tokenize(resfileLines[i]);
    if (toks[0] == "start") {
      start_parsing =true;
      continue;
    }

    if (!start_parsing) continue;
    if (toks[0][0] == '#') continue;

    // Strip I-code
    string resnum = "";
    string icode  = "";
    if (isalpha(toks[0][toks[0].size()-1])){
      icode = toks[0].substr(toks[0].size()-1,1);
      resnum = toks[0].substr(0,toks[0].size()-1);
    } else {
      resnum = toks[0].substr(0,toks[0].size());
    }
    //cout << "RESFILE: "<<resfileLines[i]<<endl;
    //cout << "Resfile PosId: "<<MslTools::getPositionId(toks[1],MslTools::toInt(resnum),icode)<<" , " <<toks[1]<<" and "<<resnum<<endl;
    _interestingPositions.push_back(MslTools::getPositionId(toks[1],MslTools::toInt(resnum),icode));
  }
}
