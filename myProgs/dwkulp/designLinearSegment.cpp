#include <iostream>
#include <cstdlib>
#include <queue>
#include "signal.h"
#include "System.h"
#include "Frame.h"
#include "Timer.h"
#include "RegEx.h"
#include "PDBFragments.h"
#include "PDBWriter.h"
#include "SasaCalculator.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "PSSMCreator.h"
#include "designLinearSegment.h"



using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("designLinearSegment");


// Utility functions
void sighandler(int sig);


struct aaDat {
  double freq;
  string aa;

  aaDat(){
    freq = 0.0;
    aa = "";
  }
  const aaDat& operator=(const aaDat &_rhs){
    if (this == &_rhs) return *this;

    freq  = _rhs.freq;
    aa    = _rhs.aa;
    
    return (*this);
  }

};
bool operator<( const aaDat& p1, const aaDat& p2 ) {  
      return p1.freq > p2.freq;
}




priority_queue<aaDat> orderMap(map<string,double> &_aMap);
void printFreqs(vector<map<string,double> > &_freqs);



int main(int argc, char *argv[]) {
  
  // Register signal handlers
  signal(SIGABRT, &sighandler);
  signal(SIGTERM, &sighandler);
  signal(SIGINT, &sighandler);

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // MslOut can suppress output, to make example output clean
  MSLOUT.turnOn("designLinearSegments");

  // Read input pdb
  System sys;
  sys.readPdb(opt.pdb);
  
  PDBFragments fragDB(opt.fragdb);
  fragDB.setPdbDir(opt.pdbDir);
  fragDB.loadFragmentDatabase();
  fragDB.setIncludeFullFile(true);

  // Do local sampling inside PDBFragment object
  cout << "Search"<<endl;
  int numMatchingFrags = 0;
  
  switch (opt.searchTypeInt){
    case PDBFragments::linear:
      numMatchingFrags = fragDB.searchForMatchingFragmentsLinear(sys,opt.stems[0],opt.stems[1],opt.regex, opt.rmsdTol);
      break;
    case PDBFragments::stemOnly:
      numMatchingFrags = fragDB.searchForMatchingFragmentsStems(sys,opt.stems,opt.num_residues,opt.regex,opt.rmsdTol);
      break;
    case PDBFragments::discreteSpots:
      numMatchingFrags = fragDB.searchForMatchingFragmentsSpots(sys,opt.stems,opt.num_residues,opt.rmsdTol);
      break;
  }
  cout << "Done search found "<<numMatchingFrags<<endl;

  // For now just dump all matches
  vector<AtomContainer *> results = fragDB.getAtomContainers();
  for (uint i = 0; i < results.size();i++){

    cout << "Writing: "<<MslTools::stringf("%s_%06d.pdb",opt.outpdb.c_str(),i)<<endl;
      System newSys;
      newSys.addAtoms(results[i]->getAtomPointers());
      newSys.writePdb(MslTools::stringf("%s_%06d.pdb",opt.outpdb.c_str(),i));

  }


  PSSMCreator pssm;
  pssm.setSequences(fragDB.getMatchedSequences());
  pssm.create(PSSMCreator::freq);
  vector<map<string,double> > freqs = pssm.getFrequencies();
  printFreqs(freqs);
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
    cout << "searchForFusions --ref PDB --list pdbs_to_search\n";
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
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }

  opt.fragdb = OP.getString("fragdb");
  if (OP.fail()){
    cerr << "ERROR 1111 fragdb not specified.\n";
    exit(1111);
  }

  opt.bbq = OP.getString("bbq");
  if (OP.fail()){
    opt.bbq = "";
  }
  
  opt.pdbDir = OP.getString("pdbdir");
  if (OP.fail()){
    opt.pdbDir = "";
  }  
  opt.stems = OP.getMultiString("stems");
  if (OP.fail()){
    cerr << "ERROR 1111 stems not specified.\n";
    exit(1111);
  }

  opt.num_residues = OP.getInt("numResidues");
  if (OP.fail()){
    opt.num_residues = -1;
  }
  opt.regex = OP.getString("regex");
  if (OP.fail()){
    opt.regex = "";
  }

  opt.outpdb = OP.getString("outpdb");
  if (OP.fail()){
    if (opt.num_residues != -1){
      opt.outpdb = MslTools::stringf("%s_%02d_segdesign",MslTools::getFileName(opt.pdb).c_str(),opt.num_residues);
    }else {
      opt.outpdb = MslTools::stringf("%s_segdesign",MslTools::getFileName(opt.pdb).c_str());
    }
  }

  opt.searchType = OP.getString("searchType");
  if (OP.fail()){
    opt.searchTypeInt = PDBFragments::stemOnly;
  } else {
    map<string,int>::iterator it;
    it = opt.searchTypeMap.find(opt.searchType);
    if (it == opt.searchTypeMap.end()){
      cerr << "ERROR 1111 searchType = '"<<opt.searchType<<"' was not found. Try one of:"<<endl;
      for (it = opt.searchTypeMap.begin(); it != opt.searchTypeMap.end();it++){
	cerr << "\t"<<it->first<<endl;
      }
      exit(1111);
    }

    opt.searchTypeInt = it->second;
  }

  opt.rmsdTol = OP.getDouble("rmsdTol");
  if (OP.fail()){
    opt.rmsdTol = 1.0;
  }
  return opt;
}



void sighandler(int sig){

     cerr << "ERROR  ****  SIGNAL "<<sig<<" caught, try to output best fusions so far:"<<endl;

}


void printFreqs(vector<map<string,double> > &_freqs){



    // For each position
    int largestVariation = 0;
    for (uint i = 0;i < _freqs.size();i++){
      if (_freqs[i].size() > largestVariation){
	largestVariation = _freqs[i].size();
      }
    }

    MSLOUT.stream() << "Make output largest variation: "<<largestVariation<<endl;
    int resNum = 1;

    vector<string> rows;
    rows.resize(21,"");
    for (uint i = 0;i < _freqs.size();i++){


      priority_queue<aaDat> orderedMap = orderMap(_freqs[i]);
      int blanks = largestVariation - orderedMap.size();
      for (uint v = 0; v < blanks;v++){
	if (v <= blanks) {
	  rows[v] += MslTools::stringf("%1s %6.2f\t", "-", 00.00);
	  continue;
	}
      }

      int v = blanks;
      while (orderedMap.size() > 0){
	  aaDat data = orderedMap.top();
	  orderedMap.pop();
	  rows[v++] += MslTools::stringf("%1s %6.2f\t", data.aa.c_str(), data.freq*100);
	  //cout << "ROW IS: "<<MslTools::stringf("%1s %8.2f ", data.aa.c_str(), data.freq*100)<<endl;
      }
      
       
      // Add residue number at very end
      //rows[20] += MslTools::stringf("=%c-%-4d=\t",seq[i],resNum++);
      
    } // FOR freqs.size()

    MSLOUT.stream() << "Make output2"<<endl;
    for (uint r = 0; r < rows.size();r++){
      if (rows[r] == "") continue;
      if (r == 20){
	fprintf(stdout,"\n");
      }
      fprintf(stdout, "%s\n",rows[r].c_str());
    }
}

priority_queue<aaDat> orderMap(map<string,double> &_aMap){

  priority_queue<aaDat> results;
  map<string,double>::iterator it;
  for (it = _aMap.begin();it != _aMap.end();it++){
    aaDat a;
    a.freq = it->second;
    a.aa   = it->first;
    results.push(a);
  }

  return results;
}
