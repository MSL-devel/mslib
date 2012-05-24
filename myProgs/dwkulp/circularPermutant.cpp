#include <iostream>
#include <cstdlib>
#include <queue>
#include "signal.h"
#include "System.h"
#include "Frame.h"
#include "Timer.h"
#include "RegEx.h"
#include "Transforms.h"
#include "PDBWriter.h"
#include "SasaCalculator.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include "MslOut.h"


#include "circularPermutant.h"



using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("circularPermutant");


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // MslOut can suppress output, to make example output clean
  MSLOUT.turnOn("circularPermutant");

  // Read in the reference pdb file
  System sys;
  sys.readPdb(opt.pdb);

  System circularPerm;
  for (uint c = 0; c < sys.chainSize();c++){
    Chain &ch = sys.getChain(c);

    int resNum = 1;

    // Get index of newStartRes
    if (!ch.positionExists(MslTools::stringf("%s,%s",ch.getChainId().c_str(), opt.newStartRes.c_str())) ||
        !ch.positionExists(MslTools::stringf("%s,%s",ch.getChainId().c_str(), opt.newEndRes.c_str()))){
      cout << "PDB: "<<MslTools::getFileName(opt.pdb)<<" Chain "<<ch.getChainId()<<" does not have start or end residues : "<<opt.newStartRes<<" , "<<opt.newEndRes<<endl;
      continue;
    }
    int startIndex = ch.getPositionIndex(&ch.getPosition(MslTools::stringf("%s,%s",ch.getChainId().c_str(), opt.newStartRes.c_str())));
    for (uint i = startIndex; i < ch.positionSize();i++){
      Position &pos = ch.getPosition(i);
      AtomContainer ac;      
      ac.addAtoms(pos.getAtomPointers());

      for (uint a = 0 ; a < ac.size();a++){
	ac(a).setResidueNumber(resNum);
      }
      circularPerm.addAtoms(ac.getAtomPointers());
      resNum++;

    }

    // Now add from original residue number 1 to endIndex.
    int endIndex   =  ch.getPositionIndex(&ch.getPosition(MslTools::stringf("%s,%s",ch.getChainId().c_str(), opt.newEndRes.c_str())));
    for (uint i = 0; i <= endIndex; i++){
      Position &pos = ch.getPosition(i);
      AtomContainer ac;      
      ac.addAtoms(pos.getAtomPointers());

      for (uint a = 0 ; a < ac.size();a++){
	ac(a).setResidueNumber(resNum);
      }
      circularPerm.addAtoms(ac.getAtomPointers());      
      resNum++;
    }
  }

  circularPerm.writePdb(opt.outPdb);


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
    cout << "circularPermutant --pdb PDB --newStartRes POS_ID --newEndRes POS_ID\n";
    exit(0);
  }


  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }
  opt.newStartRes = OP.getString("newStartRes");
  if (OP.fail()){
    cerr << "ERROR 1111 newStartRes not specified.\n";
    exit(1111);
  }
  opt.newEndRes = OP.getString("newEndRes");
  if (OP.fail()){
    cerr << "ERROR 1111 newEndRes not specified.\n";
    exit(1111);
  }
  opt.outPdb = OP.getString("outPdb");
  if (OP.fail()){
    opt.outPdb = MslTools::stringf("%s_circ.pdb",MslTools::getFileName(opt.pdb).c_str());
    cout << "Output pdb: "<<opt.outPdb<<endl;
  }
  return opt;
}

