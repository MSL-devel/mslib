#include "System.h"
#include "PDBTopology.h"
#include "AtomContainer.h"
#include "testData.h"
#include "OptionParser.h"
#include "PDBWriter.h"

using namespace MSL;
using namespace std;

#include "SysEnv.h"
static SysEnv SYSENV;

// Define charmm topology and rotamer library files.  CharmmTop is not used, but required to be specified at this point.
string charmmTop = SYSENV.getEnv("CHARMMTOP");
string rotlib    = SYSENV.getEnv("ROTLIB");

int main(int argc, char *argv[]) {

  OptionParser op;
  op.readArgv(argc,argv);
  
  writePdbFile();

  // Build a tmp system  
  System sys;
  sys.readPdb("/tmp/testPdb.pdb");

  // Get one of the residues
  Residue &res = sys.getIdentity("A,2,LEU");

  // Get backbone atoms  
  AtomPointerVector &backboneAtoms = sys.getIdentity("A,2").getAtomPointers();

  PDBTopology pdbTop;
  pdbTop.readCharmmTopology(charmmTop);
  pdbTop.readRotamerLibrary(rotlib);
  pdbTop.setAddAtomsFromRotLib(true);

  // For each amino acid
  vector<string> aas;
  //aas.push_back("ALA");
  //aas.push_back("GLY");

  aas.push_back("CYS");
  aas.push_back("ASP");
  aas.push_back("GLU");
  aas.push_back("PHE");
  aas.push_back("HSD");
  aas.push_back("HSE");
  aas.push_back("ILE");
  aas.push_back("LYS");
  aas.push_back("LEU");
  aas.push_back("MET");
  aas.push_back("ASN");
  aas.push_back("GLN");
  aas.push_back("ARG");
  aas.push_back("SER");
  aas.push_back("THR");
  aas.push_back("VAL");
  aas.push_back("TRP");
  aas.push_back("TYR");


  for (uint i = 0; i < aas.size();i++){
    string aa = aas[i];
    cout << "Building "<<aa<<endl;

    AtomContainer newResidue = pdbTop.getResidue("A,2,"+aa,backboneAtoms,3);

    PDBWriter pout;
    pout.open("/tmp/newRes-"+aa+".pdb");
    for (uint j = 0; j < newResidue[0].getNumberOfAltConformations();j++){
      for (uint a = 0; a< newResidue.size();a++){
	cout << "Atom : "<<newResidue[a].toString()<<endl;
	newResidue[a].setActiveConformation(j);
      }

      pout.write(newResidue.getAtomPointers(),false,false,true);
    }
    pout.close();

  }


  // Build one inverse-rotamer 


  // Get some side-chain atoms
  AtomPointerVector sideChainAtoms;
  sideChainAtoms.push_back(&res("CD1"));
  sideChainAtoms.push_back(&res("CD2"));
  sideChainAtoms.push_back(&res("CG"));
  sideChainAtoms.push_back(&res("CB"));

  AtomContainer newResidue = pdbTop.getResidue("A,2,LEU",sideChainAtoms,1);
  System newRes;
  newRes.addAtoms(newResidue.getAtomPointers());
  newRes.writePdb("/tmp/newRes-LEU-inverse.pdb");
}
