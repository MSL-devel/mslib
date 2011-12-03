#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "AtomSelection.h"
#include "PhiPsiStatistics.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "compareStructures.h"
#include "PDBTopology.h"
#include "VectorPair.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("compareStructures");

void alignAndCheck(System &_complex, string _selStrCom, System &_ref, string _selStrRef, Options &_opt, map<string, AtomPointerVector> &_domains, int _resNum1, int _resNum2, bool _refine);

int main(int argc, char *argv[]) {

	// Read cmdline options
	Options opt = setupOptions(argc,argv);
	System ref;
	ref.readPdb(opt.ref);
	AtomSelection selR(ref.getAtomPointers());
	map<string,AtomPointerVector> domains;
	for (uint d = 0; d < opt.domains.size();d++){

	  AtomPointerVector ats = selR.select(opt.domains[d][1]);  
	  domains[opt.domains[d][0]] = ats;
	}


	System complex;
	complex.readPdb(opt.complex);
	complex.saveCoor("orig");

	Chain &alignChain = complex.getChain(opt.chainToAlign);

	// For each of the other chains  in the complex:
	for (uint c = 0; c < opt.chainsToCheck.size();c++){
	  Chain &otherChain = complex.getChain(opt.chainsToCheck[c]);

	  if (otherChain.getChainId() == opt.chainToAlign) continue;

	  // Get Inteface Residues
	  pair<map<int,bool>,map<int,bool> > positionIndices = complex.findProteinInterfacePositions(alignChain.getChainId(),otherChain.getChainId());

	  map<int,bool>::iterator it;
	  map<int,bool>::iterator it2;

	  for (it=positionIndices.first.begin();it != positionIndices.first.end();it++){


	    for (it2=it;it2 != positionIndices.first.end();it2++){
	      if (it2 == it) continue;

	      int resNum1 = complex.getPosition(it->first).getResidueNumber();
	      int resNum2 = complex.getPosition(it2->first).getResidueNumber();

	      string selStrRef = MslTools::stringf("name CA and chain %s and resi %d-%d+%d-%d", opt.chainInRef.c_str(),resNum1-3,resNum1+3,resNum2-3,resNum2+3);
	      string selStrCom = MslTools::stringf("name CA and chain %s and resi %d-%d+%d-%d", opt.chainToAlign.c_str(),resNum1-3,resNum1+3,resNum2-3,resNum2+3);

	      alignAndCheck(complex, selStrCom, ref, selStrRef,opt, domains, resNum1, resNum2, true);


	      complex.applySavedCoor("orig");
	      

	    } // END positionIndicies.first it2



	  } // END positionIndicies.first it


	  	  
	}// END Complex chains
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
		cout << "compareStructures --ref PDB --complex PDB --chainInComplexToAlign CHAIN \n";
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


	opt.complex = OP.getString("complex");
	if (OP.fail()){
		cerr << "ERROR 1111 complex not specified.\n";
		exit(1111);
	}

	opt.ref = OP.getString("ref");
	if (OP.fail()){
		cerr << "ERROR 1111 ref not specified.\n";
		exit(1111);
	}

	opt.chainToAlign = OP.getString("chainInComplexToAlign");
	if (OP.fail()){
		cerr << "ERROR 1111 chainToAlign not specified.\n";
		exit(1111);
	}

	opt.chainsToCheck = OP.getStringVector("chainInComplexToCheck");
	if (OP.fail()){
		cerr << "ERROR 1111 chainToCheck not specified.\n";
		exit(1111);
	}

	opt.chainInRef = OP.getString("chainInRef");
	if (OP.fail()){
		cerr << "ERROR 1111 chainInRef not specified.\n";
		exit(1111);
	}

	vector<string> d = OP.getMultiString("domain");
	if (OP.fail()){
	  cerr << "WARNING no domains specified"<<endl;
	} else {
	  for (uint indexD = 0; indexD < d.size();indexD++){
	    cout << "Parsing: "<<d[indexD]<<endl;
	    vector<string> toks = MslTools::tokenize(d[indexD],";");
	    if (toks.size() % 2 != 0){
	      cerr << "ERROR 1111 domain must be paired ID;SELECTION\n";
	      exit(1111);
	    }
	    for (uint i = 0; i < toks.size();i+=2){
	      vector<string> id;
	      id.push_back(toks[i]);
	      id.push_back(toks[i+1]);
	      opt.domains.push_back(id);
	    }
	  }
	}
	

	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}
	return opt;
}


void alignAndCheck(System &_complex, string _selStrCom, System &_ref, string _selStrRef, Options &_opt, map<string, AtomPointerVector> &_domains, int _resNum1, int _resNum2, bool _refine){

  AtomSelection selR(_ref.getAtomPointers());
  AtomPointerVector refSel = selR.select(_selStrRef);
  AtomSelection selC(_complex.getChain(_opt.chainToAlign).getAtomPointers());
  AtomPointerVector comSel = selC.select(_selStrCom);


  if (refSel.size()  != comSel.size()){
    cerr << "ERROR refSel != comSel"<<endl;
    cerr <<" refSel("<<refSel.size()<<") = "<<_selStrRef<<endl;
    cerr <<" comSel("<<comSel.size()<<") = "<<_selStrCom<<endl;
    return;
  }

	      
  // Align it
  Transforms t;
  if (!t.rmsdAlignment(comSel,refSel,_complex.getAtomPointers())){
    cerr << "ERROR alignment failed."<<endl;
    return;
  }

  double rmsd = comSel.rmsd(refSel);
  //cout << "RMSD: "<<rmsd<<endl;
  // Get number of clashes with this alignment for all the "other" chains.

  for (uint c3 = 0; c3 < _opt.chainsToCheck.size();c3++){
    if (_complex.getChain(_opt.chainsToCheck[c3]).getChainId() == _opt.chainToAlign) continue;
    if (_complex.getChain(_opt.chainsToCheck[c3]).getChainId() == "") continue;
		
    map<string,AtomPointerVector>::iterator dIt;

    bool goodModel = true;
    for (dIt = _domains.begin();dIt != _domains.end();dIt++){

      int clashes = 0;
      for (uint dAt = 0; dAt != dIt->second.size();dAt++){
	for (uint cAt = 0; cAt != _complex.getChain(_opt.chainsToCheck[c3]).atomSize();cAt++){
	  if (_complex.getChain(_opt.chainsToCheck[c3]).getAtom(cAt).getName() != "CA") continue;
		      
	  double distSq = _complex.getChain(_opt.chainsToCheck[c3]).getAtom(cAt).distance2(*dIt->second[dAt]);
	  //cout << "Atoms: "<<complex.getChain(c3).getAtom(cAt).getAtomId()<<" "<<dIt->second[dAt]->getAtomId()<<endl;
	  if (distSq < 6.25){
	    clashes++;
	  }

	} // END COMPLEX ATOMS
      } // END DOMAINS ATOMS

      fprintf(stdout, "%1s, %6s, %6d, %8.3f, %50s\n", _complex.getChain(_opt.chainsToCheck[c3]).getChainId().c_str(), dIt->first.c_str(),clashes,rmsd,_selStrCom.c_str());
		  
      if (clashes > 7){
	goodModel = false;
      }
		  
    } // END DOMAINS


    if (goodModel){


      // Tweak alignment to search for better match	
      if (_refine){
	fprintf(stdout, "MATCH: %1s, %8.3f, %50s\n", _complex.getChain(_opt.chainsToCheck[c3]).getChainId().c_str(),rmsd,_selStrCom.c_str());
	for (uint r1s = -2; r1s <= 2; r1s++){
	  for (uint r1e = -2; r1e <= 2; r1e++){
	    for (uint r2s = -2; r2s <= 2; r2s++){
	      for (uint r2e = -2; r2e <= 2; r2e++){
			  
	      string selStrRef = MslTools::stringf("name CA and chain %s and resi %d-%d+%d-%d", _opt.chainInRef.c_str(),_resNum1+r1s,_resNum1+r1e,_resNum2+r2s,_resNum2+r2e);
	      string selStrCom = MslTools::stringf("name CA and chain %s and resi %d-%d+%d-%d", _opt.chainToAlign.c_str(),_resNum1+r1s,_resNum1+r1e,_resNum2+r2s,_resNum2+r2e);
			  
	      alignAndCheck(_complex, selStrCom, _ref, selStrRef, _opt, _domains,_resNum1, _resNum2, false);
	      
			  
	      } // END R2E
	    } // END R2S
	  } // END R1E
	} // END R1S

      } else { // IF _refine
	fprintf(stdout, "REFINE: %1s, %8.3f, %50s\n", _complex.getChain(_opt.chainsToCheck[c3]).getChainId().c_str(),rmsd,_selStrCom.c_str());
      }


    } // IF goodModel
		
    _complex.applySavedCoor("orig");
		
  } // END COMPLEX CHAINS
}
