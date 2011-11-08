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
#include "findCisPepBonds.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("findCisPepBonds");

int main(int argc, char *argv[]) {

	// Read cmdline options
	Options opt = setupOptions(argc,argv);

	cout << "READ LIST"<<endl;
	vector<string> pdbs;  
	ifstream fs;

	fs.open(opt.list.c_str());
	if (fs.fail()){
		cerr<<"Cannot open file "<<opt.list<<endl;
		exit(1);
	}

	while(true){
		string line;
		getline(fs, line);

		if(fs.fail()){
			//no more lines to read, quite the while.
			break;
		}

		if(line==""){
			continue;
		}
		pdbs.push_back(line);
	}

	fs.close();

	// AA type
	map<string, int> cisPepBondCounts;
	map<string, int> cisPepBondCountsTwoAA;
	map<int,int> cisChain;
	map<int,string > cisChainId;
	int totalResidues = 0;
	cout << "Process the files.."<<endl;
	for (uint i = 0; i < pdbs.size();i++){

		// Read in PDB with no Hydrogens..
		System sys;
		PDBReader rin(pdbs[i]);
		rin.open();
		rin.read(true);
		sys.addAtoms(rin.getAtomPointers());
		rin.close();
		

		for (uint c = 0; c < sys.chainSize();c++){
		  if (sys.getChain(c).positionSize() < 10) continue;

		  int previousAA_cisBond = 0;
		  for (uint p = 0; p < sys.getChain(c).positionSize()-3;p++){
		    
		    Position &pos1 = sys.getChain(c).getPosition(p);
		    Position &pos2 = sys.getChain(c).getPosition(p+1);

		    
		    double omega = PhiPsiStatistics::getOmega(pos1.getCurrentIdentity(),pos2.getCurrentIdentity());


		    // CIS Omega angle in degrees
		    if (omega < 20 && omega > -20){

		      MSLOUT.fprintf(stdout,"CIS OMEGA: %8.3f from %s, %8s\n",omega,MslTools::getFileName(pdbs[i]).c_str(),pos1.getCurrentIdentity().getIdentityId().c_str());
		      stringstream ss;
		      ss << pos1.getResidueName();
		      cisPepBondCounts[ss.str()]++;
		      ss<<"-"<<pos2.getResidueName();
		      cisPepBondCountsTwoAA[ss.str()]++;
		      
		      if (previousAA_cisBond > 0){
			
			
			cisChain[previousAA_cisBond]++;

			char tmp[80];
			sprintf(tmp," %s,%s :",MslTools::getFileName(pdbs[i]).c_str(),pos1.getCurrentIdentity().getIdentityId().c_str());

			stringstream tmpSS;
			tmpSS <<  cisChainId[previousAA_cisBond] << tmp;
			cisChainId[previousAA_cisBond] = tmpSS.str();

			
		      }
		      previousAA_cisBond++;

		    } else {
			  previousAA_cisBond = 0;
		    }
			

		    totalResidues++;
		  }//Positions
		}//Chains
	}//PDBs


	// Print out the results...
	MSLOUT.stream() << "PROCESSED "<<totalResidues<<" residues"<<endl;
	map<string,int>::iterator it;
	for (it = cisPepBondCounts.begin();it != cisPepBondCounts.end();it++){
	  MSLOUT.stream() << "CIS BOND1["<<it->first<<"]: "<<it->second<<endl;
	}

	for (it = cisPepBondCountsTwoAA.begin();it != cisPepBondCountsTwoAA.end();it++){
	  MSLOUT.stream() << "CIS BOND2["<<it->first<<"]: "<<it->second<<endl;
	}
	map<int,int>::iterator it2;
	for (it2 = cisChain.begin();it2 != cisChain.end();it2++){
	  MSLOUT.stream() << "CIS BOND3["<<it2->first<<"]: "<<it2->second<< " from "<<cisChainId[it2->first]<<endl;
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
		cout << "findCisPepBonds --list LIST \n";
		exit(0);
	}
	opt.list = OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}
	return opt;
}
