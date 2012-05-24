/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include "MslTools.h"
#include "OptionParser.h"
#include "RegEx.h"
#include "PhiPsiStatistics.h"
#include "PhiPsiReader.h"
#include "SysEnv.h"
#include "PrositeReader.h"
#include "release.h"

#include "System.h"
#include "Residue.h"

using namespace std;
using namespace MSL;
#include <fstream>

#include "designCheck.h"

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	System sys;
	sys.readPdb(opt.pdb);

	// Single body checks
	//  - glycine-rich position
	//  - proline-rich position
	//  - unpaired cys
	//  - rotamer check
	reportPositionMetrics(sys,opt);

	// Sequence Motif checks
	//  - minimially populated seq motifs
	//  - N-linked glycosylation sites
	//  - cleavage sites
	reportSequenceMotifMetrics(sys,opt);
	

	// Structure motif checks
	//  - backbone structure based sequence biases
	//  - Salt bridges
	//  - Disulfides
	//  - Loop lengths / sequence composition
	//  - Capping (ends of helices, turns, etc..)
	reportStructureMotifMetrics(sys,opt);


	ofstream pout;
	pout.open(MslTools::stringf("%s.report.py",MslTools::getFileName(opt.pdb).c_str()).c_str());
	pout << opt.pymol.toString();
	pout.close();


}
void reportStructureMotifMetrics(System &_sys, Options &_opt){ 
}

void reportSequenceMotifMetrics(System &_sys, Options &_opt){ 

  System ref;
  if (_opt.ref != ""){
    ref.readPdb(_opt.ref);
  }

  // Prosite Reader
  PrositeReader prosite;
  prosite.open(_opt.prosite);
  prosite.read();
  prosite.close();

  map<string,map<string,string> > prosite_data = prosite.getPrositeData();
  map<string,map<string,string> >::iterator proIt;
  map<string,string>::iterator proIt2;

  // Regular expression object
  RegEx re;
	//vector<pair<int,int> > matches = re.getResidueRanges(ch,"N.[ST]");	

  // For each chain
  for (uint c = 0; c < _sys.chainSize();c++){
    Chain &ch = _sys.getChain(c);

    int prositePatternFoundIndex = 1;
    for (proIt = prosite_data.begin();proIt != prosite_data.end();proIt++){
      for (proIt2 = proIt->second.begin();proIt2 != proIt->second.end();proIt2++){
	
	vector<pair<int,int> > matches;
        matches = re.getResidueRanges(ch,proIt2->second);

        for (uint m = 0; m < matches.size();m++){

	  string seq = "";
	  string selection = MslTools::stringf("chain %1s and resi %5d-%-5d",ch.getChainId().c_str(),ch.getResidue(matches[m].first).getResidueNumber(),ch.getResidue(matches[m].second).getResidueNumber());
	  string refSeq = "";

	  for (uint r = matches[m].first; r <= matches[m].second;r++){
	    seq += MslTools::getOneLetterCode(ch.getResidue(r).getResidueName());

	    if (_opt.ref != ""){
	      if (ref.positionExists(ch.getResidue(r).getPositionId())){
		refSeq += MslTools::getOneLetterCode(ref.getResidue(ch.getResidue(r).getPositionId()).getResidueName());
	      }
	    }
	  }
	  
	  // Skip if reference structure has same sequence motif
	  if (_opt.ref != "" && seq == refSeq){
	    string refStr = "**** IN REF STRUCTURE ****";
	    fprintf(stdout,"%-40s %20s %3s %s\n", proIt2->first.substr(0,40).c_str(), selection.c_str(),seq.c_str(),refStr.c_str());
	    continue;
	  }

	  fprintf(stdout,"%-40s %20s %3s\n", proIt2->first.substr(0,40).c_str(), selection.c_str(),seq.c_str());

	  string selStr = MslTools::stringf("%s and chain %1s and resi %d-%-d", MslTools::getFileName(_opt.pdb).c_str(),ch.getResidue(matches[m].first).getChainId().c_str(), ch.getResidue(matches[m].first).getResidueNumber(),ch.getResidue(matches[m].second).getResidueNumber());
	  if (ch.getResidue(matches[m].first).getChainId() == "" || ch.getResidue(matches[m].first).getChainId() == " "){
	    selStr = MslTools::stringf("%s and resi %d-%-d",MslTools::getFileName(_opt.pdb).c_str(),ch.getResidue(matches[m].first).getResidueNumber(),ch.getResidue(matches[m].second).getResidueNumber());
	  }

	  string selName = MslTools::stringf("%s-%s-%-d",MslTools::getFileName(_opt.pdb).c_str(),proIt->first.c_str(),prositePatternFoundIndex++);
	  _opt.pymol.createSelection(selName,selStr);

	} // FOR EACH MATCHES

      }// FOR EACH PROSITE2      
    } // FOR EACH PROSITE

  } // FOR EACH CHAIN

}

void reportPositionMetrics(System &_sys, Options &_opt){ 

       System ref;
       if (_opt.ref != ""){
	 ref.readPdb(_opt.ref);
       }

	SysEnv env;
	PhiPsiReader ppr(env.getEnv("MSL_PHIPSI_TABLE"));
	ppr.open();
	ppr.read();
	ppr.close();
	PhiPsiStatistics pps = ppr.getPhiPsiStatistics();
	

	int unpaired_cys_index = 1;
	for (uint c = 0; c < _sys.chainSize();c++){
	  for (uint r = 0; r < _sys.getChain(c).positionSize();r++){
	    Residue &res = _sys.getChain(c).getResidue(r);

	    // Remove non-amino acids... bad way to do this, but should work.
	    string oneLetter = MslTools::getOneLetterCode(res.getResidueName());
	    if (oneLetter == "X") continue;



	    if (r == 0 || r == _sys.positionSize()-1 ||
		MslTools::getOneLetterCode(_sys.getChain(c).getResidue(r+1).getResidueName()) == "X"  || 
		MslTools::getOneLetterCode(_sys.getChain(c).getResidue(r-1).getResidueName()) == "X") continue;

	    
	    Residue & nm1 = _sys.getChain(c).getResidue(r-1);
	    Residue & np1 = _sys.getChain(c).getResidue(r+1);

	    // Get Phi-Psi
	    double phi = PhiPsiStatistics::getPhi(nm1,res);
	    double psi = PhiPsiStatistics::getPsi(res,np1);

	    // Check glycine/proline populations
	    int resCounts  = pps.getCounts(res.getResidueName(),phi,psi);
	    int glyCounts  = pps.getCounts("GLY",phi,psi);
	    double glyFreq = pps.getFreqInBin("GLY",phi,psi);
	    int proCounts  = pps.getCounts("PRO",phi,psi);
	    double proFreq = pps.getFreqInBin("PRO",phi,psi);
	    
	    string refStr = "";
	    if (_opt.ref != ""){
	      if (ref.positionExists(res.getPositionId())){
		if (ref.getPosition(res.getPositionId()).getResidueName() == res.getResidueName()){
		  refStr = "**** IN REF STRUCTURE ****";
		}
	      }

	    }


	    // if glyFreq > 0.25 
	    if (glyFreq > 0.25){
	      fprintf(stdout, "%-40s %10s %8d %8d %8.3f %s\n","High freq. Gly",res.getIdentityId().c_str(),resCounts,glyCounts,glyFreq,refStr.c_str());
	    }
	    // if proFreq > 0.25
	    if (proFreq > 0.25){
	      fprintf(stdout, "%-40s %10s %8d %8d %8.3f %s\n","High freq. Pro",res.getIdentityId().c_str(),resCounts,proCounts,proFreq,refStr.c_str());
	    }


	    // if CYS, check for pairing
	    if (res.getResidueName() == "CYS"){

	      if (!res.atomExists("SG")){
		cerr << "ERROR CYS residue "<<res.getPositionId()<<" does not have a SG atom"<<endl;
	      }
	      bool unpairedCys = true;
	      double minDist = MslTools::doubleMax;
	      string minId = "";
	      for (uint p = 0; p < _sys.positionSize();p++){
	        Residue &res2 = _sys.getResidue(p);

		if (res2.getPositionId() == res.getPositionId()) continue;
		if (res2.getResidueName() != "CYS") continue;
		if (!res2.atomExists("SG")) continue;

		double dist = res2.getAtom("SG").distance(res.getAtom("SG"));
		if (dist < minDist){
		  minDist = dist;
		  minId   = res2.getPositionId();
		}
		if (dist < 2.2){
		  unpairedCys = false;
		  break;
		}


	      } // FOR system positions p


	      if (unpairedCys){
		fprintf(stdout, "%-40s %10s %8.3f %10s %s\n","Unpaired Cys",res.getPositionId().c_str(),minDist,minId.c_str(),refStr.c_str());
		if (refStr != "") continue;

		string selStr = MslTools::stringf("%s and chain %1s and resi %d",MslTools::getFileName(_opt.pdb).c_str(),res.getChainId().c_str(), res.getResidueNumber());
		if (res.getChainId() == "" || res.getChainId() == " "){
		  selStr = MslTools::stringf("%s and resi %d",MslTools::getFileName(_opt.pdb).c_str(),res.getResidueNumber());
		}
		string selName = MslTools::stringf("%s-%s-%-d",MslTools::getFileName(_opt.pdb).c_str(),"UNPAIRED_CYS",unpaired_cys_index++);
		_opt.pymol.createSelection(selName,selStr);
	      }


	    } // IF res == CYS
				    
	  } // Positions on chain
	} // Chains in system
}


Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option


	if (OP.countOptions() == 0){
		cout << "Usage: printSequence " << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.prosite = OP.getString("prosite");
	if (OP.fail()){
	  cerr << "ERROR 1111 no prosite database specified."<<endl;
	  exit(1111);
	}
	opt.ref = OP.getString("ref");

	return opt;
}
