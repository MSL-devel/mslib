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
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "BaselineEnergyBuilder.h"
#include "MslOut.h"
#include "CharmmEnergyCalculator.h"
#include "System.h"
#include "Residue.h"

using namespace std;
using namespace MSL;
#include <fstream>

#include "designCheck.h"

static SysEnv SYSENV;
static MslOut MSLOUT("designCheck");

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


	// Energy checks
	//   - Residues with worst energies
	reportEnergyMetrics(sys,opt);


	ofstream pout;
	pout.open(MslTools::stringf("%s.report.py",MslTools::getFileName(opt.pdb).c_str()).c_str());
	pout << opt.pymol.toString();
	pout.close();


}

void reportEnergyMetrics(System &_sys, Options &_opt){ 

  MSL::PolymerSequence pseq;
  pseq.setPDBNamesFlag(true); //Converts HIS to HSD
  pseq.setSequence(_sys);

  System ene_sys;

  CharmmSystemBuilder CSB(ene_sys, _opt.topfile, _opt.parfile);
  CSB.setDielectricConstant(_opt.dielectric);
  CSB.setUseRdielectric(_opt.distanceDependentElectrostatics);
  CSB.setVdwRescalingFactor(_opt.vdwScale);
  //CSB.setBuildNonBondedInteractions(false); // Don't build non-bonded terms.

  CSB.buildSystem(pseq);  // this builds atoms with emtpy coordinates. It also build bonds,angles and dihedral energy terms in the energyset (energyset lives inside system).


  // Apply coordinates from structure from PDB
  // Atom map to convert from PDB/ROSETTA naming to CHARMM.
  map<string,string> atom_map;
  atom_map["ILE:CD1"] = "CD";
  atom_map["ILE:1HD1"] = "HD1";
  atom_map["ILE:2HD1"] = "HD2";
  atom_map["ILE:3HD1"] = "HD3";
  atom_map["CYS:HG"] = "HG1";
  atom_map["SER:HG"] = "HG1";

  atom_map["1HA"] = "HA1";
  atom_map["2HA"] = "HA2";
  atom_map["3HA"] = "HA3";
  atom_map["1HB"] = "HB1";
  atom_map["2HB"] = "HB2";
  atom_map["3HB"] = "HB3";
  atom_map["1HD"] = "HD1";
  atom_map["2HD"] = "HD2";
  atom_map["3HD"] = "HD3";
  atom_map["1HG"] = "HG1";
  atom_map["2HG"] = "HG2";
  atom_map["3HG"] = "HG3";
  atom_map["1HE"] = "HE1";
  atom_map["2HE"] = "HE2";
  atom_map["3HE"] = "HE3";
  atom_map["1HZ"] = "HZ1";
  atom_map["2HZ"] = "HZ2";
  atom_map["3HZ"] = "HZ3";
  atom_map["H"] = "HN";
  atom_map["1HD1"] = "HD11";
  atom_map["2HD1"] = "HD12";
  atom_map["3HD1"] = "HD13";
  atom_map["1HE1"] = "HE11";
  atom_map["2HE1"] = "HE12";
  atom_map["3HE1"] = "HE13";
  atom_map["1HG1"] = "HG11";
  atom_map["2HG1"] = "HG12";
  atom_map["3HG1"] = "HG13";
  atom_map["1HD2"] = "HD21";
  atom_map["2HD2"] = "HD22";
  atom_map["3HD2"] = "HD23";
  atom_map["1HE2"] = "HE21";
  atom_map["2HE2"] = "HE22";
  atom_map["3HE2"] = "HE23";
  atom_map["1HG2"] = "HG21";
  atom_map["2HG2"] = "HG22";
  atom_map["3HG2"] = "HG23";
  atom_map["1HD3"] = "HD31";
  atom_map["2HD3"] = "HD32";
  atom_map["3HD3"] = "HD33";
  atom_map["1HE3"] = "HE31";
  atom_map["2HE3"] = "HE32";
  atom_map["3HE3"] = "HE33";
  atom_map["1HG3"] = "HG31";
  atom_map["2HG3"] = "HG32";
  atom_map["3HG3"] = "HG33";
  atom_map["1HH1"] = "HH11";
  atom_map["2HH1"] = "HH12";
  atom_map["1HH2"] = "HH21";
  atom_map["2HH2"] = "HH22";

  atom_map["1H"] = "HT1";
  atom_map["2H"] = "HT2";
  atom_map["3H"] = "HT3";

  int numAssignedAtoms = ene_sys.assignCoordinates(_sys.getAtomPointers(),&atom_map,false);
  //fprintf(stdout, "Assigned %8d atoms\n",numAssignedAtoms);

  ene_sys.buildAllAtoms();

  CSB.updateNonBonded(9.0, 10.0, 11.0);
  //CSB.updateNonBonded(opt.cuton,opt.cutoff,opt.cutnb);

  HydrogenBondBuilder hb;
  if(_opt.hbondfile != "") {
    hb.setSystem(ene_sys);
    if(!hb.readParameters(_opt.hbondfile)) {
      cerr << "ERROR 1234 Unable to read hbondfile " << _opt.hbondfile << endl;
      exit(1234);
    }
    hb.buildInteractions(-1.0);
  }

  BaselineEnergyBuilder bb;
  if(_opt.baselinefile != "") {
    bb.setSystem(ene_sys);
    if(!bb.readParameters(_opt.baselinefile)) {
      std::cerr << "Unable to read baselineenergies " <<  _opt.baselinefile << endl;
      exit(0);
    }

    if(!bb.buildInteractions()) {
      std::cerr << "Unable to build baselineInteractions " <<  endl;
      exit(0);
    }
  }

  // Write initially built file.
  std::string filename = "/tmp/initialBuild.pdb";
  //cout << "Write pdb " << filename << std::endl;
  ene_sys.writePdb(filename);


  CharmmEnergyCalculator calc(_opt.parfile);
  calc.setEnergyByType(true);
  calc.extractBondedInteractions(*ene_sys.getEnergySet());
  calc.setDielectricConstant(_opt.dielectric);
  calc.setUseRdielectric(_opt.distanceDependentElectrostatics);
  calc.setVdwRescalingFactor(_opt.vdwScale);
  calc.setNonBondedCutoffs(_opt.cuton,_opt.cutoff);

  for (uint i = 0; i < ene_sys.positionSize();i++){

    // Clear Energy-by-type in CharmmEnergyCalculator
    calc.clearEnergiesByType();

    // Calculate self and template energy
    double totalE = calc.calculateSelfEnergy(ene_sys,i,0) + calc.calculateTemplateEnergy(ene_sys,i,0,true);

    map<string,double> energiesByPosition = calc.getAllComputedEnergiesByType();

    map<string,double>::iterator it;
    if (i == 0){
      fprintf(stdout, "      %-30s ","");
      for (it = energiesByPosition.begin();it != energiesByPosition.end();it++){
	fprintf(stdout, "%15s", it->first.c_str());
      }
      fprintf(stdout,"\n");
    }
    fprintf(stdout, "ENERGY %-30s",ene_sys.getPosition(i).getPositionId().c_str());
    for (it = energiesByPosition.begin();it != energiesByPosition.end();it++){
      fprintf(stdout, "%15.3f",it->second);
    }
    fprintf(stdout,"\n");

    if (i == ene_sys.positionSize()-1){
      ene_sys.calcEnergy();
      fprintf(stdout, "ENERGY %-30s","TOTAL");
      totalE = 0.0;
      for (it = energiesByPosition.begin();it != energiesByPosition.end();it++){
	if (it->first != "TOTAL") {
	  fprintf(stdout, "%15.3f", ene_sys.getEnergySet()->getTermEnergy(it->first));
	  totalE += ene_sys.getEnergySet()->getTermEnergy(it->first);
	}
      }
      fprintf(stdout, "%15.3f", totalE);
      fprintf(stdout,"\n");
      
    }
   
  }


}

void reportStructureMotifMetrics(System &_sys, Options &_opt){ 
}

void reportSequenceMotifMetrics(System &_sys, Options &_opt){ 

  if (_opt.prosite == "") return;

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
		cout << "Usage: designCheck " << endl;
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
	  opt.prosite = SYSENV.getEnv("MSL_PROSITE");
	  if (opt.prosite == "UNDEF"){
	    cerr << "ERROR 1111 prosite undefined"<<endl;
	   exit(1111);
	  } else {
	    cerr << "WARNING prosite defaulted to: "<<opt.prosite<<endl;
	  }
	}
	opt.ref = OP.getString("ref");

	opt.topfile = OP.getString("topfile");
	if (OP.fail()){
	    opt.topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
	    if (opt.topfile == "UNDEF"){
	      cerr << "ERROR 1111 topfile not defined\n";
	      exit(1111);
	    } else {
	      cerr << "WARNING topfile defaulted to: "<<opt.topfile<<endl;
	    }
	}

	opt.parfile = OP.getString("parfile");
	if (OP.fail()){
		opt.parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
		if (opt.parfile == "UNDEF"){
		  cerr << "ERROR 1111 parfile not defined\n";
		  exit(1111);
		}else {
		  cerr << "WARNING parfile defaulted to: "<<opt.parfile<<endl;
		}
	}

	opt.hbondfile = OP.getString("hbondfile");
	if (OP.fail()){
		opt.hbondfile = SYSENV.getEnv("MSL_HBOND_PAR");
		if (opt.hbondfile == "UNDEF"){
		  cerr << "WARNING 1111 hbondfile not defined - not building hydrogen bond interactions\n";
		  opt.hbondfile = "";
		}else {
		  cerr << "WARNING hbondfile defaulted to: "<<opt.hbondfile<<endl;
		}
	}


	opt.baselinefile = OP.getString("baselinefile");
	if (OP.fail()){
		opt.baselinefile = "";
		cerr << "WARNING no baselinefile specified " <<endl;
	}

	opt.dielectric = OP.getDouble("dielectric");
	if (OP.fail()){
	  opt.dielectric = 80;
	}

	opt.distanceDependentElectrostatics = OP.getBool("distanceDependentElectrostatics");
	if (OP.fail()){
	  opt.distanceDependentElectrostatics = false;
	}

	opt.vdwScale = OP.getDouble("vdwScale");
	if (OP.fail()){
	  opt.vdwScale = 1.0;
	}

	opt.cuton = OP.getDouble("cuton");
	if (OP.fail()){
	  opt.cuton = 9.0;
	}
	opt.cutoff = OP.getDouble("cutoff");
	if (OP.fail()){
	  opt.cutoff = 15.0;
	}

	return opt;
}
