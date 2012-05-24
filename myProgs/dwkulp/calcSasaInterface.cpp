#include <iostream>
#include <cstdlib>
#include <fstream>


#include "MslTools.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Timer.h"
#include "SasaCalculator.h"
#include "Position.h"
#include "PDBTopology.h"
#include "VectorPair.h"
#include "calcSasaInterface.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <numeric>
#include <vector>

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("calcSasaInterface");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
	// MSLOUT.turnOff("calcSasaInterface");

	Timer t;
	double start = t.getWallTime();

	// Read in the pdb list
	vector<string> list;
	MslTools::readTextFile(list, opt.list);


	// first string is list[i] name, second is posId, double is tmpSasa - sasaRef.
	map<string, map<string, double> > sasaData;
	map<string, map<string, bool> > sasaResDiff;
	for (uint i = 0; i < list.size();i++){

	  System tmp;
	  tmp.readPdb(list[i]);
	  MSLOUT.stream() << "File : "<<MslTools::getFileName(list[i])<<" has "<<tmp.chainSize()<<" chains."<<endl;
	  //if (tmp.chainSize() != 2) continue;

	  fstream resfile;
	  if (opt.resfile){
	    resfile.open(MslTools::stringf("resfile_%s",MslTools::getFileName(list[i]).c_str()).c_str(),std::ios::out);
	    resfile << "NATRO"<<endl;
	    resfile << "start"<<endl;
	  }

	  SasaCalculator sasaCalc(tmp.getAtomPointers(), 1.4, 200);
	  sasaCalc.calcSasa();	  
	  double fullSasa = sasaCalc.getTotalSasa();

	  double chainSasa = 0.0;	  
	  bool doChainSasa  = true;
	  for (uint c = 0; c < tmp.chainSize();c++){

	    SasaCalculator sasaCalcMonomer(tmp.getChain(c).getAtomPointers(), 1.4, 200);
	    sasaCalcMonomer.calcSasa();
	    double monomerSasa = sasaCalcMonomer.getTotalSasa();

	    for (uint p = 0; p < tmp.getChain(c).positionSize();p++){
	      Position &pos1 = tmp.getChain(c).getPosition(p);
	      string posId = pos1.getPositionId();

	      double posSasa         = sasaCalcMonomer.getResidueSasa(posId) - sasaCalc.getResidueSasa(posId);
	      int cbDensityTotal     = 0;
	      int cbDensityInterface = 0;

	      vector<double> allTorsions;
	      vector<double> interfaceTorsions;

	      vector<double> allAngle1;
	      vector<double> interfaceAngle1;

	      vector<double> allAngle2;
	      vector<double> interfaceAngle2;

	      Atom *CB1 = NULL;

	      if (pos1.getResidueName() == "GLY"){
		CB1 = PDBTopology::getPseudoCbeta(pos1.getCurrentIdentity());
	      }
	      if (pos1.atomExists("CB")){
		CB1 = &pos1.getAtom("CB");
	      }
	      if (CB1 == NULL) continue;


	      if (!  ( pos1.atomExists("CA") && pos1.atomExists("N") && pos1.atomExists("C") ) ) continue;



	      // Number Cb within X distance in total, in other chains
	      for (uint c2 = 0; c2 < tmp.chainSize();c2++){
		if (doChainSasa){
		  SasaCalculator sasaCalcChain(tmp.getChain(c2).getAtomPointers(), 1.4, 200);
		  sasaCalcChain.calcSasa();
		  chainSasa += sasaCalcChain.getTotalSasa(); 
		}

		for (uint p2 = 0; p2 < tmp.getChain(c2).positionSize();p2++){
		  if (p == p2) continue;
		  Position &pos2 = tmp.getChain(c2).getPosition(p2);

		  Atom *CB2 = NULL;
		  if (pos2.getResidueName() == "GLY"){
		    CB2 = PDBTopology::getPseudoCbeta(pos2.getCurrentIdentity());
		  }
		  if (pos2.atomExists("CB")){
		    CB2 = &pos2.getAtom("CB");
		  }

		  // Must have CBetas defined...
		  if (CB2 == NULL) continue;
		  if (!  (pos2.atomExists("CA") && pos2.atomExists("N") && pos2.atomExists("C") ) ) continue;

		  double dist = CB1->distance(*CB2);
		  MSLOUT.stream() << "pos1,pos2,dist: "<<pos1.getPositionId()<<" ; "<<pos2.getPositionId()<<" ; "<<dist<<endl;
		  if (dist <= opt.dist) {
		    
		    cbDensityTotal++;

		    VectorPair vp(pos1.getAtom("CA").getCoor(),CB1->getCoor(), pos2.getAtom("CA").getCoor(),CB2->getCoor());
		    vp.calcAll();
		    allTorsions.push_back(vp.getTorsion());
		    allAngle1.push_back(vp.getAngle1());
		    allAngle2.push_back(vp.getAngle2());
		    if (c != c2){
		      cbDensityInterface++;
		      interfaceTorsions.push_back(vp.getTorsion());
		      interfaceAngle1.push_back(vp.getAngle1());
		      interfaceAngle2.push_back(vp.getAngle2());
		    }



		  }
		  
		  
		  if (pos2.getResidueName() == "GLY"){
		    delete(CB2);
		  }

		}// END p2		
	      } // END c2
	      doChainSasa = false;

	      if (pos1.getResidueName() == "GLY"){
		delete(CB1);
	      }


	      

	      // Print stuff out..
	      if (cbDensityInterface > 0){
		fprintf(stdout, "%s %s %s %4d%1s %8.1f %8.1f %8u %8u %8.2f %8.2f %8.2f %8.2f %8.2f  %8.2f %8.2f  %8.2f %8.2f  %8.2f %8.2f  %8.2f %8.2f \n", MslTools::getFileName(list[i]).c_str(),pos1.getChainId().c_str(), pos1.getResidueName().c_str(),pos1.getResidueNumber(), pos1.getResidueIcode().c_str(), 
			                                      chainSasa - fullSasa,
		                                              posSasa,
		                                              cbDensityTotal,
		                                              cbDensityInterface,
			                                      (double)cbDensityInterface/(double)cbDensityTotal,
		                                              mean(allTorsions), stddev(allTorsions),
		                                              mean(interfaceTorsions), stddev(interfaceTorsions),
		                                              mean(allAngle1), stddev(allAngle1),
		                                              mean(interfaceAngle1), stddev(interfaceAngle1),
		                                              mean(allAngle2), stddev(allAngle2),
		                                              mean(interfaceAngle2), stddev(interfaceAngle2));

	      }

	      if (opt.resfile){
		// Buried Interface position...
		if (cbDensityInterface > 4 || (double)cbDensityInterface/(double)cbDensityTotal > 0.25){

		    if ( ! (pos1.getResidueName() == "PRO" || pos1.getResidueName() == "GLY") ) {
		      resfile << "# Buried "<<cbDensityInterface<<endl;
		      resfile << pos1.getResidueNumber()<<pos1.getResidueIcode()<<"    "<<pos1.getChainId()<<" PIKAA "<<MslTools::getOneLetterCode(pos1.getResidueName())<<"PFWLMIVAGYT EX 1 LEVEL 4 EX 2 LEVEL 4 EX 3 EX 4 EX_CUTOFF 1 USE_INPUT_SC"<<endl;
		    }
		}

		// Edge Interface position...
//		if (cbDensityInterface > 1 && cbDensityInterface <= 4){
//		  resfile << "# Edge "<<cbDensityInterface<<endl;
//		  resfile << pos1.getResidueName()<<pos1.getResidueIcode()<<"    "<<pos1.getChainId()<<" PIKAA FWLIVYAGTNQHSRKED EX 1 LEVEL 1 EX 2 LEVEL 1 EX 3 EX 4 EX_CUTOFF 1 USE_INPUT_SC"<<endl;
//		}
	      }



	    } // END p

	  } // END c

	  if (opt.resfile){
	    resfile.close();
	  }
	}// END list


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
		cout << "calcSasaInterface --list LIST --dist DIST\n";
		exit(0);
	}
	opt.list= OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}
	opt.dist = OP.getDouble("dist");
	if (OP.fail()){
	  opt.dist = 8;
	  cerr << "WARNING distance cutoff defaulted to "<<opt.dist<<" Angstroms."<<endl;
	}
	opt.resfile = OP.getBool("resfile");

	return opt;
}


double mean(vector<double> &data){
   double tmp = accumulate( data.begin(), data.end(), 0.0f )/ data.size();

   return tmp;
}

double stddev(vector<double> &data){

   double tmp = mean(data);
   vector<double> zero_mean( data );
   transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( minus<double>(), tmp ) );

   double deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
   deviation = sqrt( deviation / ( data.size() - 1 ) );

   return deviation;
}
