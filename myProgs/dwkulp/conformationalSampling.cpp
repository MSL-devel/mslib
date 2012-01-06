#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <iterator>

#include "AtomBondBuilder.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "conformationalSampling.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("conformationalSampling");

// SysEnv
static SysEnv SYSENV;



int main(int argc, char *argv[]) {



	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnOn("conformationalSampling");

	// Read in the pdb with a set of functional groups defined
	MSLOUT.stream() << "Read in "<<opt.pdb<<endl;
	System sys;
	sys.readPdb(opt.pdb);

	MSLOUT.stream() << "AtomBondBuilder"<<endl;
	AtomBondBuilder abb;
	abb.buildConnections(sys.getAtomPointers());

	MSLOUT.stream() << "Ref Sys: "<<opt.ref<<endl;
	System ref;
	AtomPointerVector refAts;
	AtomPointerVector sysAts;
	if (opt.ref != ""){
	  ref.readPdb(opt.ref);

	  // Each atom
	  vector<string> refIds = MslTools::tokenize(opt.refSelect,";");
	  if (refIds.size() != 1){
	    for (uint i = 0 ; i < refIds.size();i++){
	      if (ref.atomExists(refIds[i])){
		refAts.push_back(&ref.getAtom(refIds[i]));
	      } else {
		MSLOUT.stream() << "ERROR ref atom "<<refIds[i]<<" does not exist"<<endl;
	      }
	    }
	  } else {

	    // Assume it is a regular selection statement
	    AtomSelection refSel(ref.getAtomPointers());
	    refAts = refSel.select(opt.refSelect);

	  }

	}

	if (opt.pdbSelect != ""){
	  vector<string> pdbIds = MslTools::tokenize(opt.pdbSelect,";");
	  if (pdbIds.size() != 1){
	    for (uint i = 0 ; i < pdbIds.size();i++){
	      if (sys.atomExists(pdbIds[i])){
		sysAts.push_back(&sys.getAtom(pdbIds[i]));
	      } else {
		MSLOUT.stream() << "ERROR pdb atom "<<pdbIds[i]<<" does not exist"<<endl;
	      }
	    }
	  } else {

	    // Assume it is a regular selection statement
	    AtomSelection pdbSel(sys.getAtomPointers());
	    sysAts = pdbSel.select(opt.pdbSelect);
	  }
	}

	if (opt.rmsd != -1.0 && refAts.size() != sysAts.size()){
	  cerr << "ERROR 2222 reference and system atom vector sizes are not equal: "<<refAts.size()<<","<<sysAts.size()<<" from "<<opt.ref<<","<<opt.pdb<<endl;
	  exit (2222);
	}

	cout << "refAts: "<<endl<<refAts<<endl;
	cout << "sysAts: "<<endl<<sysAts<<endl;

	vector<vector<DoF> > allValues;
	int maxValues = 0;
	for (uint i=0; i<opt.DoFs.size(); i++) {
	  MSLOUT.stream()<< "DoFs.size() == "<<opt.DoFs.size()<<" and i = "<<i<<endl;
	  // Make sure atoms exist..
	  for (uint a = 0; a < opt.DoFs[i].atomNames.size();a++){
	    if (!sys.atomExists(opt.DoFs[i].atomNames[a])){
	      cerr << "ERROR 2222 Atom : "<<opt.DoFs[i].atomNames[a]<<" does not exist in "<<opt.pdb<<endl;
	      exit(2222);
	    }
	  }

	  vector<DoF> allValue;
	  for (double j = opt.DoFs[i].minValue; j <= opt.DoFs[i].maxValue;j += opt.DoFs[i].stepSize){
	    DoF tmp;
	    tmp.atomNames = opt.DoFs[i].atomNames;
	    tmp.value = j;
	    allValue.push_back(tmp);
	  }
	  if (allValue.size() > maxValues){
	    maxValues = allValue.size();
	  }
	  allValues.push_back(allValue);
	}

	MSLOUT.stream()<< "Do Transforms: "<<allValues.size()<<endl;
	Transforms T;
	int comboIndex = 1;
	for (int i = 0; i < allValues[0].size();i++){
	  
	      T.setDihedral(sys.getAtom(allValues[0][i].atomNames[0]),
			  sys.getAtom(allValues[0][i].atomNames[1]),
			  sys.getAtom(allValues[0][i].atomNames[2]),
			  sys.getAtom(allValues[0][i].atomNames[3]),
			  allValues[0][i].value);	      

	  for (int j = 0; j < allValues[1].size();j++){

	      T.setDihedral(sys.getAtom(allValues[1][j].atomNames[0]),
			  sys.getAtom(allValues[1][j].atomNames[1]),
			  sys.getAtom(allValues[1][j].atomNames[2]),
			  sys.getAtom(allValues[1][j].atomNames[3]),
			  allValues[1][j].value);	      


	    for (int k = 0; k < allValues[2].size();k++){

	      T.setDihedral(sys.getAtom(allValues[2][k].atomNames[0]),
			  sys.getAtom(allValues[2][k].atomNames[1]),
			  sys.getAtom(allValues[2][k].atomNames[2]),
			  sys.getAtom(allValues[2][k].atomNames[3]),
			  allValues[2][k].value);	      

	      

	      if (opt.rmsd != -1.0) {
		double rmsd = refAts.rmsd(sysAts);
		if (rmsd < opt.rmsd){
		  string fileName =MslTools::stringf("%s_%010d.pdb",opt.outPdb.c_str(),comboIndex++);	    
		  fprintf(stdout,"%s %8.3f %8.3f %8.3f %8.3f\n",fileName.c_str(),rmsd,allValues[0][i].value,allValues[1][j].value,allValues[2][k].value);
		  sys.writePdb(fileName);
		}
	      } 

	      if (opt.clashDist != -1.0){
		// For each sysAts and refAts look for clash 
		string clashPdb = "";
		string clashRef = "";
		double clashDist = 0.0;
		bool clash = false;
		for (uint c1 = 0; c1 < sysAts.size();c1++){
		  if (sysAts[c1]->getElement() == "H") continue;
		  for (uint c2 = 0; c2 < refAts.size();c2++){
		    if (refAts[c2]->getElement() == "H") continue;
		    double dist = sysAts[c1]->distance(*refAts[c2]);
		    if (dist <= opt.clashDist){
		      clashPdb = sysAts[c1]->getAtomId();
		      clashRef = refAts[c2]->getAtomId();
		      clashDist = dist;
		      clash = true;
		      break;
		    }
		  }
		  if (clash) break;
		}

		if (clash){
		  fprintf(stdout, "clash found: %s,%s,%8.3f\n",clashPdb.c_str(),clashRef.c_str(),clashDist);
		  continue;
		} else {
		  string fileName = MslTools::stringf("%s_%010d.pdb",opt.outPdb.c_str(),comboIndex++);	    
		  fprintf(stdout,"MODEL %s %8.3f %8.3f %8.3f\n",fileName.c_str(),allValues[0][i].value,allValues[1][j].value,allValues[2][k].value);
		  sys.writePdb(fileName);
		}
	      }

	      
	    }
	  }
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
		cout << "conformationalSampling --pdb PDB --atomNames  A,23,CA A,23,CB A,23,CG A,23,CD1 --range=-180,180,5\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	vector<string> atomNames = OP.getMultiString("atomNames");
	if (OP.fail()){
	  cerr << "ERROR 1111 no atomNames specified"<<endl;
	  exit(1111);
	}
	vector<string> range = OP.getMultiString("range");
	if (OP.fail()){
	  cerr << "ERROR 1111 no range specified"<<endl;
	  exit(1111);

	}

	for (uint i = 0; i < atomNames.size();i++){
	  cout << "atomNames["<<i<<"]: "<<atomNames[i]<<endl;
	  cout << "range["<<i<<"]: "<<range[i]<<endl;
	  DoFSpec tmp;
	  tmp.atomNames = MslTools::tokenize(atomNames[i]," ");
	  vector<string> vrange = MslTools::tokenize(range[i],",");
	  cout << "Size check"<<endl;
	  if (vrange.size() != 3){
	    cerr << "ERROR 1111 range specification must be 3 values (min,max,step size) : "<<vrange.size()<<" "<<range[i]<<endl;
	    exit(1111);
	  }
	  cout << "toDouble"<<endl;
	  tmp.minValue = MslTools::toDouble(vrange[0]);
	  tmp.maxValue = MslTools::toDouble(vrange[1]);
	  tmp.stepSize = MslTools::toDouble(vrange[2]);

	  opt.DoFs.push_back(tmp);
	  cout << "Done"<<endl;
	  cout << "DoFspec: "<<tmp.toString()<<endl;
	} 

	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  vector<string> nameToks = MslTools::tokenize(MslTools::getFileName(opt.pdb),".");
	  opt.outPdb = nameToks[0]; 
	  cerr << "WARNING outPdb is defaulted to "<<opt.outPdb<<endl;
	}

	opt.ref = OP.getString("ref");
	opt.refSelect = OP.getString("refSelect");
	opt.pdbSelect = OP.getString("pdbSelect");

	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  opt.rmsd = -1.0;
	}
	opt.clashDist = OP.getDouble("clashDist");
	if (OP.fail()){
	  opt.clashDist = -1.0;
	}
	
	return opt;
}



