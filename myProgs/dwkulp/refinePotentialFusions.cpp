
#include "MslTools.h"
#include "OptionParser.h"
#include "AtomSelection.h"
#include "Transforms.h"
#include "Position.h"
#include "release.h"
#include "refinePotentialFusions.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
#include "MslOut.h"
#include "PDBFragments.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>

using namespace std;
using namespace MSL;
using namespace MslTools;

static MslOut MSLOUT("refinePotentialFusions");

int numberOfMainChainClashes(AtomPointerVector &_partner, AtomPointerVector  &_ref, int _clashTolerance=-1, int _removePartnerRes1=0, int _removeParterRes2=0,int _removeRefRes1=0,int _removeRefRes2=0);

int main(int argc, char *argv[]) {

    Options opt = setupOptions(argc, argv);

    // IO stuff
    System pdb;
    pdb.readPdb(opt.pdb);

    // Select scaffold chains that will move
    AtomSelection pdbSel(pdb.getAtomPointers());
    AtomPointerVector scaffoldChains = pdbSel.select("scaffold,chain X+Y+Z");
    AtomPointerVector notScaffoldChains = pdbSel.select("notScaffold,not chain X+Y+Z");

    // Select 4 junction positions 
    /*
      2 from chains A B or C: opt.junction_positionABC_1 , opt.junction_positionABC_2
      2 from chains X Y or Z: opt.junction_positionXYZ_1 , opt.junction_positionXYZ_2
     */
    // Test existence of junction positions
    if (!pdb.positionExists(opt.junction_positionABC_1)) cerr << "ERROR position_junctionABC_1 does not exist in "<<opt.pdb<<endl;
    if (!pdb.positionExists(opt.junction_positionABC_2)) cerr << "ERROR position_junctionABC_2 does not exist in "<<opt.pdb<<endl;
    if (!pdb.positionExists(opt.junction_positionXYZ_1)) cerr << "ERROR position_junctionXYZ_1 does not exist in "<<opt.pdb<<endl;
    if (!pdb.positionExists(opt.junction_positionXYZ_2)) cerr << "ERROR position_junctionXYZ_2 does not exist in "<<opt.pdb<<endl;
    Position &pABC_1  = pdb.getPosition(opt.junction_positionABC_1);
    Position &pABC_2  = pdb.getPosition(opt.junction_positionABC_2);
    Position &pXYZ_1  = pdb.getPosition(opt.junction_positionXYZ_1);
    Position &pXYZ_2  = pdb.getPosition(opt.junction_positionXYZ_2);
    
    Atom &pABC_1_CA  =pABC_1.getAtom("CA");
    Atom &pXYZ_1_CA = pXYZ_1.getAtom("CA");
    Atom &pABC_2_CA  =pABC_2.getAtom("CA");
    Atom &pXYZ_2_CA = pXYZ_2.getAtom("CA");

    cout << "SAMPLE DOF"<<endl;
    int minClashes = MslTools::intMax;

    double D1 = pABC_1_CA.distance2(pXYZ_1_CA);
    double D2 = pABC_2_CA.distance2(pXYZ_2_CA);

    fprintf(stdout,"Starting Conformation: [ %8.3f , %8.3f ]\n",sqrt(D1),sqrt(D2));

    fflush(stdout);
    // Transforms object 
    Transforms move;

    // For a range of Z-values 	 (+/- 5 Angstroms)
    CartesianPoint z_coord(0.0,0.0,1.0);

    std::priority_queue< std::pair<int,string>, std::vector< std::pair<int,string> > > lowClashModels;

    // Sample DOFs
    scaffoldChains.saveCoor("pre");
    for (double transZ = -opt.transZ; transZ < opt.transZ; transZ += 1.0){
      //cout << "Trans: "<<transZ<<endl;
      // Translate
      CartesianPoint transVec = z_coord * transZ;

      // Translate atoms..
      move.translate(scaffoldChains,transVec);

      scaffoldChains.saveCoor("post_transZ");
      for (double rotZ = -opt.rotZ; rotZ < opt.rotZ; rotZ += 1.0){
	//cout << "Rot: "<<rotZ<<endl;

	// Rotate
	move.Zrotate(scaffoldChains,rotZ);

	double D1 = pABC_1.getAtom("CA").distance2(pXYZ_1.getAtom("CA"));
	double D2 = pABC_2.getAtom("CA").distance2(pXYZ_2.getAtom("CA"));
	//fprintf(stdout,"Next Conformation: [ %8.3f , %8.3f ]\n%s\n%s\n",sqrt(D1),sqrt(D2),pABC_1_CA.toString().c_str(),pXYZ_1_CA.toString().c_str());
	
	// Filter for junction points being within defined limits
	if (!
	    (D1 >= opt.min_junction_distance &&
	     D1 <= opt.max_junction_distance &&
	     D2 >= opt.min_junction_distance &&
	     D2 <= opt.max_junction_distance)){
	  scaffoldChains.applySavedCoor("post_transZ");	
	  continue;
	}


	// Clash Check between scaffoldChains and notScaffoldChains
	int numClashes = numberOfMainChainClashes(scaffoldChains,notScaffoldChains,minClashes,pXYZ_1.getResidueNumber(),pXYZ_2.getResidueNumber());
	//cout << "NUM CLASHES: "<<numClashes<<" "<<minClashes<<endl;
	if (lowClashModels.size() < opt.numLowClashModels || numClashes < lowClashModels.top().first){

	  fprintf(stdout,"New conformation: %3d clashes %8.3f transZ %8.3f rotZ [ %8.3f , %8.3f ]",numClashes,transZ,rotZ,sqrt(D1),sqrt(D2));

	  if (numClashes < minClashes){
	    minClashes = numClashes;
	    scaffoldChains.saveCoor("minClashes");
	  }

	  if (lowClashModels.size() == opt.numLowClashModels){
	    fprintf(stdout, " -- previous top -- %s %d", lowClashModels.top().second.c_str(),lowClashModels.top().first);
	    scaffoldChains.clearSavedCoor(lowClashModels.top().second);
	    lowClashModels.pop();
	  } 
	  fprintf(stdout,"\n");
	  fflush(stdout);
	  string modelLabel = MslTools::stringf("%8.2f:%8.2f",transZ,rotZ);
	  scaffoldChains.saveCoor(modelLabel);
	  lowClashModels.push(pair<int,string>(numClashes, modelLabel));
	}

	scaffoldChains.applySavedCoor("post_transZ");	
      }
      scaffoldChains.applySavedCoor("pre");
    }
    
    if (minClashes == MslTools::intMax) {
      cerr << "ERROR 9999 No conformations found within junction distance\n";
      exit(9999);
    }
    /*
    scaffoldChains.applySavedCoor("minClashes");

    string fname = MslTools::stringf("minClashes_%s.pdb", MslTools::getFileName(opt.pdb).c_str());
    pdb.writePdb(fname);
    fprintf(stdout, "%-50s %6d\n",fname.c_str(),minClashes);
    */

    // Load PDB Fragment dataset
    MSLOUT.stream() << "Loading pdb fragment database..."<<endl;
    PDBFragments frag(opt.fragdb);
    frag.setPdbDir(opt.pdbdir);
    frag.loadFragmentDatabase();

    vector<string> stem1positions;
    vector<string> stem2positions;
    /*
    stem1positions.push_back(opt.junction_positionXYZ_1);
    stem1positions.push_back(opt.junction_positionXYZ_2);
    stem2positions.push_back(opt.junction_positionABC_1);
    stem2positions.push_back(opt.junction_positionABC_2);
    */

    
    int xyz_1_index = pdb.getPositionIndex(&pXYZ_1);
    int xyz_2_index = pdb.getPositionIndex(&pXYZ_2);
    int abc_1_index = pdb.getPositionIndex(&pABC_1);
    int abc_2_index = pdb.getPositionIndex(&pABC_2);
	
    int modelNum = 1;
    while (!lowClashModels.empty()){
      scaffoldChains.applySavedCoor(lowClashModels.top().second);
      lowClashModels.pop();
     

      MSLOUT.stream() << "Searching pdb fragment database... for model "<<modelNum<<endl;
      string fname = MslTools::stringf("lowClashes_%06d.pdb",modelNum);
      pdb.writePdb(fname);
      fprintf(stdout, "%-50s\n",fname.c_str());


      // +/-2 positions
      for (int p1 = -2; p1 <= 2;p1++){
	string pos1 = MslTools::getPositionId(pdb.getPosition(xyz_1_index).getChainId(),pdb.getPosition(xyz_1_index).getResidueNumber()+p1,"");
	if (! pdb.positionExists(pos1)) continue;
	if ( pdb.getPosition(pos1).getChainId() != pXYZ_1.getChainId()) continue;

	for (int p2 = -2; p2 <= 2;p2++){

	  string pos2 = MslTools::getPositionId(pdb.getPosition(xyz_2_index).getChainId(),pdb.getPosition(xyz_2_index).getResidueNumber()+p2,"");
	  if (! pdb.positionExists(pos2)) continue;
	  if ( pdb.getPosition(pos2).getChainId() != pXYZ_2.getChainId()) continue;

	  for (int p3 = -2; p3 <= 2;p3++){

	    string pos3 = MslTools::getPositionId(pdb.getPosition(abc_1_index).getChainId(),pdb.getPosition(abc_1_index).getResidueNumber()+p3,"");
	    if (! pdb.positionExists(pos3)) continue;
	    if ( pdb.getPosition(pos3).getChainId() != pABC_1.getChainId()) continue;

	    for (int p4 = -2; p4 <= 2;p4++){

	      string pos4 = MslTools::getPositionId(pdb.getPosition(abc_2_index).getChainId(),pdb.getPosition(abc_2_index).getResidueNumber()+p4,"");
	      if (! pdb.positionExists(pos4)) continue;
	      if ( pdb.getPosition(pos4).getChainId() != pABC_2.getChainId()) continue;

	      MSLOUT.stream() << "Search using stems: "<<pos1<<" "<<pos2<<" "<<pos3<<" "<<pos4<<endl;
	      stem1positions.clear();
	      stem1positions.push_back(pos1);
	      stem1positions.push_back(pos2);

	      stem2positions.clear();
	      stem2positions.push_back(pos3);
	      stem2positions.push_back(pos4);
	      
	      int numHits = 0;
	      try {
		numHits = frag.searchForMatchingDualFragments(pdb,stem1positions,
							      pdb,stem2positions,
							      opt.loop1min,opt.loop1max,opt.loop2min,opt.loop2max,opt.distanceStem1,opt.distanceStem2, opt.stemRmsdTol,opt.totalRmsdTol,false);
	      } catch (std::exception &e){
		MSLOUT.stream() << e.what()<<endl;
	      }
	      MSLOUT.fprintf(stdout, "\thits = %10d\n",numHits);

	      vector<AtomContainer *> &dualLoops = frag.getAtomContainers();
	      for (uint h = 0; h < dualLoops.size();h++){
		string loopname=MslTools::stringf("%s_%05d_%d_%d_%d_%d.pdb",MslTools::getFileName(fname).c_str(),h+1,p1,p2,p3,p4);
		dualLoops[h]->writePdb(loopname);
	      }
	    } // FOR p4
	  } // FOR p3
	} // FOR p2
      } // FOR p1





      modelNum++;
    }
}

Options setupOptions(int theArgc, char * theArgv[]){
    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: refinePotentialFusions " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdb PDB\n";
	cout << "junction_positionABC_1 PosId\n";
	cout << "junction_positionABC_2 PosId\n";
	cout << "junction_positionXYZ_1 PosId\n";
	cout << "junction_positionXYZ_2 PosId\n";
	cout << "min_junction_distance DIST\n";
	cout << "max_junction_distance DIST\n";
	cout << endl;
	exit(0);
    }

    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdb specified."<<endl;
	exit(1111);
    }

    opt.junction_positionABC_1 = OP.getString("junction_positionABC_1");
    if (OP.fail()){
	cerr << "ERROR 1111 no junction_positionABC_1 specified."<<endl;
	exit(1111);
    }    
    opt.junction_positionABC_2 = OP.getString("junction_positionABC_2");
    if (OP.fail()){
	cerr << "ERROR 1111 no junction_positionABC_1 specified."<<endl;
	exit(1111);
    }    
    opt.junction_positionXYZ_1 = OP.getString("junction_positionXYZ_1");
    if (OP.fail()){
	cerr << "ERROR 1111 no junction_positionXYZ_1 specified."<<endl;
	exit(1111);
    }    
    opt.junction_positionXYZ_2 = OP.getString("junction_positionXYZ_2");
    if (OP.fail()){
	cerr << "ERROR 1111 no junction_positionXYZ_1 specified."<<endl;
	exit(1111);
    }    

    opt.min_junction_distance = OP.getDouble("min_junction_distance");
    if (OP.fail()){
      opt.min_junction_distance = 5;
      cerr << "WARNING default min_junction_distance of "<<opt.min_junction_distance<<" is being used."<<endl;
    }

    opt.max_junction_distance = OP.getDouble("max_junction_distance");
    if (OP.fail()){
      opt.max_junction_distance = 15;
      cerr << "WARNING default max_junction_distance of "<<opt.max_junction_distance<<" is being used."<<endl;
    }

    // Do seach in Square Angstroms
    opt.min_junction_distance *= opt.min_junction_distance;
    opt.max_junction_distance *= opt.max_junction_distance;

    opt.transZ = OP.getDouble("transZ_limit");
    if (OP.fail()){
      opt.transZ = 10;
      cerr << "WARNING default transZ to "<<opt.transZ<<endl;
    }
    opt.rotZ = OP.getDouble("rotZ_limit");
    if (OP.fail()){
      opt.rotZ = 10;
      cerr << "WARNING default rotZ to "<<opt.rotZ<<endl;
    }    

    opt.numLowClashModels = OP.getInt("numLowClashModels");
    if (OP.fail()){
      opt.numLowClashModels = 1;
    }

    

	opt.distanceStem1 = OP.getDouble("distanceStem1");
	if (OP.fail()){
		cerr << "WARNING 1111 no distanceStem1 specified."<<endl;
		opt.distanceStem1 = 0.0;
	}

	opt.distanceStem2 = OP.getDouble("distanceStem2");
	if (OP.fail()){
		cerr << "WARNING 1111 no distanceStem2 specified."<<endl;
		opt.distanceStem2 = 0.0;
	}


	opt.loop1min = OP.getInt("loop1min");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop1min specified."<<endl;
		exit(1111);
	}

	opt.loop1max = OP.getInt("loop1max");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop1max specified."<<endl;
		exit(1111);
	}

	opt.loop2min = OP.getInt("loop2min");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop2min specified."<<endl;
		exit(1111);
	}

	opt.loop2max = OP.getInt("loop2max");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop2max specified."<<endl;
		exit(1111);
	}

	opt.stemRmsdTol = OP.getDouble("stemRmsdTol");
	if (OP.fail()){
	  opt.stemRmsdTol = 0.3;
	}
	opt.totalRmsdTol = OP.getDouble("totalRmsdTol");
	if (OP.fail()){
	  opt.totalRmsdTol = 0.3;
	}
	opt.fragdb = OP.getString("fragDB");
	opt.pdbdir = OP.getString("pdbdir");

    cout << OP<<endl;
    return opt;
}

int numberOfMainChainClashes(AtomPointerVector &_partner, AtomPointerVector  &_ref, int _clashTolerance,int _removePartnerRes1, int _removePartnerRes2,int _removeRefRes1, int _removeRefRes2){

  int clashes = 0;
  for (uint a = 0 ; a  < _partner.size();a++){
    if (!
	(_partner(a).getName() == "N" ||
	 _partner(a).getName() == "CA" ||
	 _partner(a).getName() == "C" ||
	 _partner(a).getName() == "CB")
	) continue;

    // Don't count clashes in _removePartnerRes1 - _removePartnerRes2 region
    if ( _removePartnerRes1 != 0 && _partner(a).getResidueNumber() > _removePartnerRes1 && _removePartnerRes2 != 0 && _partner(a).getResidueNumber() < _removePartnerRes2) continue;
	

    for (uint b = 0 ; b  < _ref.size();b++){
    if (!
	(_ref(b).getName() == "N" ||
	 _ref(b).getName() == "CA" ||
	 _ref(b).getName() == "C" ||
	 _ref(b).getName() == "CB")
	) continue;

      // Don't count clashes in _removeRefRes1 - _removeRefRes2 region
      //if ( _removeRefRes1 != 0 && _ref(b).getResidueNumber() > _removeRefRes1 && _removeRefRes2 != 0 && _ref(b).getResidueNumber() < _removeRefRes2) continue;

      double distSq = _partner(a).distance2(_ref(b));

      // < 3.0 Angstroms is a clash
      if (distSq < 9.00){
	clashes++;
      }

      if (_clashTolerance != -1 && clashes > _clashTolerance){
	return clashes;
      }
    }
  }

  return clashes;
}
