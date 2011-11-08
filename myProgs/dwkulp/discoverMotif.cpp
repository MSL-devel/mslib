#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "PDBWriter.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "discoverMotif.h"
#include "VectorHashing.h"

using namespace std;
using namespace MSL;

// Build Inverse Rotamers from input PDB and collect the bounding boxes
void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers, vector<map<string,double> > &_boundingBoxes);
string extractFileName(string _input);
pair<string,string> extractPositions(string _input);

// Compute the geometry between bounding boxes;
void computeInterBoxDistances(vector<double> _minDistances, vector<double> _maxDistances,vector<map<string,double> > &_boundingBoxes);

// MslOut 
static MslOut MSLOUT("discoverMotif");

int main(int argc, char *argv[]) {

	// Parse commandline options
	Options opt = setupOptions(argc,argv);

        // MslOut can suppress output, to make example output clean
        MSLOUT.turnOn("discoverMotif");

	Timer t;
	double start = t.getWallTime();

	// Read in the pdb with a set of functional groups defined
	System sys;
	sys.readPdb(opt.pdb);

	// For each residue, build inverse rotamers
	vector<AtomContainer *> inverseRotamers;
	vector<map<string,double> > boundingBoxes;

	// Build inverse rotamers, get bounding boxes
	buildInverseRotamers(opt,sys,inverseRotamers,boundingBoxes);

	// Compute Inter bounding box distances
	vector<double> minDistances;
	vector<double> maxDistances;
	computeInterBoxDistances(minDistances,maxDistances,boundingBoxes);

	// Filter rotamer pairs
	// Keys = RotamerId = PositionId,CONF_NUM where CONF_NUM = rotamer number at this position.
	//map<string,map<string,bool> > rotamerPairsExlcudeList;
	//filterRotamerPairs(inverseRotamers,rotamerPairsExcludeList);

	// Compute rotamer pair vector geometry.
	//computeRotamerPairVectorGeometry(inverseRotamers);
	System inverseRotamerSystem;
	for (uint i = 0; i < inverseRotamers.size();i++){
		MSLOUT.stream() << "Adding "<<inverseRotamers[i]->size()<<" atoms at position "<<i<<endl;
		inverseRotamerSystem.addAtoms(inverseRotamers[i]->getAtomPointers());
	}

	if (opt.debug){
		MSLOUT.stream() << "System of inverse rotamers:"<<endl;
		MSLOUT.stream() << inverseRotamerSystem.toString()<<endl;
		for (uint p = 0; p < inverseRotamerSystem.positionSize();p++){
			Position &pos = inverseRotamerSystem.getPosition(p);
			MSLOUT.stream()<< "Position "<<pos.toString()<< " has "<<pos.getTotalNumberOfRotamers()<<" rotamers."<<endl;
			for (uint r = 0; r < pos.getTotalNumberOfRotamers();r++){
				stringstream ss;
				ss << "/tmp/inverseRotamer_"<<pos.getResidueNumber()<<"_"<<r<<".pdb";
				pos.setActiveRotamer(r);

				PDBWriter pdbW;
				pdbW.open(ss.str());
				pdbW.write(pos.getAtomPointers());
				pdbW.close();

			}
		}
	}

	double endTimeBuilding = t.getWallTime();
	MSLOUT.fprintf(stdout,"Time %8.3f for building rotamers\n",(endTimeBuilding - start));

	MSLOUT.turnOff("VectorHashing");
	// Write out VectorData ... end program .. new program read in vector data
	VectorHashing rotamerHash;
	MSLOUT.stream() << "Build VectorHash from inverse rotamers\n";
	rotamerHash.addToVectorHash(inverseRotamerSystem,"testingChainA",true);
	MSLOUT.stream() << "done building hash\n";
	double endTimeHash1 = t.getWallTime();
	MSLOUT.fprintf(stdout,"Time %8.3f for building hash with inverse rotamers\n",(endTimeHash1 - endTimeBuilding));


	VectorHashing testHash;
	testHash.load_checkpoint(opt.vh);

	MSLOUT.turnOn("VectorHashing");
	double endTimeHash2 = t.getWallTime();
	vector<map<string, vector<string> > > results = testHash.searchForVectorMatchAll(rotamerHash, inverseRotamerSystem.positionSize()-1);
	double endTimeSearch = t.getWallTime();
	MSLOUT.fprintf(stdout,"Time %8.3f for searching\n",(endTimeSearch - endTimeHash2));


	// Iterate over results, align onto input PDB
	int resultCount = 1;
	for (uint v = 0; v < results.size();v++){

	  map<string, vector<string> >::iterator cycleIt = results[v].begin();
	  for (;cycleIt != results[v].end();cycleIt++){
	    stringstream ss;
	    ss << "CYCLE: "<<cycleIt->first;

	    vector<string> inputPositionIds;
	    pair<string, vector<string> > resultPositionIds;
	    
	    string fileName = extractFileName(cycleIt->first);
	    pair<string,string>  positions = extractPositions(cycleIt->first);

	    resultPositionIds.first = fileName;
	    resultPositionIds.second.push_back(positions.first);
	    inputPositionIds.push_back(positions.second);

	    for (uint c = 0; c < cycleIt->second.size();c++){
	      ss << " AND "<< cycleIt->second[c];

	      pair<string,string> positions = extractPositions(cycleIt->second[c]);
	      resultPositionIds.second.push_back(positions.first);
	      inputPositionIds.push_back(positions.second);

	    }
	    cout << ss.str()<<endl;

	    // Get Input Positions from 'sys' object
	    AtomPointerVector inputAtoms;
	    for (uint p = 0; p < inputPositionIds.size();p++){
	      inputAtoms += sys.getPosition(inputPositionIds[p]).getAtomPointers();
	    }

	    AtomSelection selInput(inputAtoms);
	    AtomPointerVector inputCa = selInput.select("name N+CA+C");
	    
	    // Get Result Positions from PDB file
	    System resultSys;
	    resultSys.readPdb(resultPositionIds.first);
	    AtomPointerVector resultAtoms;
	    for (uint p = 0; p < resultPositionIds.second.size();p++){
	      resultAtoms += resultSys.getPosition(resultPositionIds.second[p]).getAtomPointers();
	    }

	    AtomSelection selResult(resultAtoms);
	    AtomPointerVector resultCa = selResult.select("name N+CA+C");

	    //cout << "USING ATOMS TO TRANSFORM1: \n"<<resultCa<<endl;
	    //cout << "USING ATOMS TO TRANSFORM2: \n"<<inputCa<<endl;
	    Transforms t;
	    if (!t.rmsdAlignment(resultCa, inputCa, resultSys.getAtomPointers())){
	      MSLOUT.stream() << "ERROR aligning resultCa and inputCa"<<endl;
	    }

	    double rmsd = resultCa.rmsd(inputCa);

	    if (rmsd < opt.rmsd){
	      MSLOUT.fprintf(stdout, "RMSD %8.3f\n",rmsd);

	      char tmp[100];
	      sprintf(tmp,"/tmp/motif-%010d-%s.pdb",resultCount,MslTools::getFileName(resultPositionIds.first).c_str());
	      stringstream ssTmp;
	      ssTmp << tmp;
	    
	      PDBWriter pout;
	      pout.open(ssTmp.str());
	      pout.write(resultSys.getAtomPointers());
	      pout.close();

	      resultCount++;
	    }

	  }

	}


	//map<string, map<string, vector<pair<string,string> > > > possibleCombinations;
	//filterPositionPairsForVectorDataMatch(VectorData,testSys, possibleCombinations);

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
		cout << "discoverMotif --pdb PDB --rotlib ROTLIB\n";
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.rotlib = OP.getString("rotlib");
	if (OP.fail()){
		cerr << "ERROR 1111 rotlib not specified.\n";
		exit(1111);
	}
	opt.numRotamers = OP.getInt("numRotamers");
	if (OP.fail()){
	  opt.numRotamers = 10;
	}
	opt.list = OP.getString("list");
	if (OP.fail()){
	  opt.list = "";
	}

	opt.vh = OP.getString("vh");
	if (OP.fail()){
	  opt.vh = "";
	}

	if (opt.list == opt.vh){
	  cerr << "ERROR 1111 Use EITHER --list OR --vh\n";
	  exit(1111);
	}

	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  cerr << "WARNING RMSD Tolerance set to 1.0\n";
	  opt.rmsd = 1.0;
	}
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}
	return opt;
}

void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers, vector<map<string,double> > &_boundingBoxes){
  
	// Create a pdb topology object
	PDBTopology pdbTop;

	// Set a rotamer library and use atoms defined in the rotamer library to build
	pdbTop.readRotamerLibrary(_opt.rotlib);
	pdbTop.setAddAtomsFromRotLib(true);

	// Each position build opt.numRotamer rotamers
	for (uint i = 0; i < _sys.positionSize();i++){

	  MSLOUT.stream() << "I: "<<i<<endl;
	  Residue &res = _sys.getResidue(i);
	  string id = res.getIdentityId();

	  AtomContainer ac = pdbTop.getResidue(id,res.getAtomPointers(),_opt.numRotamers);
	  _inverseRotamers.push_back(new AtomContainer(ac.getAtomPointers()));

	  AtomContainer *lastRes = _inverseRotamers.back();

	  
	  // Compute the bounding box of the C-alpha atoms...
	  map<string,double> box;
	  box["minX"] = MslTools::doubleMax;
	  box["minY"] = MslTools::doubleMax;
	  box["minZ"] = MslTools::doubleMax;
	  box["maxX"] = -MslTools::doubleMax;
	  box["maxY"] = -MslTools::doubleMax;
	  box["maxZ"] = -MslTools::doubleMax;

	  Atom &ca = (*lastRes).getAtom(id+",CA");
	  for (uint j = 0; j < (*lastRes)(0).getNumberOfAltConformations();j++){

	    ca.setActiveConformation(j);

	    if (ca.getX() < box["minX"]){
	      box["minX"] = ca.getX();
	    }

	    if (ca.getY() < box["minY"]){
	      box["minY"] = ca.getY();
	    }

	    if (ca.getZ() < box["minZ"]){
	      box["minZ"] = ca.getZ();
	    }

	    if (ca.getX() >= box["maxX"]){
	      box["maxX"] = ca.getX();
	    }

	    if (ca.getY() >= box["maxY"]){
	      box["maxY"] = ca.getY();
	    }

	    if (ca.getZ() >= box["maxZ"]){
	      box["maxZ"] = ca.getZ();
	    }
	   
	    
	  }
	  // Reset Ca conformation
	  ca.setActiveConformation(0);

	  // Add bounding box
	  _boundingBoxes.push_back(box);
	  MSLOUT.stream() << "BOX["<<_boundingBoxes.size()-1<<"]: "<<box["minX"]<<","<<box["minY"]<<","<<box["minZ"]<<" and "<<box["maxX"]<<","<<box["maxY"]<<","<<box["maxZ"]<<endl;

	  // DEBUG SECTION
	  if (_opt.debug){
	    AtomContainer *lastRes = _inverseRotamers.back();
	    for (uint j = 0; j < (*lastRes)(0).getNumberOfAltConformations();j++){

	      // Change each atoms active conformation
	      for (uint k =0;k <(*lastRes).size();k++){
		(*lastRes)(k).setActiveConformation(j);
	      }


	      // write it out.
	      char name[80];
	      sprintf(name,"/tmp/position-%04d_rotamer-%04d.pdb",i,j);
	      PDBWriter pout;
	      pout.open((string)name);
	      pout.write((*lastRes).getAtomPointers());
	      pout.close();
	    }

	    PyMolVisualization pyviz;

	    for (uint b = 0; b < _boundingBoxes.size();b++){
	      
	      // Make a box...
	      /*
                     ______________
                 ___|___________  |
                |   |         |   |
		1 - 2         5 - 6
                |   |         |   |
		4 - 3         8 - 7
		|___|_________|   |
		    |_____________|
	       */
	  CartesianPoint pt1(
				 _boundingBoxes[b]["minX"],
				 _boundingBoxes[b]["minY"],
				 _boundingBoxes[b]["minZ"]);

	  CartesianPoint pt2(
				 _boundingBoxes[b]["minX"],
				 _boundingBoxes[b]["minY"],
				 _boundingBoxes[b]["maxZ"]);

	  CartesianPoint pt3(
				 _boundingBoxes[b]["maxX"],
				 _boundingBoxes[b]["minY"],
				 _boundingBoxes[b]["maxZ"]);

	  CartesianPoint pt4(
				 _boundingBoxes[b]["maxX"],
				 _boundingBoxes[b]["minY"],
				 _boundingBoxes[b]["minZ"]);

	  CartesianPoint pt5(
				 _boundingBoxes[b]["minX"],
				 _boundingBoxes[b]["maxY"],
				 _boundingBoxes[b]["minZ"]);

	  CartesianPoint pt6(
				 _boundingBoxes[b]["minX"],
				 _boundingBoxes[b]["maxY"],
				 _boundingBoxes[b]["maxZ"]);

	  CartesianPoint pt7(
				 _boundingBoxes[b]["maxX"],
				 _boundingBoxes[b]["maxY"],
				 _boundingBoxes[b]["maxZ"]);

	  CartesianPoint pt8(
				 _boundingBoxes[b]["maxX"],
				 _boundingBoxes[b]["maxY"],
				 _boundingBoxes[b]["minZ"]);



	  pyviz.createCylinder(pt1,pt2,"random",0.2);
	  pyviz.createCylinder(pt2,pt3,"random",0.2);
	  pyviz.createCylinder(pt3,pt4,"random",0.2);
	  pyviz.createCylinder(pt4,pt1,"random",0.2);

	  pyviz.createCylinder(pt5,pt6,"random",0.2);
	  pyviz.createCylinder(pt6,pt7,"random",0.2);
	  pyviz.createCylinder(pt7,pt8,"random",0.2);
	  pyviz.createCylinder(pt8,pt5,"random",0.2);

	  pyviz.createCylinder(pt1,pt5,"random",0.2);
	  pyviz.createCylinder(pt2,pt6,"random",0.2);
	  pyviz.createCylinder(pt3,pt7,"random",0.2);
	  pyviz.createCylinder(pt4,pt8,"random",0.2);

	  
	    }

	  ofstream fout;
	  fout.open("/tmp/boxes.py");
	  fout << pyviz<<endl;
	  fout.close();
	    
	  }
	  // END DEBUG

	}

}

void computeInterBoxDistances(vector<double> _minDistances, vector<double> _maxDistances,vector<map<string,double> > &_boundingBoxes){

	// Compute distances between bounding boxes...
	vector<vector<pair<double,double> > > coarseDistRange;
	coarseDistRange.resize(_boundingBoxes.size());


	for (uint i = 0; i < _boundingBoxes.size();i++){
	  // Box center + diagnol length
	  CartesianPoint cent1(
			       (_boundingBoxes[i]["minX"] +_boundingBoxes[i]["maxX"]) / 2,
			       (_boundingBoxes[i]["minY"] +_boundingBoxes[i]["maxY"]) / 2,
			       (_boundingBoxes[i]["minZ"] +_boundingBoxes[i]["maxZ"]) / 2);
	  CartesianPoint corner1(
				 _boundingBoxes[i]["minX"],
				 _boundingBoxes[i]["minY"],
				 _boundingBoxes[i]["minZ"]);
	  
	  double boxDiag1 = cent1.distance(corner1);

	  for (uint j = i+1; j < _boundingBoxes.size();j++){

	    // For now average max,min points for distance calculation
	    CartesianPoint cent2(
				 (_boundingBoxes[j]["minX"] +_boundingBoxes[j]["maxX"]) / 2,
				 (_boundingBoxes[j]["minY"] +_boundingBoxes[j]["maxY"]) / 2,
				 (_boundingBoxes[j]["minZ"] +_boundingBoxes[j]["maxZ"]) / 2);

	    CartesianPoint corner2(
				   _boundingBoxes[j]["minX"],
				   _boundingBoxes[j]["minY"],
				   _boundingBoxes[j]["minZ"]);

	    double dist = cent1.distance(cent2);
	    double boxDiag2 = cent2.distance(corner2);

	    double maxDist = dist + boxDiag1 + boxDiag2;
	    double minDist = dist - boxDiag1 - boxDiag2;
	    if (minDist < 3.5) minDist = 3.5;
	  
	    coarseDistRange[i].push_back(pair<double,double>(minDist,maxDist));

	    MSLOUT.stream() << "I: "<<i<<" J: "<<j<<" range: "<<minDist<<","<<maxDist<<endl;

	    _minDistances.push_back(minDist);
	    _maxDistances.push_back(maxDist);

	  }
	}



}



string extractFileName(string _input){
  vector<string> toks = MslTools::tokenize(_input,":");
  return toks[0];
}
pair<string,string> extractPositions(string _input){
  //FORMAT: /Users/dwkulp/work/VaccineDesign_PGT128/tertFragSearch/pgt128_noGlycans.pdb:,A 102 PRO 0--X,102

  vector<string> toks = MslTools::tokenize(_input,"--");
  vector<string> toks2 = MslTools::tokenize(toks[0],",");
  
  vector<string> pos1 = MslTools::tokenize(toks2[1]," ");
  stringstream pos1ss;
  pos1ss << pos1[0]<<","<<pos1[1];


  return pair<string,string>(pos1ss.str(),toks[1]);
}

void filterRotamerPairs(){

	// All Rotamers from each position vs All Rotamers for All other positions
	//  1.  clashes between backbone of one rotamer and functional (Cb onward) atoms should result in removal of the rotamer (one attached to backbone atoms which clashed)
	//  2.  backbone-backbone clashes means this pair can not exist together , maybe an exclude-rotamer-pair list?


}



	/*
	  VectorData format.

	  Option1: vector<VectorData> sorted by Distance,then Dot Product
	  Option2: map<string, vector<VectorData> > DistanceBin:DotProductBin
	 */
	// Collect edges.  
	//   Find a position pair that has VectorData = to FunctionalRotamerVectorData
	//       Add edge....

	  /*
	         VectorDataHash object
		 
		   Store VectorData for each pair of residues in 2 formats:
		       1) map<string, vector<VectorData *> > geometricHash;     //  key is DistanceBin:DotProductBin
		       2) map<string, VectorData>  pairPositionHash;            //  key is PDB;PosId1:PosId2 , value VectorData
		   
		   Inverse rotamer geo hash must contain information such as PosId,RotamerId and index into rotamer list?

	           void searchForVectorMatch (map<string, vector<VectorData> > &invereRotamerGeoHash, string pairPosition, int uniqueVectorMatches);
		   void searchForAllVectorMatches(map<string, vector<VectorData> > &inverseRotamerGeoHash, int uniqueVectorMatches);

		   --> Foreach position1:
		            keep track of position1 hits. If hits >= uniquePositions, then store as a "match"

			    map<string, vector<pair<string,string> > > , key = "FuncPosX", value = "FuncPosY", "PosJ"
		        --> Foreach position2:
			       call searchForVectorMatch()
			       store VectorData match... in above format?
			       if matched inputPos1Pos2 is unique, increment counter.


			    if enough uniquePosition matches, then put into completeSystems.


		   Create :
		       map<string, DATA> completeSystems;
		   
		       DATA = vector<    vector<FunctionalPosition mapping to CandidatePosition> > , inner vector is uniqueVectorMatches long. outer vector is the number of complete matches.
	   */


	// Position1, Position2 = vector< pair<RotamerPosN, RotamerPosM> >;
	/*

	  If a position has edges to all functional positions:
	       Position1 -> [ FuncPos1 to FuncPos2 (Position3,Position9...) , FuncPos1 to FuncPos3 (), FuncPos1 to FuncPos4() ] = Can be considered a FuncPos1 = CandidateFunc1
	       Position2 -> [ FuncPos1 to FuncPos2 (Position3,Position9...) , FuncPos1 to FuncPos3 (), FuncPos1 to FuncPos4() ] = Can be considered a FuncPos1 = CandidateFunc1
	       ....
	       Position2 -> [ FuncPos2 to FuncPos1 (Position3,Position9...) , FuncPos2 to FuncPos3 (), FuncPos2 to FuncPos4(..) ] = Can be considered a FuncPos2 = CandidateFunc2

	       If a position does not have ALL proper edges, remove it.

	       Now we need a search for the combinations...... 

	          For all candidateFunc1:
		        Is the FuncPos2 of candiateFunc1 a candidateFunc2 ?
			Is the FuncPos3 of candiateFunc1 a candiateFunc3 ?
			Is the FuncPos4 of candiateFunc1 a candiateFunc4 ?

			If all true, then print 
			   CanidiateFunc1 = PosX with RotamerZ
			   CanidiateFunc2 = PosX with RotamerZ
			   CanidiateFunc3 = PosX with RotamerZ
			   CanidiateFunc4 = PosX with RotamerZ
		done.

	       
	 */

