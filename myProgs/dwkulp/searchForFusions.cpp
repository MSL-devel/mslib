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
#include "MslExceptions.h"


#include "searchForFusions.h"



using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("searchForFusions");


// Utility functions
void sighandler(int sig);
void fillTerminalFrames(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points);
void fillTerminalFramesAve(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points);
void fillLoopFrames(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points);
vector<double> getAngles(Frame &_frame, vector<CartesianPoint> &_points);
double getRadius(AtomPointerVector &_av);
double getSasa(System &_partner, System &_ref, pair<string,string> _closestChains, double _totalUnboundSasa);
int numberOfCalphaClashes(AtomPointerVector &_partner, AtomPointerVector  &_ref);
pair<string,string> getClosestChains(System &_partner, System &_ref);



// Fusion data 
struct fusionDat {
  AtomContainer ac;
  double complexSasa;
  string name;
  fusionDat(){
    complexSasa = 0.0;
    name = "";
  }

  ~fusionDat(){
    ac.removeAllAtoms();
  }

  fusionDat(const fusionDat &_rhs){
    fusionDat *_rhs_nonconst = const_cast<fusionDat *>(&_rhs);
    this->ac.removeAllAtoms();
    this->ac.addAtoms(_rhs_nonconst->ac.getAtomPointers());
    _rhs_nonconst = NULL;

    this->complexSasa = _rhs.complexSasa;
    this->name = _rhs.name;

  }

  const fusionDat& operator=(const fusionDat &_rhs){
    if (this == &_rhs) return *this;


    fusionDat *_rhs_nonconst = const_cast<fusionDat *>(&_rhs);
    this->ac.removeAllAtoms();
    this->ac.addAtoms(_rhs_nonconst->ac.getAtomPointers());
    _rhs_nonconst = NULL;


    this->complexSasa = _rhs.complexSasa;

    this->name = _rhs.name;


    return (*this);
  }

  
};

bool operator<( const fusionDat& p1, const fusionDat& p2 ) {  
      return p1.complexSasa < p2.complexSasa;
}


priority_queue<fusionDat>  bestFusions;



int main(int argc, char *argv[]) {
  
  // Register signal handlers
  signal(SIGABRT, &sighandler);
  signal(SIGTERM, &sighandler);
  signal(SIGINT, &sighandler);

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // MslOut can suppress output, to make example output clean
  //MSLOUT.turnOn("searchForFusions");

  // Read in the reference pdb file
  System sys;
  sys.readPdb(opt.ref);
  double refRadius = getRadius(sys.getChain(0).getAtomPointers());

  //  Get AtomPointerVector for all but last two positions of refChain
  AtomPointerVector refChainNoEnds;
  for (uint p = 2; p < sys.getChain(0).positionSize()-2;p++){
    refChainNoEnds += sys.getChain(0).getPosition(p).getAtomPointers();
  }

  // Store unbound sasa
  SasaCalculator sas(refChainNoEnds,1.4,100);
  sas.calcSasa();
  double refSasa   = sas.getTotalSasa();


  // Define a set of frames on the reference pdb
  MSLOUT.stream() << "Fill Terminal Frames for reference pdb"<<endl;
  map<string,Frame> refFrames;
  map<string,vector<CartesianPoint> > refCoors;//each frame has 3+ coordinates 
  //fillTerminalFrames(sys,refFrames,refCoors);
  fillTerminalFramesAve(sys,refFrames,refCoors);
  sys.getAtomPointers().saveCoor("orig");

  MSLOUT.stream() << "Read text file"<<endl;
  vector<string> lines;
  if (opt.list != ""){
    MslTools::readTextFile(lines, opt.list);
  }

  if (opt.testPdb != ""){
    lines.push_back(opt.testPdb);
  }


  // Transforms object helps with rmsd alignments
  Transforms move;

  // For each file..
  for (uint i = 0; i < lines.size();i++){
    string fileName = MslTools::trim(lines[i]);
    
    cout << "Reading : "<<fileName<<endl;

    System partner;
    if (!partner.readPdb(fileName)){
      cerr << "ERROR 34245 in reading PDB: "<<fileName<<endl;
      exit(34245);
    }
    partner.getAtomPointers().saveCoor("orig");
    MSLOUT.stream() << MslTools::getFileName(fileName)<<" has "<<partner.getSizes()<<endl;
    double partnerRadius = getRadius(partner.getChain(0).getAtomPointers());

    // Store unbound sasa
    SasaCalculator sas2(partner.getChain(0).getAtomPointers(),1.4,100);
    sas2.calcSasa();
    double totalUnboundSasa   = sas2.getTotalSasa() + refSasa;


    // Rename chains to X,Y,Z
    if (partner.chainSize() != 3) {
      cerr << "PDB must be a homotrimer ! this file: "<<fileName<<" has "<<partner.chainSize()<<" chains!"<<endl;
      continue;
    }

    partner.getChain(0).setChainId("X");
    partner.getChain(1).setChainId("Y");
    partner.getChain(2).setChainId("Z");

    map<string,Frame> partnerFrames;
    map<string,vector<CartesianPoint> > partnerCoors;
    //MSLOUT.stream() << "Fill Terminal Frames for "<<MslTools::getFileName(fileName)<<" pdb"<<endl;
    //fillTerminalFrames(partner,partnerFrames,partnerCoors);
    //fillTerminalFramesAve(partner,partnerFrames,partnerCoors);
    MSLOUT.stream() << "Fill Loop Frames for "<<MslTools::getFileName(fileName)<<" pdb"<<endl;
    fillLoopFrames(partner,partnerFrames,partnerCoors);


    // Store best fusions by SASA.
    //priority_queue<fusionDat>  bestFusions;
    bestFusions	  = priority_queue<fusionDat>();

    // For each reference frame
    map<string,Frame>::iterator refIt;
    MSLOUT.stream() << "Number of reference frames: "<<refFrames.size()<<endl;
    for (refIt = refFrames.begin();refIt != refFrames.end();refIt++){

      MSLOUT.stream() << "Working on frame: "<<refIt->first<<endl;

      // Make sure the reference system is in its original orientation + position.
      sys.getAtomPointers().applySavedCoor("orig");

      double refDistance = refIt->second.getCenter().distance(refCoors[refIt->first][0]);
      vector<double> refAngles = getAngles(refIt->second,refCoors[refIt->first]);
      double maxRefAngle = *(std::max_element(refAngles.begin(),refAngles.end()));
      double minRefAngle = *(std::min_element(refAngles.begin(),refAngles.end()));


      refIt->second.transformToGlobalBasis(sys.getAtomPointers());

      char name2[80];
      sprintf(name2, "ref_%s_%s.pdb",MslTools::getFileName(fileName).c_str(),refIt->first.c_str());
      //sys.writePdb((string)name2);
	    
      // For each partner frame
      map<string,Frame>::iterator partnerIt;
      MSLOUT.stream() << "Number of partner frames: "<<partnerFrames.size()<<endl;


      for (partnerIt = partnerFrames.begin();partnerIt != partnerFrames.end();partnerIt++){

	// Skip in compatible ends N-N or C-C linkages..
	if ( (refIt->first.substr(0,9) == "terminalN" && partnerIt->first.substr(0,9) == "terminalN") ||
	     (refIt->first.substr(0,9) == "terminalC" && partnerIt->first.substr(0,9) == "terminalC")) continue;


	// Skip in compatible ends : terminal average on ref shouldn't be paired with terminalN or terminalC, only internal fusions or "loop" frames
	if ( (refIt->first.substr(0,9) == "terminals" || partnerIt->first.substr(0,5) == "loops") ) {

	  if (refIt->first.substr(0,9) == "terminalN" || refIt->first.substr(0,9) == "terminalC"){
	    continue;
	  }
	  if (partnerIt->first.substr(0,9) == "terminalN" || partnerIt->first.substr(0,9) == "terminalC") {
	    continue;
	  }

	}


	MSLOUT.stream() << "\tChecking frame: "<<partnerIt->first<<endl;

	// check distances, angles
	double partnerDistance = partnerIt->second.getCenter().distance(partnerCoors[partnerIt->first][0]);

	if (abs(refDistance - partnerDistance) > 10) {
	  //MSLOUT.stream() << "\t\tDistance check failed for "<<MslTools::getFileName(fileName)<<" "<<refIt->first<<","<<partnerIt->first<<" "<<refDistance<<" vs "<<partnerDistance<<endl;
	  continue;
	}
	vector<double> partnerAngles = getAngles(partnerIt->second,partnerCoors[partnerIt->first]);
	double maxPartnerAngle = *(std::max_element(partnerAngles.begin(),partnerAngles.end()));
	double minPartnerAngle = *(std::min_element(partnerAngles.begin(),partnerAngles.end()));

	      
	if (abs(maxRefAngle - maxPartnerAngle) > 30 || abs(minRefAngle - minPartnerAngle) > 30){
	  //MSLOUT.stream() << "\t\tAngle check failed for "<<MslTools::getFileName(fileName)<<refIt->first<<","<<partnerIt->first<<endl;
	  continue;
	}


	// Make sure the partner system is in its original orientation + position.
	partner.getAtomPointers().applySavedCoor("orig");

	// align partner frame to global basis
	partnerIt->second.transformToGlobalBasis(partner.getAtomPointers());


	// Keep track of best fusion
	fusionDat bestFusion;
	bestFusion.complexSasa = -MslTools::doubleMax;
	bestFusion.name = "";


	partner.getAtomPointers().saveCoor("pre_invert");
	for (uint invert=0; invert < 2;invert++){
	  partner.getAtomPointers().applySavedCoor("pre_invert");

	  // Invert the partner (rotate x 180, then z 180);
	  if (invert == 1){
	    move.Xrotate(partner.getAtomPointers(),180);
	    move.Zrotate(partner.getAtomPointers(),180);
	    //partner.writePdb("fooInvert.pdb");
	  }

	  // Discover closest chains
	  pair<string,string> closestChains = getClosestChains(partner,sys);

	  partner.getAtomPointers().saveCoor("pre_translate");

	  // For a range of Z-values 	 (+/- 5 Angstroms)
	  CartesianPoint z_coord(0.0,0.0,1.0);

	  for (int z=-6; z <= 6;z+=2){
	    partner.getAtomPointers().applySavedCoor("pre_translate");
		  
	    CartesianPoint transVec = z_coord * z;
		
	    // Translate atoms..
	    move.translate(partner.getAtomPointers(),transVec);


	    // Rotate +/- around Z-axis...
	    partner.getAtomPointers().saveCoor("pre_zrot");
	    for (int zrot=-20; zrot <= 20;zrot+=10){
	      partner.getAtomPointers().applySavedCoor("pre_zrot");
	      
	      // Rotate atoms
	      move.Zrotate(partner.getAtomPointers(),zrot);


		  /* **** TEST POTENTIAL FUSION ****** */

		  // Do simple sphere collision test
		  double distGC = partner.getChain(closestChains.first).getAtomPointers().getGeometricCenter().distance(sys.getChain(closestChains.second).getAtomPointers().getGeometricCenter());
		  double touchingDistance = refRadius + partnerRadius;

		  if (distGC < touchingDistance-5){
		    //MSLOUT.stream() << "\t\tSpherical collision detected "<<distGC<<" , "<<touchingDistance-5<<" for "<<MslTools::getFileName(fileName)<<refIt->first<<","<<partnerIt->first<<endl;
		    continue;
		  }
	    
		  if (distGC > touchingDistance+5){
		    //MSLOUT.stream() << "\t\tProteins are too far away "<<distGC<<" , "<<touchingDistance+5<<" for "<<MslTools::getFileName(fileName)<<refIt->first<<","<<partnerIt->first<<endl;
		    continue;
		  }
		
		  // Do simple bb-atom clash check
		  int clashes = numberOfCalphaClashes(partner.getChain(closestChains.first).getAtomPointers(),sys.getChain(closestChains.second).getAtomPointers());
		  if (clashes >= opt.maxCaClashes){
		    //MSLOUT.stream() << "\t\tProteins have too many C-alpha clashes: "<<clashes<< " >= "<<opt.maxCaClashes<<endl;
		    continue;
		  }
	    
	    
		  // Get Sasa between chains A,B,C and X,Y,Z.
		  double deltaSasa = getSasa(partner,sys,closestChains,totalUnboundSasa); // remove ending positions of sys those shouldn't count in SASA calculation

		  if (deltaSasa > bestFusion.complexSasa){

		    bestFusion.complexSasa = deltaSasa;

		    bestFusion.ac.removeAllAtoms();
		    bestFusion.ac.addAtoms(partner.getAtomPointers());
		    bestFusion.ac.addAtoms(sys.getAtomPointers());

		    char name[80];
		    sprintf(name, "fusion_%s_%s_%s_%1d_%03d_%03d.pdb",MslTools::getFileName(fileName).c_str(),refIt->first.c_str(),partnerIt->first.c_str(),invert,z,zrot);
		    bestFusion.name = (string)name;

		    MSLOUT.stream() << "NEW BEST FUSION: "<<bestFusion.name<<" "<<bestFusion.complexSasa<<endl;
		  }

	    } // ZROT
	  }// Z trans

	  
	} // INVERT

	// Keep track of best fusion for this ref-partner frame pair.
	if (bestFusion.name != ""){

	  if (bestFusions.size() < 10){
	    bestFusions.push(bestFusion);
	  } else {
	    priority_queue<fusionDat>  tmpFusions;
	    while (bestFusions.size() != 1){
	      tmpFusions.push(bestFusions.top());
	      bestFusions.pop();
	    }
	    bestFusions = priority_queue<fusionDat>();
	    while (!tmpFusions.empty()){
	      bestFusions.push(tmpFusions.top());
	      tmpFusions.pop();
	    }
	    bestFusions.push(bestFusion);
	    tmpFusions = priority_queue<fusionDat>();
	  }
	}



      } // PARTNER FRAMES


    } // REF FRAMES


	  
    // Take fusions that bury the most surface area and dump to pdb.
    for (uint top=0; top < 10 && !bestFusions.empty(); top++){

      fusionDat item = bestFusions.top();
      bestFusions.pop();

      PDBWriter pout;
      pout.open(item.name);
      AtomPointerVector &ac = item.ac.getAtomPointers();
      for (uint a = 0; a < ac.size();a++){
	ac(a).setSegID("");
      }
      pout.write(ac);
      pout.close();
    
      MSLOUT.fprintf(stdout,"POTENTIAL FUSION(%s)[%03d]: %s %8.3f\n",MslTools::getFileName(fileName).c_str(),top,item.name.c_str(),item.complexSasa);
      item.ac.removeAllAtoms();
    }
  } // FOR LINES
	

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
    cout << "searchForFusions --ref PDB --list pdbs_to_search\n";
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


  opt.ref = OP.getString("ref");
  if (OP.fail()){
    cerr << "ERROR 1111 ref not specified.\n";
    exit(1111);
  }
  opt.list = OP.getString("list");
  opt.testPdb = OP.getString("testPdb");

  if (opt.list == "" && opt.testPdb == ""){
    cerr << "ERROR 1111 --list OR --testPdb must be used"<<endl;
    exit(1111);
  }

  if (opt.list != "" && opt.testPdb != ""){
    cerr << "ERROR 1111 --list OR --testPdb must be used"<<endl;
    exit(1111);
  }

  opt.maxCaClashes = OP.getDouble("maxCaClashes");
  if (OP.fail()){
    opt.maxCaClashes = 10;
    cerr << "WARNING 1111 maxCaClashes defaulted to "<<opt.maxCaClashes<<endl;
  }

  return opt;
}



/*
  Add a terminal frame for C-terminal and N-Terminal
*/
void fillTerminalFrames(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points){

  /* NTERM */

  // All combinations of +/- 10 residues from each terminal
  vector<CartesianPoint> points;
  for (uint c = 0; c < _sys.chainSize();c++){    

    CartesianPoint thisChainGC(0.0,0.0,0.0);

    for (uint nterm=0; nterm < 10; nterm++){
      Position &pos1 = _sys.getChain(c).getPosition(nterm);
      if (!pos1.atomExists("CA")) continue;

      thisChainGC += pos1.getAtom("CA").getCoor();
    }

    thisChainGC /= 10;

    points.push_back(thisChainGC);
  }
  

  if (points.size() != 3){
    cerr << "ERROR 3 points were not found"<<endl;
  } else {

    // Compute a frame please!
    Frame aframe;
    aframe.computeFrameFrom3Points(points[0],points[1],points[2],true);


    char name[80];
    sprintf(name, "terminalN_%s%d_%d",
	    _sys.getChain(0).getChainId().c_str(),
	    _sys.getChain(0).getPosition(0).getResidueNumber(),
	    _sys.getChain(0).getPosition(9).getResidueNumber());

    string fname = (string)name;
    _frames[fname] = aframe;
    _points[fname] = points;
  }

  /* CTERM */

  // All combinations of +/- 10 residues from each terminal
  points.clear();
  for (uint c = 0; c < _sys.chainSize();c++){    

    CartesianPoint thisChainGC(0.0,0.0,0.0);

    for (uint cterm=1; cterm <= 10; cterm++){
      Position &pos1 = _sys.getChain(c).getPosition(_sys.getChain(c).positionSize()-cterm);
      if (!pos1.atomExists("CA")) continue;

      thisChainGC += pos1.getAtom("CA").getCoor();
    }

    thisChainGC /= 10;

    points.push_back(thisChainGC);
  }
  

  if (points.size() != 3){
    cerr << "ERROR 3 points were not found"<<endl;
  } else {

    // Compute a frame please!
    Frame aframe;
    aframe.computeFrameFrom3Points(points[0],points[1],points[2],true);


    char name[80];
    sprintf(name, "terminalC_%s%d_%d",
	    _sys.getChain(0).getChainId().c_str(),
	    _sys.getChain(0).getPosition(_sys.getChain(0).positionSize()-10).getResidueNumber(),
	    _sys.getChain(0).getPosition(_sys.getChain(0).positionSize()-1).getResidueNumber());

    string fname = (string)name;
    _frames[fname] = aframe;
    _points[fname] = points;
  }


}
void fillTerminalFramesAve(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points){

  // All combinations of +/- 10 residues from each terminal

  for (uint nterm=0; nterm < 5; nterm++){


    for (uint cterm=1; cterm < 5; cterm++){

      // Collect the points (geometric center of coordinates of CA atoms from terminal residues)
      vector<CartesianPoint> points;
      for (uint c = 0; c < _sys.chainSize();c++){
	Position &pos1 = _sys.getChain(c).getPosition(nterm);

	if (!pos1.atomExists("CA")) continue;

	Position &pos2 = _sys.getChain(c).getPosition(_sys.getChain(c).positionSize()-cterm);
	if (!pos2.atomExists("CA")) continue;

	// Create geometric center.
	CartesianPoint gc(0.0,0.0,0.0);
	gc += pos1.getAtom("CA").getCoor();
	gc += pos2.getAtom("CA").getCoor();
	gc /= 2.0;



	points.push_back(gc);

      } // FOR chains

      if (points.size() != 3){
	cerr << "ERROR 3 points were not found"<<endl;
	continue;
      }

      // Compute a frame please!
      Frame aframe;
      aframe.computeFrameFrom3Points(points[0],points[1],points[2],true);



      
      char name[80];
      sprintf(name, "terminals_%s%d_%d",
	      _sys.getChain(0).getChainId().c_str(),
	      _sys.getChain(0).getPosition(nterm).getResidueNumber(),
	      _sys.getChain(0).getPosition(_sys.getChain(0).positionSize()-cterm).getResidueNumber());

      string fname = (string)name;
      _frames[fname] = aframe;
      _points[fname] = points;
      //MSLOUT.stream() << "TERMINAL FRAME: "<<fname<<endl;
      //MSLOUT.stream() << "\tPOINTS: "<<points[0].toString()<<" "<<points[1].toString()<<" "<<points[2].toString()<<endl;
      //MSLOUT.stream() << "\tGC: "<<aframe.getCenter().toString()<<endl;

    } // FOR cterm

  } // FOR nterm
    
  //MSLOUT.stream() << "Done fillTerminalLoops"<<endl;
}
void fillLoopFrames(System &_sys, map<string,Frame> &_frames,map<string,vector<CartesianPoint> > &_points){

  // Loops between two regular secondary structure

  // Look at chain A for SSE pattern we want..
  RegEx re;
  re.setStringType(RegEx::SegID); // SSE is kept in SegID.
  string regex = "[H,E]{3}([L,T,S]{4,10})[H,E]{3}";


  // Now do a sequence search...
  vector<pair<int,int> > matchingResidueIndices = re.getResidueRanges(_sys.getChain(0),regex);
 
 
  // Loop over each match.
  for (uint m = 0; m < matchingResidueIndices.size();m++){
 
    // .. do something cool with matched residues ...
 

    vector<CartesianPoint> points; 
    for (uint c = 0; c < _sys.chainSize();c++){    

      CartesianPoint gc;
      int loopLength = 0; 

      // +/- 3 for [H,E]{3} on either side. the '(' ')' are not being obeyed at the moment.
      for (uint r = matchingResidueIndices[m].first+3; r <= matchingResidueIndices[m].second-3;r++){
 
	// Get the position
	Position *pos1 = NULL;
	try {
	  pos1 = &_sys.getChain(0).getPosition(r);
	  if (!pos1->atomExists("CA")) continue;
	} catch (MslNotFoundException){
	  continue;
	}

	// Get equivalent position on the other chain
	if (c > 0){
	  char idStr[80];
	  sprintf(idStr,"%s,%d",
		  _sys.getChain(c).getChainId().c_str(),
		  pos1->getResidueNumber());

	  try{
	    pos1 = &_sys.getChain(c).getPosition((string)idStr);
	  } catch (MslNotFoundException){
	    continue;
	  }
	}
	//MSLOUT.stream() << "MATCH RESIDUE: "<<pos1->toString()<<" "<<pos1->getAtom("CA").getSegID()<<endl;
	gc += pos1->getAtom("CA").getCoor();
	loopLength++;

      }
      if (loopLength == 0) continue;
      gc /= (double)loopLength; 
      points.push_back(gc);

    } // FOR CHAIN

    if (points.size() != 3){
      cerr << "ERROR 3 points were not found"<<endl;
      continue;
    }
    // Compute a frame please!
    Frame aframe;
    aframe.computeFrameFrom3Points(points[0],points[1],points[2],true);
    char name[80];
    try{
      sprintf(name, "loops_%s%d_%d",
	    _sys.getChain(0).getPosition(matchingResidueIndices[m].first+3).getChainId().c_str(),
	    _sys.getChain(0).getPosition(matchingResidueIndices[m].first+3).getResidueNumber(),
	    _sys.getChain(0).getPosition(matchingResidueIndices[m].second-3).getResidueNumber());
    } catch (MslNotFoundException){
      continue;
    }

    string fname = (string)name;
    _frames[fname] = aframe;
    _points[fname] = points;

  } // FOR MATCHES

  //MSLOUT.stream() << "Done fillLoopFrames"<<endl;
}
vector<double> getAngles(Frame &_frame, vector<CartesianPoint> &_points){

  vector<double> angles;
  for (uint i = 0; i < _points.size();i++){
    for (uint j = i+1; j < _points.size();j++){
      angles.push_back(_points[i].angle(_frame.getCenter(),_points[j]));
    }
  }

  return angles;
}

double getRadius(AtomPointerVector &_av){
  double minX = MslTools::doubleMax;
  double maxX = -MslTools::doubleMax;

  double minY = MslTools::doubleMax;
  double maxY = -MslTools::doubleMax;

  double minZ = MslTools::doubleMax;
  double maxZ = -MslTools::doubleMax;


  for (uint i = 0 ; i < _av.size();i++){
    if (_av(i).getX() >= maxX) maxX = _av(i).getX();
    if (_av(i).getX() <= minX) minX = _av(i).getX();

    if (_av(i).getY() >= maxY) maxY = _av(i).getY();
    if (_av(i).getY() <= minY) minY = _av(i).getY();

    if (_av(i).getZ() >= maxZ) maxZ = _av(i).getZ();
    if (_av(i).getZ() <= minZ) minZ = _av(i).getZ();
  }


  double maxDelta = maxX-minX;
  if (maxDelta < maxY-minY) maxDelta = maxY-minY;
  if (maxDelta < maxZ-minZ) maxDelta = maxZ-minZ;

  double minDelta = maxX-minX;
  if (minDelta > maxY-minY) minDelta = maxY-minY;
  if (minDelta > maxZ-minZ) minDelta = maxZ-minZ;

  return minDelta/2;

}


pair<string,string> getClosestChains(System &_partner, System &_ref){

  double minDistSq = MslTools::doubleMax;
  pair<string,string> result;
  for (uint c1 = 0; c1 < _partner.chainSize();c1++){
    for (uint c2 = 0; c2 < _ref.chainSize();c2++){

      double distSq = _ref.getChain(c2).getAtomPointers().getGeometricCenter().distance2(_partner.getChain(c1).getAtomPointers().getGeometricCenter());
      if ( distSq < minDistSq){
	minDistSq = distSq;
	result.first  = _partner.getChain(c1).getChainId();
	result.second = _ref.getChain(c2).getChainId();
      }
      
    }
  }
  
  return result;
}

int numberOfCalphaClashes(AtomPointerVector &_partner, AtomPointerVector  &_ref){

  int clashes = 0;
  for (uint a = 0 ; a  < _partner.size();a++){
    if (_partner(a).getName() != "CA") continue;

    for (uint b = 0 ; b  < _ref.size();b++){
      if (_ref(b).getName() != "CA") continue;


      double distSq = _partner(a).distance2(_ref(b));

      // < 3.0 Angstroms is a clash
      if (distSq < 9.00){
	clashes++;
      }
    }
  }

  return clashes;
}

double getSasa(System &_partner, System &_ref, pair<string,string> _closestChains, double _totalUnboundSasa){

  Chain &partnerChain = _partner.getChain(_closestChains.first);
  Chain &refChain = _ref.getChain(_closestChains.second);


  // Make more expensive, but more accurate..
  // Parse id , remove internal loop residues or terminal residues...
  

  //  Get AtomPointerVector for all but last two positions of refChain
  AtomPointerVector refChainNoEnds;
  for (uint p = 2; p < refChain.positionSize()-2;p++){
    refChainNoEnds += refChain.getPosition(p).getAtomPointers();
  }

  AtomPointerVector complex;
  complex += refChainNoEnds + partnerChain.getAtomPointers();
  
  SasaCalculator complexSasa(complex,1.4,100);
  complexSasa.calcSasa();

  return (_totalUnboundSasa - complexSasa.getTotalSasa());

}

void sighandler(int sig){

     cerr << "ERROR  ****  SIGNAL "<<sig<<" caught, try to output best fusions so far:"<<endl;

      // Take fusions that bury the most surface area and dump to pdb.
    for (uint top=0; top < 10 && !bestFusions.empty(); top++){

      fusionDat item = bestFusions.top();
      bestFusions.pop();

      PDBWriter pout;
      pout.open(item.name);
      AtomPointerVector &ac = item.ac.getAtomPointers();
      for (uint i = 0; i < ac.size();i++){
	ac(i).setSegID("");
      }
      pout.write(ac);
      pout.close();
    
      MSLOUT.fprintf(stdout,"POTENTIAL FUSION[%03d]: %s %8.3f\n",top,item.name.c_str(),item.complexSasa);
      item.ac.removeAllAtoms();
    }
}
