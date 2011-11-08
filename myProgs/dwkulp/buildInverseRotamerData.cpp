#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "DistanceHashing.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "discoverMotif.h"

using namespace std;
using namespace MSL;

// Build Inverse Rotamers from input PDB and collect the bounding boxes
void buildInverseRotamers(Options &_opt, System &_sys,vector<AtomContainer *> &_inverseRotamers, vector<map<string,double> > &_boundingBoxes);


// MslOut 
static MslOut MSLOUT("buildInverseRotamers");

int main(int argc, char *argv[]) {

        // MslOut can suppress output, to make example output clean
        //MSLOUT.turnAllOff();

	Options opt = setupOptions(argc,argv);

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
	for (uint i = 0; i < sys.positionSize();i++){

	  cout << "I: "<<i<<endl;
	  Residue &res = sys.getResidue(i);
	  string id = res.getIdentityId();

	  AtomContainer ac = pdbTop.getResidue(id,res.getAtomPointers(),opt.numRotamers);
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
	  cout << "BOX["<<_boundingBoxes.size()-1<<"]: "<<box["minX"]<<","<<box["minY"]<<","<<box["minZ"]<<" and "<<box["maxX"]<<","<<box["maxY"]<<","<<box["maxZ"]<<endl;

	  // DEBUG SECTION
	  if (opt.debug){
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

