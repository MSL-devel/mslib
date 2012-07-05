
// MSL Includes
#include "System.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "ResidueSelection.h"
#include "CartesianGeometry.h"
#include "CartesianPoint.h"
#include "AtomSelection.h"
#include "File.h"
#include "moveProteinCenter.h"

// STL Includes
#include<iostream>
#include<string>
using namespace std;
using namespace MSL;

int main(int argc, char *argv[]){
	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	// Read PDB
	System sys;
	sys.readPdb(opt.pdb);
	if (!sys.getPDBReader()->doesFileExist()) {
		cerr << "Pdb file does not exist, failing.\n";
		exit(1111);
	}
	
	AtomPointerVector centerAtoms;
	if (opt.sele != "") {
	    AtomSelection as(sys.getAtomPointers());
	    char selStr[80];
	    sprintf(selStr,"%s",opt.sele.c_str());
	    centerAtoms = as.select(selStr);
	}
	else {
	    centerAtoms = sys.getAtomPointers();
	}
	AtomPointerVector allAtoms = sys.getAtomPointers();
	
	CartesianPoint total(0.,0.,0.);
	for (uint i = 0; i < centerAtoms.size(); i++) {
	    total += centerAtoms[i]->getCoor();
	}
	CartesianPoint avg = total/centerAtoms.size();
	cout << "Average: " << avg.getX() << "\t" << avg.getY() << "\t" << avg.getZ() << endl;
	
	CartesianPoint newCenter(opt.x, opt.y, opt.z);
	for (uint i = 0; i < allAtoms.size(); i++) {
	    allAtoms[i]->getCoor() -= avg;
	    allAtoms[i]->getCoor() += newCenter;
	}
	
	PDBWriter writer;
	writer.open(opt.outfile);
	if (!writer.write(allAtoms)) {
	   cerr << "Problem writing " << opt.outfile << endl;
	}
	writer.close();
	
	exit(0); 
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0) {
		cout << "Usage:" << endl;
		cout << endl;
		cout << "moveProteinCenter --pdb PDB [--x X --y Y --z Z --sele SELE --outfile OUTFILE]\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "Pdb not specified, failing.\n";
		exit(1111);
	}

	opt.x = OP.getDouble("x");
	if (OP.fail()){
		cerr << "Warning, x coordinate set to 0.0\n";
		opt.x = 0.0;
	}

	opt.y = OP.getDouble("y");
	if (OP.fail()){
		cerr << "Warning, y coordinate set to 0.0\n";
		opt.y = 0.0;
	}

	opt.z = OP.getDouble("z");
	if (OP.fail()){
		cerr << "Warning, z coordinate set to 0.0\n";
		opt.z = 0.0;
	}

	opt.sele = OP.getString("sele");
	if (OP.fail()){
		opt.sele = "";
	}

	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		cerr << "Warning, outfile set to test.pdb" << endl;
		opt.outfile = "test.pdb";
	}

	return opt;
}

