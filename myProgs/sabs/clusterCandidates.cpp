#include "OptionParser.h"
#include "release.h"
#include "MslTools.h"
#include "Transforms.h"
#include "System.h"
#include "AtomSelection.h"

using namespace std;
using namespace MSL;

string programName = "clusterCandidates";
string programDescription = "This program creates gly helix dimers acoording to parameters in a file and clusters them based on rmsd";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "19 February 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


System sys;
System helicalSys;
map<string,string> model2GeomMap;

struct Options {
	string parFile;
	double rmsd;
	string pdbFile;
	string helicalAxisFile;
	string clusterFile;
	string centroidFile;

	/***** MANAGEMENT VARIABLES ******/
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> required; //list of required options
	vector<string> allowed; //list of allowed options
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string configfile;
};

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void usage() {
	cout << "Usage: clusterCandidates --configfile <conf>" << endl << endl;
	cout << "parfile <file> (contains lines of format \"00001 00 -10 -0.75 7.6 A 35 35 HA2 7.7 A 39 39 HA2 7.7 B 35 35 HA2 7.7 B 39 39 HA2 7.7\")" << endl << endl;
	cout << "pdbfile <templatefile> " << endl << endl;
	cout << "helicalaxisfile <file> " << endl << endl;
	cout << "rmsd <double> " << endl << endl;
	cout << "outcentroidfile <filename> (output file contains lines in same format as parfile - one for each centroid)\n" << endl;
	cout << "outclusterfile <filename> (one line for each cluster containing the ids of all models in the cluster  with the id of centroid at the start) \n" << endl;
	cout << endl;
}

Options parseOptions(int theArgc, char * theArgv[], Options defaults) {
	Options opt;
	vector<string> required;

	opt.required.push_back("parfile");
	opt.required.push_back("pdbfile");
	opt.required.push_back("helicalaxisfile");
	opt.required.push_back("rmsd");
	opt.required.push_back("outclusterfile");
	opt.required.push_back("outcentroidfile");

	opt.allowed.push_back("help");
	opt.allowed.push_back("h");
	opt.allowed.push_back("version");
	opt.allowed.push_back("v");
	opt.allowed.push_back("configfile");


	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	// Parse the options
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.setAllowed(opt.allowed);
	OP.setRequired(opt.required);	
	OP.autoExtendOptions();

	OP.readArgv(theArgc, theArgv);

	opt.errorMessages = "";
	opt.warningMessages = "";


	if (OP.countOptions() == 0){
		usage();
		exit(0);
	}

	opt.version = OP.getBool("version");
	if (opt.version) {
		version();
		exit(0);
	}

	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.help = OP.getBool("help");
	if (opt.help) {
		usage();
		exit(0);
	}


	opt.configfile = OP.getString("configfile");
	
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
		}
	}


	opt.parFile = OP.getString("parfile");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "parfile not specified \n";
	}
	opt.pdbFile = OP.getString("pdbfile");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "pdbfile not specified \n";
	}
	opt.helicalAxisFile = OP.getString("helicalaxisfile");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "helicalaxisfile not specified \n";
	}
	opt.rmsd = OP.getDouble("rmsd");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "rmsd not specified \n";
	}
	opt.centroidFile = OP.getString("outcentroidfile");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "outcentroidfile not specified \n";
	}
	opt.clusterFile = OP.getString("outclusterfile");
	if(OP.fail()) {
		opt.errorFlag = true;
		opt.errorMessages += "outclusterfile not specified \n";
	}

	return opt;
}

class Structure {

	public:
	Structure(string _line, System& _sys, System& _helicalSys) {
		parLine = _line;
		vector<string> toks = MslTools::tokenizeAndTrim(_line);	
		modelNum = toks[0];
		double axialR = MslTools::toDouble(toks[1]);
		double crossA = MslTools::toDouble(toks[2]);
		double zShift = MslTools::toDouble(toks[3]);
		double xShift = MslTools::toDouble(toks[4]);
		build(axialR,crossA,zShift,xShift,_sys,_helicalSys);

	}

	~Structure() {
		for(AtomPointerVector::iterator it = atoms.begin(); it != atoms.end(); it++) {
			delete *it;
		}
		atoms.clear();
	}
	double rmsd(Structure& _struct) {
		return atoms.rmsd(_struct.atoms);
	}

	void build(double _axialRotate, double _crossingAngle, double _zShift, double _xShift, System & _sys, System & _helicalSys ) {
		AtomSelection sel1(_sys.getAtomPointers());
		AtomPointerVector& chainA = sel1.select("A,chain A");
		AtomPointerVector& chainB = sel1.select("B,chain B");

		//====== Z Shift (Crossing Point) ======
		AtomSelection sel2(_helicalSys.getAtomPointers());
		AtomPointerVector& axisA = sel2.select("chain A");
		AtomPointerVector& axisB = sel2.select("chain B");

		
		CartesianPoint origin(0,0,0);
		CartesianPoint zAxis(0,0,1); 
		CartesianPoint xAxis(1,0,0); 

		//====== Z Shift (Crossing Point) ======
		CartesianPoint zShiftCP(0.0, 0.0, _zShift);

		Transforms trans;
		trans.translate(chainA, zShiftCP);

		//===== Axial Rotation ======
		trans.rotate(chainA, _axialRotate, origin, zAxis);

		//====== Local Crossing Angle ======
		trans.rotate(chainA, (_crossingAngle/2.0), origin, xAxis);
		trans.rotate(axisA,  (_crossingAngle/2.0), origin, xAxis);

		//====== X shift (Interhelical Distance) =======
		CartesianPoint interDistVect;
		interDistVect.setCoor((_xShift/-2.0), 0.0, 0.0);
		trans.translate(chainA, interDistVect);
		trans.translate(axisA, interDistVect);

		c2Symmetry(chainA, chainB);
		c2Symmetry(axisA, axisB);

		for(int i = 0; i < chainA.size(); i++) {
			atoms.push_back(new Atom(*(chainA[i])));
		}
		for(int i = 0; i < chainB.size(); i++) {
			atoms.push_back(new Atom(*(chainB[i])));
		}

	}

	void writePdb(string filename) {
		AtomContainer ac(atoms);
		if(!ac.writePdb(filename)) {
			cerr << "Unable to write " << filename << endl;
		}
	}

	void c2Symmetry(AtomPointerVector& _apvA,AtomPointerVector& _apvB) {
		for(int i = 0; i < _apvA.size(); i++) {
			_apvB[i]->setCoor(_apvA[i]->getCoor());
		}	
		// Rotation matrix for 180 degrees
		Matrix m(3,3,0.0);
		m[0][0] = -1.0;
		m[0][1] = 0.0;
		m[0][2] = 0.0;
		m[1][0] = 0.0;
		m[1][1] = -1.0;
		m[1][2] = 0.0;
		m[2][0] = 0.0;
		m[2][1] = 0.0;
		m[2][2] = 1.0;

		// Rotate chain B around Z axis
		Transforms trans; 
		trans.rotate(_apvB, m);
	}

	void print(ofstream& out) {
		out << parLine << endl;
	}

	AtomPointerVector atoms;
	string parLine;
	string modelNum;

};

Structure* getStructure(string _modelNum) {
	Structure* s = NULL;
	if(model2GeomMap.find(_modelNum) != model2GeomMap.end()) {
		 s = new Structure(model2GeomMap[_modelNum],sys,helicalSys);
	} else {
		cerr << "Unable to find params for modelNum " <<  _modelNum << endl;
		exit(0);
	}
	sys.applySavedCoor("origState");
	helicalSys.applySavedCoor("origState");
	return s;
}

class Cluster {
	public:
	void add(Structure * _s) {
		if(members.size() == 0) {
			for(int i = 0; i < _s->atoms.size(); i++) {
				centre.push_back(new Atom(*(_s->atoms[i])));
			}
		} else {
			for(int i = 0; i < _s->atoms.size(); i++) {
				CartesianPoint old = centre[i]->getCoor();
				centre[i]->setCoor( (old * members.size()+ _s->atoms[i]->getCoor()) / (members.size() + 1.0));
			}
		}
		members.push_back(_s->modelNum);

	}

	double rmsd(Structure* _s) {
		return centre.rmsd(_s->atoms);
	}

	int findCentroid() {
		int best = 0;
		Structure* s = getStructure(members[0]) ;
		double bestRmsd = rmsd(s);
		delete s;
		for(int i = 1; i < members.size(); i++) {
			s = getStructure(members[i]);
			double thisRmsd = rmsd(s);
			if(thisRmsd < bestRmsd) {
				bestRmsd = thisRmsd;
				best = i;
			}
			delete s;
		}	
		centroid = best;
		return centroid;
	}

	
	void printCentroid(ofstream& out) {

		if(model2GeomMap.find(members[centroid]) != model2GeomMap.end()) {
			out << model2GeomMap[members[centroid]] << endl;
		}
	}

	void printCluster(ofstream& out) {
		out << members[centroid];
		for(int i = 0; i < members.size(); i++) {
			if(i != centroid) {
				out <<  " " << members[i];
			}
		}
		out << endl;
	}

	vector<string> members;
	AtomPointerVector centre;
	int centroid;
};

int main(int argc, char* argv[]) {
	Options defaults;
	Options opt = parseOptions(argc,argv,defaults);
	if(opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;

		usage();
		exit(1);
	}	

	if(opt.warningFlag) {
		cerr << opt.warningMessages << endl;
	}

	// open the par file and the output files and check permissions
	ifstream par;
	par.open(opt.parFile.c_str());
	if(!par.is_open()) {
		cerr << "Unable to read " << opt.parFile << endl;
		exit(0);
	}

	ofstream centroidFile;
	centroidFile.open(opt.centroidFile.c_str());
	if(!centroidFile.is_open()) {
		cerr << "Unable to open " << opt.centroidFile << endl;
		exit(0);
	}
	ofstream clusterFile;
	clusterFile.open(opt.clusterFile.c_str());
	if(!clusterFile.is_open()) {
		cerr << "Unable to open " << opt.clusterFile << endl;
		exit(0);
	}

	if(!sys.readPdb(opt.pdbFile)) {
		cerr << "Unable to read " << opt.pdbFile << endl;
		exit(0);
	}
	sys.saveCoor("origState");
	
	if(!helicalSys.readPdb(opt.helicalAxisFile)) {
		cerr << "Unable to read " << opt.helicalAxisFile << endl;
		exit(0);
	}
	helicalSys.saveCoor("origState");

	vector<string> structures;
	string line;
	while(par.good()) {
		getline(par,line);
		cout << line << endl;
		if(line.substr(0,1) == "#" || line.length() < 8) {
			continue;
		}
		vector<string> toks = MslTools::tokenizeAndTrim(line);
		structures.push_back(toks[0]);
		model2GeomMap[toks[0]] = line;
	}

	vector<Cluster*> clusters;
	clusters.push_back(new Cluster);
	Structure* s = getStructure(structures[0]);
	clusters[0]->add(s);
	delete s;
	for(int i = 0; i < structures.size(); i++) {
		bool newCluster = true;
		cout << "Processing " << i << "th structure. Found " << clusters.size() << " clusters so far. " << endl;
		s = getStructure(structures[i]);
		for(int j = 0; j < clusters.size(); j++) {
			if(clusters[j]->rmsd(s) < opt.rmsd) {
				clusters[j]->add(s);
				newCluster = false;
				break;
			}
		}
		if(newCluster) {
			clusters.push_back(new Cluster);
			clusters.back()->add(s);
		}
		delete s;
	}

	cout << "Found " << clusters.size() << " clusters" << endl;
	for(int i = 0; i < clusters.size(); i++) {
		clusters[i]->findCentroid();
		clusters[i]->printCluster(clusterFile);
		clusters[i]->printCentroid(centroidFile);
	}

}
