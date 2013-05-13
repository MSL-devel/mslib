#include <iostream>
#include <fstream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "OptionParser.h"
#include "SelfPairManager.h"
#include "MslTools.h"
#include "DeadEndElimination.h"
#include "MonteCarloManager.h"
#include "SelfConsistentMeanField.h"
#include "HydrogenBondBuilder.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "MonteCarloManager.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "FormatConverter.h"
#include "CRDReader.h"
#include "CRDWriter.h"
#include "SysEnv.h"


using namespace MSL;
using namespace std;

string programName = "CAHTM";
string programDescription = "This program repacks two helices from a set starting position";
string programAuthor = "Benjamin K. Mueller, Sabareesh Subramaniam";
string programVersion = "0.0.10";
string programDate = "13 May 2013";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

static SysEnv ENV;

/******************************************************************************************************************************************************************************/


class HelixDimer : public AtomContainer {
	public:
	HelixDimer(string _name, double _energy, int _thread);

	void print();
	void setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds);

	void setDeltaEnergyByTerm(map<string,double> _deltaEbyTerm) {deltaEByTerm = _deltaEbyTerm;}

	map<string,double>& getDeltaEnergyByTerm(){return deltaEByTerm;}
	void addAxisAtom(Atom* _atom){axisAtoms.push_back(_atom);}
	AtomPointerVector& getAxes(){return axisAtoms;}

	bool operator< (HelixDimer& _b) const {
		if(energy < _b.energy) {
			return true;
		} else {
			return false;
		}
	}

	void printHelixDimerDetails(ofstream & _fout);
	int getNumHbonds() {return hbonds.size();}
	int getThread() {return thread;}
	
	string getName() {return name;}
	double getEnergy() {return energy;}
	double getXShift() {return xshift;}
	double getZShift() {return zshift;}
	double getAxialRotation() {return axialrot;}
	double getCrossingAngle() {return crossang;}
	string getInterface() {return interface;}
	string getProlineMask() {return prolineMask;}

	// Format: A,11,CA:B,12,O=2.6;
	vector<string>& getHbonds() {return hbonds;}
	

	private:

	string name;
	double energy;
	double xshift;
	double zshift;
	double axialrot;
	double crossang;
	string interface;
	int thread;
	string prolineMask;
	vector<string> hbonds;
	map<string,double> deltaEByTerm;
	AtomPointerVector axisAtoms; 
};

struct compareHelixDimers {
	bool operator () (HelixDimer* lhs, HelixDimer* rhs) {return *lhs < *rhs;}
};


HelixDimer::HelixDimer(string _name,double _energy, int _thread) : AtomContainer() {
	name = _name;
	energy = _energy;
	thread = _thread;
}

void HelixDimer::print() {
	cout << name << " " << energy << " " << hbonds.size() << " "  << thread << endl;
}

void HelixDimer::setHelixDimerDetails(double _x, double _z, double _a, double _c, string _interface, string _prolineMask, vector<string> _hbonds) {
	xshift = _x;
	zshift = _z;
	axialrot = _a;
	crossang = _c;
	interface = _interface;
	prolineMask = _prolineMask;
	hbonds = _hbonds;
}

void HelixDimer::printHelixDimerDetails(ofstream & _fout) {
	_fout << "interfacial " << interface << endl;
	_fout << "prolineMask " << prolineMask << endl;
	_fout << "name " << name << endl;
	_fout << "thread " << thread << endl;
	_fout << "Xshift " << xshift << endl;
	_fout << "Zshift " << zshift << endl;
	_fout << "axialRot " << axialrot << endl;
	_fout << "crossAng " << crossang << endl;
	_fout << "energy " << energy << endl;
	string hbondsList = "";
	for(int h = 0; h < hbonds.size(); h++) {
		if(h == 0) {
			hbondsList += hbonds[h];
		} else {
			hbondsList += ";" + hbonds[h];
		}
	}
	_fout << "hBonds " << hbondsList << endl;
}


class HelixDimerCluster {
	public:
	HelixDimerCluster(HelixDimer * _structure);

	void addHelixDimer(HelixDimer * _structure);

	AtomPointerVector& getAtomPointers();
	
	void setDetails(string _origSeq, string _seq, string _uName, string _uAccession,string _outputDir,int _resStart, int _resEnd, int _tmStart, int _tmEnd);

	void convertToPdbNames();
	void printHelixDimerClusterPdbs(int id); 
	void printHelixDimerClusterCrds(int id, bool allStructures, bool writeAxis); 

	void printDetails(int id);

	void makePse(int id); 

	vector<HelixDimer*>& getMembers() {return members;}
	
	string getOrigSequence() {return origSeq;}
	string getModeledSequence() {return seq;}

	string getUniprotName() {return uniprotName;}
	string getUniprotAccession(){return uniprotAccession;}

	int getResStart() {return resStart;}
	int getResEnd() {return resEnd;}

	int getTmStart() {return tmStart;}
	int getTmEnd() {return tmEnd;}

	void printTermEnergies(int id);

	private:
	vector<HelixDimer*> members; // should be added in sorted order
	string origSeq;
	string seq;
	string uniprotName;
	string uniprotAccession;
	string outputDir;
	int resStart;
	int resEnd;
	int tmStart;
	int tmEnd;

};


HelixDimerCluster::HelixDimerCluster(HelixDimer * _structure) {
	members.push_back(_structure);
}

void HelixDimerCluster::addHelixDimer(HelixDimer * _structure) {
	if(members.size() < 1) {
		members.push_back(_structure);
		return;
	}
	if(_structure->getEnergy() < 0) {
		members.push_back(_structure);
	}
}
AtomPointerVector& HelixDimerCluster::getAtomPointers() {
	if(members.size() > 0) {
		return members[0]->getAtomPointers();
	} else {
		cerr << "ERROR 13543: No members in cluster" << endl;
		exit(0);
	}
}
void HelixDimerCluster::setDetails(string _origSeq, string _seq, string _uName, string _uAccession, string _outputDir, int _resStart, int _resEnd, int _tmStart, int _tmEnd) {
	uniprotName = _uName;
	uniprotAccession = _uAccession;
	seq = _seq;
	origSeq = _origSeq;
	outputDir = _outputDir;
	resStart = _resStart;
	resEnd = _resEnd;
	tmStart = _tmStart;
	tmEnd = _tmEnd;
}

void HelixDimerCluster::printTermEnergies(int id) {
	char name[1000];
	sprintf(name,"%s/%s_%02d.energy",outputDir.c_str(),uniprotAccession.c_str(),id);
	ofstream eOut;
	eOut.open(name);
	map<string,double>& deltaEByTerm = members[0]->getDeltaEnergyByTerm();
	for(map<string,double>::iterator it = deltaEByTerm.begin(); it != deltaEByTerm.end(); it++) {
		eOut << it->first << " " << it->second << endl;
	}
	eOut.close();
}

void HelixDimerCluster::printHelixDimerClusterCrds(int id, bool allStructures,bool writeAxis) {
	//members[0]->print();
	CRDWriter writer;
	char name[1000];
	sprintf(name,"%s/%s_%02d.crd",outputDir.c_str(),uniprotAccession.c_str(),id);
	writer.open(string(name));
	writer.addRemark("Energy " + MslTools::doubleToString(members[0]->getEnergy()) );
	writer.addRemark("Name - " + members[0]->getName() + " HydrogenBonds " + MslTools::intToString(members[0]->getNumHbonds()) );
	AtomPointerVector& ats = members[0]->getAtomPointers();
	writer.write(ats);
	writer.close();

	if(writeAxis) {
		sprintf(name,"%s/%s_%02d_axis.crd",outputDir.c_str(),uniprotAccession.c_str(),id);
		writer.open(string(name));
		writer.write(members[0]->getAxes());
		writer.close();
	}

	if(allStructures) {
		for(int i = 0; i < members.size(); i++) {
			sprintf(name,"%s/%s_%02d_%03d.crd",outputDir.c_str(),uniprotAccession.c_str(),id,i);
			writer.open(string(name));
			writer.clearRemarks();
			writer.addRemark("Energy " + MslTools::doubleToString(members[i]->getEnergy()) );
			writer.addRemark("Name - " + members[i]->getName() + " HydrogenBonds " + MslTools::intToString(members[i]->getNumHbonds()) );
			AtomPointerVector& ats = members[i]->getAtomPointers();
			writer.write(ats);
			writer.close();
		}
	}
}

void HelixDimerCluster::convertToPdbNames() {
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");

	for(int i  = 0; i < members.size(); i++) {
		AtomPointerVector& ats = members[i]->getAtomPointers();
		fc.convert(ats);
	}
}

void HelixDimerCluster::printHelixDimerClusterPdbs(int id) {
	//members[0]->print();
	PDBWriter writer;
	writer.setConvertFormat("CHARMM22","PDB2");
	char name[1000];
	sprintf(name,"%s/%s_%02d.pdb",outputDir.c_str(),uniprotAccession.c_str(),id);
	writer.open(string(name));
	char helixA[1000];
	char helixB[1000];

	
	// Format is here http://www.wwpdb.org/documentation/format33/sect5.html
	//COLUMNS        DATA  TYPE     FIELD         DEFINITION
	//-----------------------------------------------------------------------------------
	// 1 -  6        Record name    "HELIX "
	// 8 - 10        Integer        serNum        Serial number of the helix. This starts
	//                                            at 1  and increases incrementally.
	//12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
	//                                            number, each helix is given an 
	//                                            alphanumeric character helix identifier.
	//16 - 18        Residue name   initResName   Name of the initial residue.
	//20             Character      initChainID   Chain identifier for the chain containing
	//                                            this  helix.
	//22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
	//26             AChar          initICode     Insertion code of the initial residue.
	//28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
	//32             Character      endChainID    Chain identifier for the chain containing
	//                                            this  helix.
	//34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
	//38             AChar          endICode      Insertion code of the terminal residue.
	//39 - 40        Integer        helixClass    Helix class (see below).
	//41 - 70        String         comment       Comment about this helix.
	//72 - 76        Integer        length        Length of this helix.
	//HELIX    1  HA ALA A   19  ILE A   41  1                                  23
	//HELIX    2  HA ALA B   19  ILE B   41  1                                  23


	
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");

	for(int i = 0; i < members.size(); i++) {
		AtomPointerVector& ats = members[i]->getAtomPointers();

		if(i == 0) {
			// BEWARE: hacky code - will work only for Homo
			Atom& startAtom = *ats[0];
			Atom& endAtom = *ats[ats.size() -1];

			int helixSize = resEnd - resStart + 1;

			sprintf(helixA,"HELIX    1  HA %3s A %4d  %3s A %4d  1                              %5d\n",fc.getResidueName(startAtom.getResidueName()).c_str(),startAtom.getResidueNumber(),fc.getResidueName(endAtom.getResidueName()).c_str(),endAtom.getResidueNumber(),helixSize);
			sprintf(helixB,"HELIX    2  HB %3s B %4d  %3s B %4d  1                              %5d\n",fc.getResidueName(startAtom.getResidueName()).c_str(),startAtom.getResidueNumber(),fc.getResidueName(endAtom.getResidueName()).c_str(),endAtom.getResidueNumber(),helixSize);

			string line = string(helixA) + string(helixB);
			writer.writeln(line); 

		}
		writer.clearRemarks();
		writer.addRemark("Energy " + MslTools::doubleToString(members[i]->getEnergy()) );
		writer.addRemark("Name - " + members[i]->getName() + " HydrogenBonds " + MslTools::intToString(members[i]->getNumHbonds()) );

		writer.writeREMARKS();
		writer.write(ats,true,false,true);
	}
	writer.close();
}

void HelixDimerCluster::printDetails(int id) {
	ofstream fout;
	char filename[1000];
	sprintf(filename,"%s/%s_%02d.txt",outputDir.c_str(),uniprotAccession.c_str(),id);
	fout.open(filename);
	fout << "originalSeq " << origSeq << endl;
	fout << "modelledSeq " << seq << endl;
	fout << "interfacial " << members[0]->getInterface() << endl;
	fout << "prolineMask " << members[0]->getProlineMask() << endl;
	fout << "resStart " << resStart << endl;
	fout << "resEnd " << resEnd << endl;
	fout << "tmStart " << tmStart << endl;
	fout << "tmEnd " << tmEnd << endl;
	fout << "numModels " << members.size() << endl;
	fout << "thread " << members[0]->getThread() << endl;
	fout << "XShift " << members[0]->getXShift() << endl;
	fout << "ZShift " << members[0]->getZShift() << endl;
	fout << "axialRot " << members[0]->getAxialRotation() << endl;
	fout << "crossAng " << members[0]->getCrossingAngle() << endl;
	fout << "energy " << members[0]->getEnergy() << endl;
	fout << "hbonds ";
	vector<string> & hbonds = members[0]->getHbonds();
	for(int i = 0; i < hbonds.size(); i++) {
		fout << hbonds[i] << ";";
	}
	fout << endl;
	// print details of cluster members
	fout << "MODEL num thread   XShift   ZShift axialRot crossAng numhb energy " << endl;
	for(int i = 0; i < members.size(); i++) {
		char line[1000];
		sprintf(line,"MODEL %3d %6d %8.3f %8.3f %8.3f %8.3f %5d %f",i,members[i]->getThread(),members[i]->getXShift(), members[i]->getZShift(), members[i]->getAxialRotation(),members[i]->getCrossingAngle(),members[i]->getNumHbonds(),members[i]->getEnergy());
		//fout << "MODEL " << i << " " << members[i]->getThread() << " " << members[i]->getXShift() << " " << members[i]->getZShift() << " " << members[i]->getAxialRotation() << " " << members[i]->getCrossingAngle() << " " << members[i]->getEnergy() << " " << members[i]->getNumHbonds() << endl; 
		fout << line << endl;
	}
	fout.close();
}

void HelixDimerCluster::makePse(int id) {

	char filename[1000];
	///sprintf(filename,"%s/%s_%02d",outputDir.c_str(),uniprotAccession.c_str(),id);
	sprintf(filename,"%s_%02d",uniprotAccession.c_str(),id);


	ofstream fout;

	char scriptfilename[1000];
	sprintf(scriptfilename,"%s/%s.inp",outputDir.c_str(),filename);
	
	char pdbfilename[1000];
	sprintf(pdbfilename,"%s.pdb",filename);

	char psefilename[1000];
	sprintf(psefilename,"%s.pse",filename);

	// PyMOL script
	/*
	bg_color white
	Load MODEL_000.pdb
	Rotate X, 90
	Color green, chain B and element C
	Color cyan, chain B and element C
	Show cartoon
	Show stick, resi 28+31+32+35+36+39+40
	Distance hb000, ///A/35/HA1, ///B/32/O
	Color black, hb000
	Distance hb001, ///B/35/HA1, ///A/32/O
	Color black, hb001
	Distance hb002, ///A/36/HA, ///B/35/O
	Color black, hb002
	Distance hb003, ///B/36/HA, ///A/35/O
	Color black, hb003
	Distance hb004, ///A/39/HG1, ///B/36/O
	Color black, hb004
	Distance hb005, ///B/39/HG1, ///A/36/O
	Color black, hb005
	Zoom all, 3
	Set label_color, black
	Save MODEL_000.pse
	Hide labels
	Ray 400,600
	Png MODEL_000_view1.png
	Rotate Y, 90
	Ray 500,600
	Png MODEL_000_view2.png
*/

	fout.open(scriptfilename);
	fout << "bg_color white" << endl;
	fout << "load " << pdbfilename << endl;
	fout << "rotate X, 90, all, 0" << endl;
	fout << "color green, chain B and element C" << endl;
	fout << "color cyan, chain B and element C" << endl;
	fout << "show cartoon" << endl;

	// find the interface residue numbers from the seq, resStart and interface mask
	string interface = members[0]->getInterface();

	vector<int> interfaceRes;
	for(int i = 0; i < interface.length();i++) {
		if(interface[i] == '1' || interface[i] == '2') {
			interfaceRes.push_back(resStart + i);
		}
	}
	//fout << "show stick, resi " << interfaceRes[0];

	// form the strings
	/*
	select interfaceA, chain A and resi 6+9+10+13+14+16+18+20+21
	select interfaceB, chain B and resi 6+9+10+13+14+16+18+20+21
	show stick, interfaceA interfaceB
	*/

	stringstream ss;
	ss <<  interfaceRes[0];
	for(int i = 1; i < interfaceRes.size(); i++) {
		ss << "+" << interfaceRes[i];
	}

	fout << "select interfaceA, chain A and resi " << ss.str() << endl;
	fout << "select interfaceB, chain B and resi " << ss.str() << endl;
	fout << "show stick, interfaceA interfaceB" << endl;

	// A,12,ALA,HA1:B,34,LEU,O=2.6;A,12,ALA,HA1:B,34,LEU,O=3.5
	vector<string> hbondTokens  = members[0]->getHbonds();
	
	FormatConverter fc;
	fc.setNamespaces("CHARMM22","PDB2");
	for(int i = 0; i < hbondTokens.size(); i++) {
		vector<string> bondData = MslTools::tokenize(hbondTokens[i],":=");
		string chain1,resNum1,resName1,atom1;
		vector<string> donorData = MslTools::tokenize(bondData[0],",");
		chain1 = donorData[0];
		resNum1 = donorData[1];
		resName1 = donorData[2];
		if(resName1 == "HSD" || resName1 == "HSE" || resName1 == "HSP") {
			resName1 = "HIS";
		}
		atom1 = donorData[3];
		atom1 = fc.getAtomName(atom1,resName1);
		string chain2,resNum2,resName2,atom2;
		vector<string> acceptorData = MslTools::tokenize(bondData[1],",");
		chain2 = acceptorData[0];
		resNum2 = acceptorData[1];
		resName2 = acceptorData[2];
		if(resName2 == "HSD" || resName2 == "HSE" || resName2 == "HSP") {
			resName2 = "HIS";
		}
		atom2 = acceptorData[3];
		atom2 = fc.getAtomName(atom2,resName2);
		char distLine[500];
		sprintf(distLine,"distance hb%03d, ///%s/%s/%s, ///%s/%s/%s\ncolor black, hb%03d",i,chain1.c_str(),resNum1.c_str(),atom1.c_str(),chain2.c_str(),resNum2.c_str(),atom2.c_str(),i);
		fout << distLine << endl;
	}

	fout << "zoom all, 3" << endl;
	fout << "set label_color, black" << endl;
	fout << "save " << psefilename << endl;
	fout << "hide labels" << endl;
	fout << "ray 400,600" << endl;
	fout << "png " <<  filename << "_view1.png" << endl;
	fout << "rotate Y, 90, all, 0" << endl;
	fout << "ray 500,600" << endl;
	fout << "png " << filename << "_view2.png" << endl;
//TODO : once pymol is fixed uncomment this part - update the correct path before invoking the .inp
/*
	char command [1000];
	sprintf(command,"pymol -cqd @%s",scriptfilename);
	system(command);
	*/
}

/******************************************************************************************************************************************************************************/










struct Options {
	// Required
	string fullSequence;

	// optional
	string backboneCrd;
	string logFile;
	string pdbOutputDir;

	int tmStart;
	int tmEnd;

	string helixGeoFile;

	string rulesFile;

	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;
	string rotLibFile;

	// side-chain repack variable (optional)
	int MCCycles;
	int MCMaxRejects;

	bool verbose;
	bool greedy;
	int greedyCycles;
	int seed;

	// protein information (optional)
	string uniprotName;
	string uniprotAccession;

	// clustering options (optional)
	double rmsdCutoff;
	bool clusterSolutions;
	bool printAllCrds;
	bool printAxes;
	bool printTermEnergies;

	int fullSequenceStart;

	// the actual AAs being modeled
	int startResNum; 
	int endResNum;

	int threadStart;
	int threadEnd;

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser
	string rerunConf; // data for a configuration file that would rerun the job as the current run

	string configfile;
};

struct HbondInfo{
	string chain;
	string resA;
	string resB;
	string atomName;
	double distBreak;

	HbondInfo(string _c, string _resA, string _resB, string _aName, string _dist) {
		chain = _c;
		resA = _resA;
		resB = _resB;
		atomName = _aName;
		distBreak = MslTools::toDouble(_dist);
	}

	bool isValid(System & _sys) {
		string posA = chain + "," + resA;
		string posB = "";
		if(chain == "A") {
			posB = "B," + resB;
		} else {
			posB = "A," + resB;
		}
		if(!_sys.positionExists(posA) || !_sys.positionExists(posB) ) {
			return false;
		}

		if(atomName == "HA1") {
			if(!_sys.identityExists(posA + ",GLY" )) {
				return false;
			}
		}
		return true;
	}

	void print(ofstream & fout) {
		fout << chain << " " << resA << " " << resB << " " << atomName << " " << distBreak << endl;
	}
};

Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);

// Just add 10 U(0,1) uniform random variables, offset by 0.5 to make mean = 0 and divide by variance = (10 * var(U(0,1))) 
double getStandardNormal(RandomNumberGenerator& RNG) {
	double retVal = 0.0;
	for(int i = 0; i < 10; i ++) {
		retVal += RNG.getRandomDouble();
	}
	return (retVal/10.0 - 0.5) * 1.2;
}

void printOptions(Options & _op, ofstream & _fout) {
	_fout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;

	_fout << "backboneCrd " << _op.backboneCrd << endl;
	_fout << "logFile " << _op.logFile << endl;
	_fout << "pdbOutputDir " << _op.pdbOutputDir << endl;

	_fout << "fullSequence " << _op.fullSequence << endl;
	_fout << "tmStart " << _op.tmStart << endl;
	_fout << "tmEnd " << _op.tmEnd << endl;

	_fout << "helixGeoFile " << _op.helixGeoFile << endl;
	_fout << "rulesFile " << _op.rulesFile << endl;

	_fout << "topFile " << _op.topFile << endl;
	_fout << "parFile " << _op.parFile << endl;
	_fout << "solvFile " << _op.solvFile << endl;
	_fout << "hBondFile " << _op.hBondFile << endl;
	_fout << "rotLibFile " << _op.rotLibFile << endl;

	_fout << "MCCycles " << _op.MCCycles << endl;
	_fout << "MCMaxRejects " << _op.MCMaxRejects << endl;

	_fout << "verbose " << _op.verbose << endl;
	_fout << "greedyOptimizer " << _op.greedy << endl;
	_fout << "greedyCycles " << _op.greedyCycles << endl;
	_fout << "seed " << _op.seed << endl;

	_fout << "uniprotName " << _op.uniprotName << endl;
	_fout << "uniprotAccession " << _op.uniprotAccession << endl;

	_fout << "rmsdCutoff " << _op.rmsdCutoff << endl;
	_fout << "clusterSolutions " << _op.clusterSolutions << endl;
	_fout << "printAllCrds " << _op.printAllCrds << endl;
	_fout << "printAxes " << _op.printAxes << endl;
	_fout << "printTermEnergies " << _op.printTermEnergies << endl;

	_fout << "fullSequenceStart " << _op.fullSequenceStart << endl;

	_fout << "startResNum " << _op.startResNum << endl;
	_fout << "endResNum " << _op.endResNum << endl;

	_fout << "threadStart " << _op.threadStart << endl;
	_fout << "threadEnd " << _op.threadEnd << endl;

	if(_op.configfile != "") {
		_fout << "configfile " << _op.configfile << endl;
	}

	_fout << endl;

}

void c2Symmetry(AtomPointerVector & _apvA, AtomPointerVector & _apvB) {
	
	/* Faster code 
	for (uint i=0; i < _apvA.size(); i++) {
			_apvB[i]->copyAllCoor(*_apvA[i]);
			vector<CartesianPoint*>& bCoors = _apvB[i]->getAllCoor();

			for(uint j = 0; j < bCoors.size(); j++) {
				bCoors[j]->setX(0 - bCoors[j]->getX());
				bCoors[j]->setY(0 - bCoors[j]->getY());
			}
					
		}
	*/

	// Set coordinates of chain A to chain B
	for (uint i=0; i < _apvA.size(); i++) {
		_apvB[i]->copyAllCoor(*_apvA[i]);
	}

	// Rotation matrix for 180 degrees
	// flips the sign on the x and y coordinates
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


unsigned int CAOfClosestApproach(Chain & _chainA, Chain & _chainB) {
	double minCAtoCAdist = MslTools::doubleMax;
	string closestCA_A = "";
	string closestCA_B = "";
	
	unsigned int closestAtomIndex = 0;
	for (unsigned int i=0; i < _chainA.positionSize(); i++) {
		Position& posA = _chainA.getPosition(i);
		Position& posB = _chainB.getPosition(i);
		if (posA.atomExists("CA") && posB.atomExists("CA")) {
			Atom& posACA = posA.getLastFoundAtom();
			Atom& posBCA = posB.getLastFoundAtom();
			double dist = posACA.distance(posBCA);
			if (dist < minCAtoCAdist) {
				minCAtoCAdist = dist;
				closestAtomIndex = i;
			}
		}
	}

	return _chainA(closestAtomIndex).getResidueNumber();
}

map<string, unsigned int> interfaceResidueCheck(AtomPointerVector & _chainA, AtomPointerVector & _chainB) {
	map<string, unsigned int> atomsWithIn4AOfPosition;
	for (uint i=0; i < _chainA.size(); i++) {
		if (_chainA[i]->getName() != "CA" && _chainA[i]->getName() != "C" && _chainA[i]->getName() != "N" && _chainA[i]->getName() != "HN" && _chainA[i]->getName() != "O") {
			for (uint j=0; j < _chainB.size(); j++) {
				if (_chainB[j]->getName() != "CA" && _chainB[j]->getName() != "C" && _chainB[j]->getName() != "N" && _chainB[j]->getName() != "HN" && _chainB[j]->getName() != "O") {
					if (_chainA[i]->distance(*_chainB[j]) < 4.0) {
						if (atomsWithIn4AOfPosition.find(_chainA[i]->getIdentityId()) != atomsWithIn4AOfPosition.end()) {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()]++;
						}
						else {
							atomsWithIn4AOfPosition[_chainA[i]->getIdentityId()] = 1;
						}
					}
				}
			}
		}
	}
	
	//for(map<string, unsigned int>::iterator it=atomsWithIn4AOfPosition.begin(); it != atomsWithIn4AOfPosition.end(); it++) {
	//	_fout << "pos: " << it->first << " count: " << it->second << endl;
	//}

	return atomsWithIn4AOfPosition;
}

string convertToPolymerSequence(string _seq, int _startResNum) {
	// convert a 1 letter _sequence like AIGGG and startResNum = 32 to 
	// A:{32}ALA ILE GLY GLY GLY
	// B:{32}ALA ILE GLY GLY GLY
	string ps = "";
	for(string::iterator it = _seq.begin(); it != _seq.end();it++ ) {
		stringstream ss;
		ss << *it;
		string resName = MslTools::getThreeLetterCode(ss.str());
		if(resName == "HIS") {
			ps = ps + " HSE";
		} else {
			ps = ps + " " + resName;
		}
	}
	ps = ":{" + MslTools::intToString(_startResNum) + "} " + ps;
	return "A" + ps + "\nB" + ps;
}

void reThreadResidues(vector<Position*> & positions, int offset) {
	for(int i = 0; i < positions.size(); i++) {
		int resNum = positions[i]->getResidueNumber();
		positions[i]->setResidueNumber(resNum + offset);
	}		
}

bool hydrogenBondCheck(System & _sys, vector<string> _parsedGeoInformation, double & _xShiftStart) {
	
	vector<HbondInfo*> hbondList;
	for(unsigned k = 5; k < _parsedGeoInformation.size(); k+=5 ) {
		HbondInfo * t = new HbondInfo(_parsedGeoInformation[k],_parsedGeoInformation[k+1],_parsedGeoInformation[k+2],_parsedGeoInformation[k+3],_parsedGeoInformation[k+4]);
		if(t->isValid(_sys)) {
			hbondList.push_back(t);
		} else {
			delete t;
		}
	}


	/*
	for(int i = 0; i < hbondList.size(); i++) {
		hbondList[i]->print(fout);
	}
	*/
	if(hbondList.size() < 4) {
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return false;
	} else {
		_xShiftStart = hbondList[hbondList.size() - 4]->distBreak - 0.1;
		for(int i = 0; i < hbondList.size(); i++) {
			delete hbondList[i];
		}
		return true;
	}

}

void renumberResidues(System& _sys, int _startResNum) {
	// just set the index + startResNum of each position as its residue number
	vector<Position*>& positions = _sys.getPositions();
	for(int i = 0; i < positions.size(); i++) {
		int resNum = positions[i]->getIndexInChain();
		positions[i]->setResidueNumber(resNum + _startResNum);
	}	
}

bool rulesCheck(System & _sys, string _geoIndex, map<int, string> _rulesMap) {

	int idx = MslTools::toInt(_geoIndex);
	if (_rulesMap.find(idx) == _rulesMap.end()) {
		return true;
	}

	vector<string> parsedRules = MslTools::tokenize(_rulesMap[idx], ",");
	// read over all the rules for a given model
	for (uint m = 1; m < parsedRules.size(); m++) {
		int delimiter  = parsedRules[m].find_first_of(":");
		string position = parsedRules[m].substr(0,1)+","+parsedRules[m].substr(1,delimiter-1); // makes string A,35
		string rule = parsedRules[m].substr(delimiter+1,(parsedRules[m].length()-delimiter-1)); 
		if(!_sys.positionExists(position)) {
			continue;
		}
		//fout << "position " << position << " must be " << rule << " , " << " is actually " << _sys.getPosition(position).getResidueName() << endl;

		// Check if there is a ! "not"
		size_t findNot = rule.find("!");

		// if there is a "not"
		if (findNot == 0) { // ! will only appear in position 0
			rule = rule.substr(2);
			rule = MslTools::trim(rule, "]"); // residues that cannot exist in given position
			// look at position given in rule	
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					return false;
				}
			}
			
		}
		else { // no "not" found
			rule = rule.substr(1);
			rule = MslTools::trim(rule, "]"); // residues that must exist in given position
			// look at position given in rule	
			string resName = MslTools::getOneLetterCode(_sys.getPosition(position).getResidueName());
			bool found = false;
			for (uint n=0; n < rule.length(); n++) {
				if (rule.substr(n, 1) == resName) {
					found = true;
					break;
				}
			}
			if(!found) {
				return false;
			}
		}
	}
	
	return true;
}

void transformation(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, CartesianPoint & _ori, CartesianPoint & _xAxis, CartesianPoint & _zAxis, double _zShift, double _axialRotation, double _crossingAngle, double _xShift, Transforms & _trans) {
	
	//====== Z Shift (Crossing Point) ======
	CartesianPoint zShiftCP(0.0, 0.0, _zShift);
	_trans.translate(_chainA, zShiftCP);

	//===== Axial Rotation ======
	_trans.rotate(_chainA, _axialRotation, _ori, _zAxis);

	//====== Local Crossing Angle ======
	_trans.rotate(_chainA, (_crossingAngle/2.0), _ori, _xAxis);
	_trans.rotate(_axisA, (_crossingAngle/2.0), _ori, _xAxis);

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*_xShift/2.0), 0.0, 0.0);
	_trans.translate(_chainA, interDistVect);
	_trans.translate(_axisA, interDistVect);

	c2Symmetry(_chainA, _chainB);
	c2Symmetry(_axisA, _axisB);
}

void repackSideChains(System & _sys, SelfPairManager & _spm, bool _greedy, int _greedyCycles) {

	_spm.calculateEnergies();

	if(!_greedy) {
		_spm.setMCOptions(1000.0, 0.5, 20000, MonteCarloManager::EXPONENTIAL, 2000, 100, 0.01);
		_spm.setRunDEE(false);
		_spm.setRunEnum(false);
		_spm.setRunSCMF(true);
		_spm.setRunSCMFBiasedMC(true);
		_spm.runOptimizer();
	} else {
		_spm.runGreedyOptimizer(_greedyCycles);
	}
}

void backboneMovement(AtomPointerVector & _chainA, AtomPointerVector & _chainB, AtomPointerVector & _axisA, AtomPointerVector & _axisB, Transforms _trans, double _deltaMove, unsigned int moveType) {

	 if (moveType == 0) {
		// Z Shift
		CartesianPoint translateA = _axisA(1).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove; // unit vector of helical _axis times the amount to shift by

		_trans.translate(_chainA, translateA);

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 1) {
		// Axial Rotation
		_trans.rotate(_chainA, (_deltaMove), _axisA(0).getCoor(), _axisA(1).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else 	if (moveType == 2) {
		// Crossing Angle 
		_trans.rotate(_chainA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());
		_trans.rotate(_axisA, (_deltaMove * 0.5), _axisA(0).getCoor(), _axisB(0).getCoor());

		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else if (moveType == 3) {
		// XShift
		// Helix A interhelical distance
		CartesianPoint translateA = _axisB(0).getCoor() - _axisA(0).getCoor(); // vector minus helical center 
		translateA = translateA.getUnit() * _deltaMove * -0.5; // unit vector of helical axis times the amount to shift by

		_trans.translate(_chainA, translateA);
		_trans.translate(_axisA, translateA);

		// Helix B interhelical distance
		c2Symmetry(_chainA, _chainB);
		c2Symmetry(_axisA, _axisB);

	} else {
		cerr << "Unknown moveType " << moveType << " in backboneMovement. Should be 0-3 " << endl;
	}
}

double computeMonomerEnergy(System & _sys, Transforms & _trans, RandomNumberGenerator & _RNG, SystemRotamerLoader & _sysRot, SelfPairManager & _spm) {

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector &chainA = _sys.getChain("A").getAtomPointers();
	AtomPointerVector &chainB = _sys.getChain("B").getAtomPointers();

	/******************************************************************************
	 *                     === SHIFT HELICES 500 A APART ===
	 ******************************************************************************/
	double xShift = 500.0;

	//====== X shift (Interhelical Distance) =======
	CartesianPoint interDistVect;
	interDistVect.setCoor((-1.0*xShift/2.0), 0.0, 0.0);
	_trans.translate(chainA, interDistVect);

	c2Symmetry(chainA, chainB);

	// Repack side chains
	repackSideChains(_sys, _spm, true, 1);

	return _spm.getMinBound()[0];

}

void readRulesFile(string _fileName, map<int, string> & _rulesFileMap) {
	ifstream file;
	file.open(_fileName.c_str());
	if(!file.is_open()) {
		cerr << "Unable to open " << _fileName << endl;
		exit(0);
	}

	string tmpRulesLine;

	while(file) {
		getline(file, tmpRulesLine);
		if (tmpRulesLine.length() > 1) {
			vector<string> token = MslTools::tokenizeAndTrim(tmpRulesLine,",");
			_rulesFileMap[MslTools::toInt(token[0])] = tmpRulesLine;
		}
	}
	file.close();
}

void moveZCenterOfCAMassToOrigin(AtomPointerVector& _apV) {

	AtomSelection sel(_apV);
	AtomPointerVector & caApV = sel.select("name CA");
	double zShift = 0.0;
	for(int i = 0; i < caApV.size(); i++) {
		zShift += (caApV[i]->getCoor()).getZ();
	}
	zShift = -1.0 * zShift/double(caApV.size());
	//fout << x << " " << y << " " << pt << " " << caApV.size() << endl;

	for(int i = 0; i < _apV.size(); i++) {
		CartesianPoint& pt = _apV[i]->getCoor();
		pt.setZ(pt.getZ() +  zShift);
	}

}

vector<string> getInterHelicalHbonds(EnergySet* & _ESet) {
	unsigned int numHbonds = 0;
	// Why are we doing this?
	//_ESet->setAllTermsInactive();
	//_ESet->setTermActive("SCWRL4_HBOND", true);
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];

	vector<string> hbonds;
	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		if(atoms[0]->getChainId() == atoms[2]->getChainId()) {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		if(e < 0) {
			//_fout << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
			hbonds.push_back(atoms[0]->getAtomOfIdentityId() + ":" + atoms[2]->getAtomOfIdentityId() + "=" + MslTools::doubleToString(atoms[0]->distance(*atoms[2])));
			numHbonds++;
		}
	}
	// Why are we doing this?
	//_ESet->setAllTermsActive();
	//_ESet->setTermActive("CHARMM_ELEC", false);
	return hbonds;
}

void readGeometryFile(string _filename, vector<string>& _fileVec) {
	ifstream file;
	file.open(_filename.c_str()); 
	if(!file.is_open()) {
		cerr << "Unable to open " << _filename << endl;
		exit(0);
	}

	string parameterList;

	while(file) {
		getline(file, parameterList);
		_fileVec.push_back(parameterList);
	}
	file.close();
}


void clusterSolutions(vector<HelixDimerCluster*>& clusters,vector<HelixDimer*>& _structures, double _rmsdCutoff, string _origSeq, string _builtSeq, Options& opt) {
	for(int i  = 0; i < _structures.size(); i++) {
		bool diff = true;
		AtomPointerVector& a1 = _structures[i]->getAtomPointers();
		AtomSelection sel1(a1);
		AtomPointerVector& ca1 = sel1.select("name CA");
		//cout << i << endl;
		for(int j = 0; j < clusters.size(); j++) {
			AtomPointerVector& a2 = clusters[j]->getAtomPointers();
			AtomSelection sel2(a2);
			AtomPointerVector& ca2 = sel2.select("name CA");
			//cout << _structures[i]->getName() << " " << _structures[j]->getName() << " ";
			double r = ca1.rmsd(ca2);
			if(r < _rmsdCutoff) {
				clusters[j]->addHelixDimer(_structures[i]);
				diff = false;
				break;
			}
		}
		if(diff) {
			clusters.push_back(new HelixDimerCluster(_structures[i]));
			clusters.back()->setDetails(_origSeq, _builtSeq, opt.uniprotName,opt.uniprotAccession,opt.pdbOutputDir,opt.startResNum,opt.endResNum,opt.tmStart,opt.tmEnd);
		}
	}

}

map<string,double> getEnergyByTerm(EnergySet* _eSet) {
	// get all terms
	map<string,double> eByTerm;
	map<string,vector<Interaction*> > * allTerms = _eSet->getEnergyTerms();
	for(map<string,vector<Interaction*> >::iterator it = allTerms->begin(); it != allTerms->end(); it++) {
		if(_eSet->isTermActive(it->first)) {
			eByTerm[it->first] =  _eSet->getTermEnergy(it->first);
		}
	}
	return eByTerm;
}

int main(int argc, char *argv[]) {

	time(&startTime);	

	/******************************************************************************
	 *                 === PARSE THE COMMAND LINE OPTIONS ===
	 ******************************************************************************/
	Options defaults;
  
  	Options opt = parseOptions(argc, argv, defaults);
  	if (opt.errorFlag) {
  		cerr << endl;
  		cerr << "The program terminated with errors:" << endl;
  		cerr << endl;
  		cerr << opt.errorMessages << endl;
  		cerr << endl;
  		cerr << opt.OPerrors << endl;
  
  		usage();
  		exit(1);
  	}	

	ofstream fout;
	fout.open(opt.logFile.c_str());

	if(!fout.is_open()) {
		cerr << "Unable to open " << opt.logFile << endl;
		exit(0);
	}

	printOptions(opt, fout);

	/******************************************************************************
	 *                     === READ IN GEOMETRY FILE ===
	 ******************************************************************************/
	vector<string> fileVec;
	readGeometryFile(opt.helixGeoFile,fileVec);

	// Import sequence, determine length
	string originalTMSeq = opt.fullSequence.substr(opt.startResNum-opt.fullSequenceStart,opt.endResNum - opt.startResNum + 1);
	unsigned int sequenceLength = originalTMSeq.length(); 
	// cant do sequences of size less than 4
	if(sequenceLength < 4) {
		cerr << "Sequence " << originalTMSeq << "is too small (should be >= 4 AA long)" << endl; 
		exit(0);
	}

	string modelledTMSeq ="" ; // replace P with A
	string prolineMask = "";
	for(int i = 0; i < originalTMSeq.length(); i++) {
		char AA = originalTMSeq[i];
		if(AA == 'P') {
			modelledTMSeq += "A";
			prolineMask += "1";
		} else {
			modelledTMSeq += AA;
			prolineMask += "0";
		}
	}
	/**********************************************************************************
	*
	*    printProteinOutFile
	*
	**********************************************************************************/

	ofstream pout;
	string poutName  = opt.pdbOutputDir + "/" + opt.uniprotAccession + ".out";
	pout.open(poutName.c_str());
	if(!pout.is_open()) {
		cerr << "Unable to open " << poutName << endl;
		exit(0);
	}

	pout << "uniprotName " << opt.uniprotName << endl;
	pout << "uniprotAccession " << opt.uniprotAccession << endl; 
	pout << "sequence " << opt.fullSequence << endl;
	pout << "sequenceStart " << opt.fullSequenceStart << endl;
	pout << "originalSeq " << originalTMSeq << endl;
	pout << "modelledSeq " << modelledTMSeq << endl;
	pout << "prolineMask " << prolineMask << endl;
	pout << "resStart " << opt.startResNum << endl;
	pout << "resEnd " << opt.endResNum << endl;
	pout << "tmStart " << opt.tmStart << endl;
	pout << "tmEnd " << opt.tmEnd << endl;
	pout << "CaH-TMversion " << programVersion << endl;


	//cout << modelledTMSeq << endl;

	string sequence = convertToPolymerSequence(modelledTMSeq,1); // so that the 4th residue will be the middle one (35th) on the GLY 69 backbone
	PolymerSequence PS(sequence); 

	// Create system with sequence - a string with the following format
	// A:{startingResNum} ALA ILE ...\n
	// B:{startingResNum} ALA ILE ...

	/******************************************************************************
	 *                     === DECLARE SYSTEM ===
	 ******************************************************************************/
	System sys;
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	CSB.setSolvent("CHEX");

	if(!CSB.buildSystem(PS)) {
		cerr << "Unable to build system from " << sequence << endl;
		exit(0);
	} else {
		//fout << "CharmmSystem built for sequence" << endl;
	}

	Chain & chainA = sys.getChain("A");
	Chain & chainB = sys.getChain("B");

	// Set up chain A and chain B atom pointer vectors
	AtomPointerVector & apvChainA = chainA.getAtomPointers();
	AtomPointerVector & apvChainB = chainB.getAtomPointers();

	vector<Position*>& positions = sys.getPositions();


	// Read in Gly-69 to use as backbone coordinate template
	CRDReader cRead;
	cRead.open(opt.backboneCrd);
	if(!cRead.read()) {
		fout << "Unable to read " << opt.backboneCrd << endl;
		exit(0);
	}
	cRead.close();
	AtomPointerVector& glyAPV = cRead.getAtomPointers();

	// Reference points for Helices
	CartesianPoint ori(0.0,0.0,0.0);
	CartesianPoint zAxis(0.0,0.0,1.0);
	CartesianPoint xAxis(1.0,0.0,0.0);

	// Objects used for transformations
	Transforms trans; 
	trans.setTransformAllCoors(true); // transform all coordinates (non-active rotamers)
	trans.setNaturalMovements(true); // all atoms are rotated such as the total movement of the atoms is minimized

	// Random Number Generator
	RandomNumberGenerator RNG1;
	RNG1.setSeed(opt.seed);

	// Read Rotamer Library File
	SystemRotamerLoader sysRot(sys, opt.rotLibFile);
	sysRot.defineRotamerSamplingLevels();

	/******************************************************************************
	 *                     === COPY BACKBONE COORDINATES ===
	 ******************************************************************************/
	sys.assignCoordinates(glyAPV,false);
	sys.buildAtoms();

	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	hb.buildInteractions(30);

	/******************************************************************************
	 *                     === INITIAL VARIABLE SET UP ===
	 ******************************************************************************/
	EnergySet* Eset = sys.getEnergySet();

	// Set all terms active, besides Charmm-Elec
	Eset->setAllTermsActive();
	Eset->setTermActive("CHARMM_ELEC", false);

	/******************************************************************************
	 *                     === HELICAL AXIS SET UP ===
	 ******************************************************************************/

	string axis = "\
ATOM      1  O   DUM A   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      2  Z   DUM A   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
ATOM      3  O   DUM B   1       0.000   0.000   0.000  1.00  0.00           P\n\
ATOM      4  Z   DUM B   1       0.000   0.000   1.000  1.00  0.00           O\n\
TER\n\
END";
	

	PDBReader readAxis;
	if(!readAxis.read(axis)) {
		cerr << "Unable to read axis" << endl;
		exit(0);
	}

	System helicalAxis;
	//helicalAxis.readPdb(opt.helicalAxisPdbFile);
	helicalAxis.addAtoms(readAxis.getAtomPointers());

	AtomPointerVector &axisA = helicalAxis.getChain("A").getAtomPointers();
	AtomPointerVector &axisB = helicalAxis.getChain("B").getAtomPointers();

	helicalAxis.saveCoor("originState");

	// Declare SelfPairManager and Set Seed
	SelfPairManager spm;
	spm.seed(RNG1.getSeed()); 

	spm.setOnTheFly(opt.greedy);
	spm.setVerbose(opt.verbose);

	/******************************************************************************
	 *                  === LOAD ROTAMERS & SET-UP SPM ===
	 ******************************************************************************/
	for (uint k=0; k < sys.positionSize(); k++) {
		Position &pos = sys.getPosition(k);

		if (pos.getResidueName() != "GLY" && pos.getResidueName() != "ALA" && pos.getResidueName() != "PRO") {
			if (!sysRot.loadRotamers(&pos, pos.getResidueName(), "SL95.00")) { 
				cerr << "Cannot load rotamers for " << pos.getResidueName() << endl;
			}
		}
	}

//	spm.saveEnergiesByTerm(true);
	spm.setSystem(&sys);

	map<string,double> monomerEnergyByTerm;
	double monomerEnergy;

	// Store monomer energy by term
	if(opt.printTermEnergies) {
		monomerEnergy = computeMonomerEnergy(sys, trans, RNG1, sysRot, spm);
		sys.setActiveRotamers(spm.getMinStates()[0]);
		sys.calcEnergy();
		monomerEnergyByTerm = getEnergyByTerm(sys.getEnergySet());

		ofstream meOut;
		string meOutName  = opt.pdbOutputDir + "/" + opt.uniprotAccession + "_monomer.energy";
		meOut.open(meOutName.c_str());
		if(!meOut.is_open()) {
			cerr << "Unable to open " << meOutName << endl;
			exit(0);
		}
		for(map<string,double>::iterator it = monomerEnergyByTerm.begin(); it != monomerEnergyByTerm.end(); it++) {
			meOut << it->first << " " << it->second << endl;
		}
		meOut.close();

	}

//	fout <<  "Monomer Energy " << monomerEnergy << endl;
//	cout <<  "Monomer Energy " << computeMonomerEnergy(sys, trans, RNG1, sysRot, spm) << endl;
//	cout << spm.getSummary(spm.getMinStates()[0]) << endl;

	/******************************************************************************
	 *                  === READ IN RULES ===
	 ******************************************************************************/
	map<int, string> rulesFileMap;
	readRulesFile(opt.rulesFile, rulesFileMap);


	// To store the structures in case we are clustering
	vector<HelixDimer*> structures;

	// We will need a format converter to convert to PDB names
	FormatConverter fc;



	map<string,double> energyByTerm;
	/******************************************************************************
	 *              === LOOP OVER ALL POSSIBLE INTERFACE POSITIONS ===
	 ******************************************************************************/
	for(int j = opt.threadStart; j <= opt.threadEnd; j++) {

		//bool computedMonomerEnergy = false;
		chainA.renumberChain(j);
		chainB.renumberChain(j);
		//renumberResidues(sys, 32 - j);

		/******************************************************************************
		 *                     === COPY BACKBONE COORDINATES ===
		 ******************************************************************************/
		sys.wipeAllCoordinates();
		sys.assignCoordinates(glyAPV,false);		
		sys.buildAllAtoms();
		sys.saveCoor("currentTMFragment");

		/******************************************************************************
		 *              === LOOP OVER ALL LINES IN HELIX GEOMETRY FILE ===
		 ******************************************************************************/
		for (uint i=0; i<fileVec.size()-1; i++) {

			chainA.renumberChain(j);
			chainB.renumberChain(j);
			//renumberResidues(sys, 32 - j);

			helicalAxis.applySavedCoor("originState");
			sys.applySavedCoor("currentTMFragment");

			// Parse parameter file line
			fileVec[i] = MslTools::trim(fileVec[i], "\t\n\r");
			vector<string> parsedGeoInformation = MslTools::tokenize(fileVec[i], " ");

			fout << parsedGeoInformation[0] << " A:" << j << " B:" << j << endl;

			if (parsedGeoInformation.size() % 5 != 0) {
				cerr << "File line: " << i << " is of incompatible length" << endl;
				exit(1);
			}
			
			// index axialRot crossingAngle zShift xShift ChainDonor DonorResNum AcceptorResNum DonorAtom DistBondBreaks
			// 00001 75 -35.0 1.0 6.8 A 7 8 HA2 9.6 
			// vector length = 5, 0 hbonds
			// vector length = 10, 1 hbonds
			// vector length = 15, 2 hbonds

			/******************************************************************************
			 *                     === HYDROGEN BOND COUNT CHECK ===
			 ******************************************************************************/
			double xShiftStart = 0;
			if(!hydrogenBondCheck(sys, parsedGeoInformation, xShiftStart)) {
				fout << "less than 4 hbonds" << endl;

				// renumber residues - rethreading
				//reThreadResidues(positions,-1);
				fout << endl;
				continue;
			}

			/******************************************************************************
			 *                     === CHECK SEQUENCE RULES ===
			 ******************************************************************************/
			if(!rulesCheck(sys, parsedGeoInformation[0], rulesFileMap)) {
				fout << "sequence does not conform with given rules" << endl;

				// renumber residues - rethreading
				//reThreadResidues(positions,-1);
				fout << endl;
				continue;
			}
			/******************************************************************************
			 *                     === INITIAL STARTING POSITION ===
			 ******************************************************************************/
			double xShift = xShiftStart;
			double crossingAngle = MslTools::toDouble(parsedGeoInformation[2]);
			double axialRotation = MslTools::toDouble(parsedGeoInformation[1]);
			double zShift = MslTools::toDouble(parsedGeoInformation[3]);

			fout << "Starting parameters: crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift;
			fout << " xShift: " << xShift << endl;

			// Transform helices to initial starting position
			transformation(apvChainA, apvChainB, axisA, axisB, ori, xAxis, zAxis, zShift, axialRotation, crossingAngle, xShift, trans);
			// Optimizatize Initial Starting Position
			repackSideChains(sys, spm, opt.greedy, opt.greedyCycles);

			sys.setActiveRotamers(spm.getMinStates()[0]);
			double currentEnergy = spm.getMinBound()[0];
			sys.saveAltCoor("savedBestState");
			helicalAxis.saveAltCoor("BestAxis");

			/******************************************************************************
			 *              === CALCULATE MONOMER ENERGY ===
			 ******************************************************************************/
		/*	
			if(!computedMonomerEnergy) {
				monomerEnergy = computeMonomerEnergy(sys, trans, RNG1, sysRot, spm);
				fout << "The Monomer Energy is: " << monomerEnergy << endl;
				sys.applySavedCoor("savedBestState");
				//fout << spm.getSummary(spm.getMinStates()[0]) << endl;
				computedMonomerEnergy = true;
			} 
		*/

			fout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;

			// filter based on previous run
			// if the outermost energy is still > 100 it is likely we will not get <0 energy with this model

			if(currentEnergy-monomerEnergy > 100) {
				fout << "discard model since deltaE Exceeds 100	at dOut" << endl;
				fout << endl;
				continue;
			}

			/******************************************************************************
			 *                     === X SHIFT REPACKS ===
			 ******************************************************************************/
			double bestEnergy = currentEnergy;
			double savedXShift = xShift;
			double previousEnergy = monomerEnergy;
			double deltaXShift = -0.1;
			double xShiftEnd = MslTools::toDouble(parsedGeoInformation[4]);
			
			while (xShift >= xShiftEnd) {
			
				xShift += deltaXShift;

				// Move the helix
				backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, 3 );
			
				// Run Optimization
				repackSideChains(sys, spm, opt.greedy, opt.greedyCycles);
			
				vector<unsigned int> MCOFinal;
				MCOFinal = spm.getMinStates()[0];
				sys.setActiveRotamers(MCOFinal);
			
				currentEnergy = spm.getMinBound()[0];
			
				if (currentEnergy < bestEnergy) {
					bestEnergy = currentEnergy;
					savedXShift = xShift;
					sys.saveAltCoor("savedBestState");
					helicalAxis.saveAltCoor("BestAxis");
				}
			
				fout << "xShift: " << xShift << " energy: " << currentEnergy-monomerEnergy << endl;
			
				// If energy increase twice in a row, and it is above the monomer energy, quit
				if (currentEnergy > (monomerEnergy+10.0) && previousEnergy > (monomerEnergy+10.0) && currentEnergy > previousEnergy) {
					//fout << "Energy increasing above monomer energy..." << endl;
					break;	
				}
				else {
					previousEnergy = currentEnergy;
				}
			
			}
			

			/******************************************************************************
			 *               === LOCAL BACKBONE MONTE CARLO REPACKS ===
			 ******************************************************************************/
			//fout << endl << "Monte Carlo repack from best position." << endl;

      			xShift = savedXShift;

			if (opt.MCCycles > 0) {
				MonteCarloManager MCMngr(1000.0, 0.5, opt.MCCycles, MonteCarloManager::EXPONENTIAL, opt.MCMaxRejects);

				MCMngr.setEner(bestEnergy);
				
				while(!MCMngr.getComplete()) {

					sys.applySavedCoor("savedBestState");
					helicalAxis.applySavedCoor("BestAxis");

					int moveToPreform = RNG1.getRandomInt(3);

					double deltaXShift = 0.0;
					double deltaZShift = 0.0;
					double deltaCrossingAngle = 0.0;
					double deltaAxialRotation = 0.0; 

					//======================================
					//====== Z Shift (Crossing Point) ======
					//======================================
					if (moveToPreform == 0) {
						deltaZShift = getStandardNormal(RNG1) * 0.1;
						backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaZShift, moveToPreform);
					} else if (moveToPreform == 1) {
					//===========================
					//===== Axial Rotation ======
					//===========================
						deltaAxialRotation = getStandardNormal(RNG1) * 1.0;
						backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaAxialRotation, moveToPreform);
					} else if (moveToPreform == 2) {
					//==================================
					//====== Local Crossing Angle ======
					//==================================
						deltaCrossingAngle = getStandardNormal(RNG1) * 1.0;
						backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaCrossingAngle, moveToPreform);
					} else if (moveToPreform == 3) {
					//==============================================
					//====== X shift (Interhelical Distance) =======
					//==============================================
						deltaXShift = getStandardNormal(RNG1) * 0.1;
						backboneMovement(apvChainA, apvChainB, axisA, axisB, trans, deltaXShift, moveToPreform);
					}

				
					// Run Optimization
					repackSideChains(sys, spm, opt.greedy, opt.greedyCycles);
			
					vector<unsigned int> MCOFinal = spm.getMinStates()[0];
					sys.setActiveRotamers(MCOFinal);
					currentEnergy = spm.getMinBound()[0];

					if (!MCMngr.accept(currentEnergy)) {
						//fout << "state rejected   energy: " << currentEnergy << endl;
					}
					else {
						// check number of hydrogen bonds and if it is less than 4 dont accept
						// TODO: hardcoded 4 
						vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);
						if(interHelicalHbonds.size() < 4) {
							// dont accept if we lost hydrogen bonds
						} else {
							bestEnergy = currentEnergy;
							sys.saveAltCoor("savedBestState");
							helicalAxis.saveAltCoor("BestAxis");

							xShift = xShift + deltaXShift;
							crossingAngle = crossingAngle + deltaCrossingAngle;
							axialRotation = axialRotation + deltaAxialRotation;
							zShift = zShift +  deltaZShift;

							fout << "MCAccept   xShift: " << xShift << " crossingAngle: " << crossingAngle << " axialRotation: " << axialRotation << " zShift: " << zShift << " energy: " << currentEnergy-monomerEnergy << endl;
						}
					}

				}
				//fout << "Monte Carlo repack complete. " << endl << endl;
			}

			sys.applySavedCoor("savedBestState");


			double finalEnergy = sys.calcEnergy()-monomerEnergy;
			if(opt.printTermEnergies) {
				energyByTerm = getEnergyByTerm(sys.getEnergySet());

				for(map<string,double>::iterator it = monomerEnergyByTerm.begin(); it != monomerEnergyByTerm.end(); it++) {
					if(energyByTerm.find(it->first) != energyByTerm.end()) {
						energyByTerm[it->first] -= monomerEnergyByTerm[it->first];
					} else {
						// impossible
						fout << "ERROR 1242 Missing Energy term " << it->first << " Exiting" << endl; 
						exit(0);
					}
				}
			}
			if (finalEnergy > 0) { // TODO Check number of Hbonds???
				// renumber residues - rethreading
				fout << "Energy " << finalEnergy << " above zero - not written" << endl;
				fout << endl;
				//reThreadResidues(positions,-1);

				continue;
			}

			
			moveZCenterOfCAMassToOrigin(sys.getAtomPointers());

			// Renumber structures (all should go 1 - N, or from starting point determined by user)
			//int currStartNum = positions[0]->getResidueNumber();
			chainA.renumberChain(opt.startResNum);
			chainB.renumberChain(opt.startResNum);
			//renumberResidues(sys,opt.startResNum);

			// get vector of Hbonds, check size
			vector<string> interHelicalHbonds =  getInterHelicalHbonds(Eset);
			//fout << "number of CAH bonds: " << interHelicalHbonds.size() << endl; 

			// Store everything
			// create class structure to store data
			HelixDimer * st = new HelixDimer("",finalEnergy,j);
			
			st->addAtoms(sys.getAtomPointers());
			// save the structure interface residues, energy, hbond list, CA of closest approach mark it with 2 in the mask 
			map<string, unsigned int> interfaceMap = interfaceResidueCheck(apvChainA, apvChainB);

			unsigned int closestCA = CAOfClosestApproach(chainA,chainB);
			string interfaceString = "";
			for (uint r=0; r < chainA.positionSize(); r++) {
				sequence = sequence + MslTools::getOneLetterCode(chainA.getPosition(r).getResidueName());

				if (interfaceMap.find(chainA.getIdentity(r).getIdentityId()) != interfaceMap.end()) {
					interfaceString = interfaceString + "1";
				}
				else {
					interfaceString = interfaceString + "0";
				}
			}

			//mark residue with CA of closest approach with 2
			interfaceString.replace(closestCA - opt.startResNum,1,"2");

			st->setHelixDimerDetails(xShift,zShift,axialRotation,crossingAngle,interfaceString,prolineMask,interHelicalHbonds);
			if(opt.printTermEnergies) {
				st->setDeltaEnergyByTerm(energyByTerm);
			}

			if (opt.clusterSolutions) {
				// Save structure to vector of structures
				fout << "clustering structure with Energy: " << finalEnergy << endl;
				fout << endl;
				structures.push_back(st);
				
				AtomPointerVector& axisAtoms = helicalAxis.getAtomPointers();
				for(AtomPointerVector::iterator it =  axisAtoms.begin(); it != axisAtoms.end(); it++) {
					Atom* a = new Atom(* *it);
					st->addAxisAtom(a);
				}
				
			}
			else { 
				char tmp[1000];


				sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.crd", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
				CRDWriter crd;
				string crdName = opt.pdbOutputDir+"/"+string(tmp);
				crd.open(crdName);
				crd.write(sys.getAtomPointers());
				fout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
				crd.close();

				// TODO convert to PDB names before printing
				sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.pdb", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
				fout << "wrote file: " << string(tmp) << " Energy: " << finalEnergy << endl;
				sys.writePdb(opt.pdbOutputDir+"/"+string(tmp));

				sprintf(tmp,"%s_A_%02d_B_%02d_x%05.3f_z%05.3f_a%06.3f_c%05.3f.txt", opt.uniprotAccession.c_str(), positions[0]->getResidueNumber(), positions[positions.size()/2]->getResidueNumber(), xShift, zShift, axialRotation, crossingAngle);
				fout << "wrote text file: " << string(tmp) << endl;
				
				ofstream txt;
				txt.open(tmp);
				if(txt.is_open()) {
					txt << "resStart " << opt.startResNum << endl;
					txt << "resEnd " << opt.endResNum << endl;
					txt << "TMStart " << opt.startResNum << endl;
					txt << "TMEnd " << opt.endResNum << endl;
					txt << "originalSeq " << originalTMSeq << endl;
					txt << "modelledSeq " << modelledTMSeq << endl;
					st->printHelixDimerDetails(txt);
					delete st;

				} else {
					cerr << "Unable to open " << tmp << endl;
				}
				

				fout << sys.getEnergySummary();
				fout << endl;
			}
		} // Geo Loop End
	} // Threading Loop End

	if (opt.clusterSolutions) {
		sort(structures.begin(),structures.end(),compareHelixDimers());
		// convert all structures to PDB names
		/*
		for(int i = 0; i < structures.size(); i++) {
			How do we handle changed hydrogen bond names?
			AtomPointerVector & a = structures[i]->getAtomPointers();
			for(int j = 0; j < a.size(); j++) {
				int resNum = a[j]->getResidueNumber();
				fc.setPdbFromCharmm(*a[j],"22",resNum == opt.startResNum,resNum == opt.endResNum);
			}
		}
		*/
		vector<HelixDimerCluster*> clusters;
		clusterSolutions(clusters,structures,opt.rmsdCutoff,originalTMSeq,modelledTMSeq,opt);
		//cout << "Output " << clusters.size() << " HelixDimerClusters at rmsd " << _rmsdCutoff <<  endl;
		for(int i = 1; i <= clusters.size(); i++) {
			clusters[i-1]->printHelixDimerClusterCrds(i,opt.printAllCrds,opt.printAxes);
			clusters[i-1]->printDetails(i);
			if(opt.printTermEnergies) {
				clusters[i-1]->printTermEnergies(i);
			}
		}


		for(int i = 1; i <= clusters.size(); i++) {
			//clusters[i-1]->convertToPdbNames();
			clusters[i-1]->printHelixDimerClusterPdbs(i);
			clusters[i-1]->makePse(i);
		}

		// TODO delete HelixDimerCluster* form the vector clusters

	}

	time(&endTime);
	diffTime = difftime (endTime, startTime);
	fout << endl << "Total Time: " << setiosflags(ios::fixed) << setprecision(0) << diffTime << " seconds" << endl;
}


Options parseOptions(int _argc, char * _argv[], Options defaults) {

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	vector<string> allowed;

	opt.required.push_back("fullSequence");

	opt.allowed.push_back("backboneCrd");
	opt.allowed.push_back("logFile");
	opt.allowed.push_back("pdbOutputDir");
	opt.allowed.push_back("tmStart");
	opt.allowed.push_back("tmEnd");

	opt.allowed.push_back("helixGeoFile");
	opt.allowed.push_back("rulesFile");

	opt.allowed.push_back("topFile");
	opt.allowed.push_back("parFile");
	opt.allowed.push_back("solvFile");
	opt.allowed.push_back("hBondFile");
	opt.allowed.push_back("rotLibFile");

	opt.allowed.push_back("MCCycles");
	opt.allowed.push_back("MCMaxRejects");

	opt.allowed.push_back("verbose");
	opt.allowed.push_back("greedyOptimizer");
	opt.allowed.push_back("greedyCycles");
	opt.allowed.push_back("seed");
	
	opt.allowed.push_back("uniprotName");
	opt.allowed.push_back("uniprotAccession");

	opt.allowed.push_back("rmsdCutoff");
	opt.allowed.push_back("clusterSolutions");
	opt.allowed.push_back("printAllCrds");
	opt.allowed.push_back("printAxes");
	opt.allowed.push_back("printTermEnergies");

	opt.allowed.push_back("fullSequenceStart");

	opt.allowed.push_back("startResNum");
	opt.allowed.push_back("endResNum");

	opt.allowed.push_back("threadStart");
	opt.allowed.push_back("threadEnd");

	opt.allowed.push_back("configfile");

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
	if (OP.fail()) {
		opt.help = OP.getBool("h");
	}

	if (opt.help) {
		help(defaults);
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;

	opt.errorMessages = "";
	opt.warningMessages = "";

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.fullSequence = OP.getString("fullSequence");
	if (OP.fail()) {
		opt.errorMessages += "fullSequence (1 letter aa) not specified\n";
		opt.errorFlag = true;
	}

	opt.backboneCrd = OP.getString("backboneCrd");
	if (OP.fail()) {
		opt.warningMessages += "backboneCrd file not specified using /data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd\n";
		opt.warningFlag = true;
		opt.backboneCrd = "/data01/sabs/tmRepacks/pdbFiles/69-gly-residue-helix.crd";
	}
	opt.helixGeoFile = OP.getString("helixGeoFile");
	if (OP.fail()) {
		opt.helixGeoFile = "/data01/sabs/tmRepacks/centroids_0.5_2/all_centroids.txt";
		opt.warningMessages += "helixGeoFile file not specified using " + opt.helixGeoFile + "\n";
		opt.warningFlag = true;
	}

	opt.uniprotName = OP.getString("uniprotName");
	if (OP.fail()) {
		opt.uniprotName = "PROTEIN_UNK";
		opt.warningMessages += "uniprotName not specified using " + opt.uniprotName + "\n";
		opt.warningFlag = true;
	}
	opt.uniprotAccession = OP.getString("uniprotAccession");
	if (OP.fail()) {
		opt.uniprotAccession = "P00000";
		opt.warningMessages += "uniprotAccession not specified using " + opt.uniprotAccession + "\n";
		opt.warningFlag = true;
	}

	opt.pdbOutputDir = OP.getString("pdbOutputDir");
	if (OP.fail()) {
		opt.pdbOutputDir = "CAHTM-" + MslTools::getRandomAlphaNumString(6);
		opt.warningMessages += "pdbOutputDir not specified using " + opt.pdbOutputDir + "\n";
		string cmd = "mkdir -p " +  opt.pdbOutputDir;
		if(system(cmd.c_str()) ) {
			cerr << "Unable to make directory " << opt.pdbOutputDir << endl;
			exit(0);
		}
		opt.warningFlag = true;
	}

	opt.logFile = OP.getString("logFile");
	if (OP.fail()) {
		opt.logFile = opt.pdbOutputDir + "/" + opt.uniprotAccession + ".log";
		opt.warningMessages += "logFile not specified using " + opt.logFile + "\n";
		opt.warningFlag = true;
	}

	
	opt.tmStart = OP.getInt("tmStart");
	if(OP.fail()) {
		opt.warningMessages += "tmStart not specified using 1\n";
		opt.warningFlag = true;
		opt.tmStart = 1;
	}
	opt.tmEnd = OP.getInt("tmEnd");
	if(OP.fail()) {
		opt.tmEnd = opt.fullSequence.length();
		opt.warningMessages += "tmEnd not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
	}

	opt.rulesFile = OP.getString("rulesFile");
	if (OP.fail()) {
		opt.rulesFile = "/data01/sabs/tmRepacks/GLY_69_Homo_2/tmRules/rules_10kcals_vdw_only/tmRules.out";
		opt.warningMessages += "rulesFile not specified using " + opt.rulesFile + "\n";
		opt.warningFlag = true;
	}


	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_TOP";
		if(ENV.isDefined(envVar)) {
			opt.topFile = ENV.getEnv(envVar);
			opt.warningMessages += "topFile not specified using " + opt.topFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine topFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}

	}
	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_PAR";
		if(ENV.isDefined(envVar)) {
			opt.parFile = ENV.getEnv(envVar);
			opt.warningMessages += "parFile not specified using " + opt.parFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine parFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		string envVar = "MSL_CHARMM_SOLV";
		if(ENV.isDefined(envVar)) {
			opt.solvFile = ENV.getEnv(envVar);
			opt.warningMessages += "solvFile not specified using " + opt.solvFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine solvFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		string envVar = "MSL_HBOND_CA_PAR";
		if(ENV.isDefined(envVar)) {
			opt.hBondFile = ENV.getEnv(envVar);
			opt.warningMessages += "hBondFile not specified using " + opt.hBondFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine hBondFile - MSL_HBOND_CA_PAR - not set\n"	;
			opt.errorFlag = true;
		}
	}
	opt.rotLibFile = OP.getString("rotLibFile");
	if (OP.fail()) {
		string envVar = "MSL_EBL";
		if(ENV.isDefined(envVar)) {
			opt.rotLibFile = ENV.getEnv(envVar);
			opt.warningMessages += "rotLibFile not specified using " + opt.rotLibFile + "\n";
			opt.warningFlag = true;
		} else {
			opt.errorMessages += "Unable to determine rotLibFile - " + envVar + " - not set\n"	;
			opt.errorFlag = true;
		}
	}

	opt.MCCycles = OP.getInt("MCCycles");
	if (OP.fail()) {
		opt.warningMessages += "MCCycles not specified using 10\n";
		opt.warningFlag = true;
		opt.MCCycles = 10;
	}
	opt.MCMaxRejects = OP.getInt("MCMaxRejects");
	if (OP.fail()) {
		opt.warningMessages += "MCMaxRejects not specified using 2\n";
		opt.warningFlag = true;
		opt.MCMaxRejects = 2;
	}

	opt.verbose = OP.getBool("verbose");
	if (OP.fail()) {
		opt.warningMessages += "verbose not specified using false\n";
		opt.warningFlag = true;
		opt.verbose = false;
	}
	opt.greedy = OP.getBool("greedyOptimizer");
	if (OP.fail()) {
		opt.warningMessages += "greedyOptimizer not specified using true\n";
		opt.warningFlag = true;
		opt.greedy = true;
	}
	opt.greedyCycles = OP.getInt("greedyCycles");
	if (OP.fail()) {
		opt.warningMessages += "greedyCycles not specified using 1\n";
		opt.warningFlag = true;
		opt.greedyCycles = 1;
	}
	opt.seed = OP.getInt("seed");
	if (OP.fail()) {
		opt.warningMessages += "seed not specified using 1\n";
		opt.warningFlag = true;
		opt.seed = 1;
	}

	
	opt.clusterSolutions = OP.getBool("clusterSolutions");
	if (OP.fail()) {
		opt.clusterSolutions = true;
		opt.warningMessages += "clusterSolutions not specified using true\n";
		opt.warningFlag = true;
	}
	opt.printAllCrds = OP.getBool("printAllCrds");
	if (OP.fail()) {
		opt.warningMessages += "printAllCrds not specified using false\n";
		opt.warningFlag = true;
		opt.printAllCrds = false;
	}
	opt.printAxes = OP.getBool("printAxes");
	if (OP.fail()) {
		opt.warningMessages += "printAxes not specified using false\n";
		opt.warningFlag = true;
		opt.printAxes = false;
	}
	opt.printTermEnergies = OP.getBool("printTermEnergies");
	if (OP.fail()) {
		opt.printTermEnergies = true;
		opt.warningMessages += "printTermEnergies not specified using true\n";
		opt.warningFlag = true;
	}
	opt.rmsdCutoff = OP.getDouble("rmsdCutoff");
	if (OP.fail()) {
		opt.rmsdCutoff = 2.0;
		opt.warningMessages += "rmsdCutoff not specified using 2.0\n";
		opt.warningFlag = true;
	}

	opt.fullSequenceStart = OP.getInt("fullSequenceStart");
	if (OP.fail()) {
		opt.warningMessages += "fullSequenceStart not specified using 1\n";
		opt.warningFlag = true;
		opt.fullSequenceStart = 1;
	}

	opt.startResNum = OP.getInt("startResNum");
	if (OP.fail()) {
		opt.warningMessages += "startResNum not specified using " + MslTools::intToString(opt.tmStart) + "\n";
		opt.warningFlag = true;
		opt.startResNum = opt.tmStart;
	}
	opt.endResNum = OP.getInt("endResNum");
	if (OP.fail()) {
		opt.warningMessages += "endResNum not specified using " + MslTools::intToString(opt.tmEnd) + "\n";
		opt.warningFlag = true;
		opt.endResNum = opt.tmEnd;
	}

	opt.threadStart = OP.getInt("threadStart");
	if (OP.fail()) {
		// excluded 2 residues in the beginning
		// opt.threadStart = 35 - ( (tmEnd-tmStart + 1) - 3 );  
		 opt.threadStart = opt.startResNum + 37 - opt.endResNum ;  
		opt.warningMessages += "threadStart not specified using " + MslTools::intToString(opt.threadStart) + "\n";
		opt.warningFlag = true;
	}

	opt.threadEnd = OP.getInt("threadEnd");
	if (OP.fail()) {
		opt.threadEnd = 33;
		opt.warningMessages += "threadEnd not specified using " + MslTools::intToString(opt.threadEnd) + "\n";
		opt.warningFlag = true;
	}

	// return the Options structure

	return opt;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help(Options defaults) {
	cout << "This program runs as:" << endl;
	cout << " % CAHTMPredictionProgram " << endl;
	cout << "   --fullSequence <1-letter sequence> " << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --backboneCrd <backboneCrd> --logFile <logFile> --pdbOutputDir <dir> " << endl;
	cout << "   --tmStart <int> --tmEnd <int> --helixGeoFile <geometryFile> " << endl;
	cout << "   --rulesFile <file> --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file> --rotLibFile <file>" << endl;
	cout << "   --MCCycles <int> --MCMaxRejects=<int>" << endl;
	cout << "   --greedyOptimizer=<true/false> --greedyCycles=<int>  --seed <int> --verbose <true/false>" << endl;
	cout << "   --uniprotName <name> --uniportAccession <name> " << endl;
	cout << "   --clusterSolutions=<true/false> --rmsdCutoff=<double> --printAllCrds<true/false> --printAxes <true/false> --printTermEnergies <true/false>" << endl;
	cout << "   --fullSequenceStart <int> " << endl;
	cout << "   --startResNum <int> --endResNum <int>" << endl;
	cout << "   --threadStart <int> --threadEnd <int>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << endl;
}

