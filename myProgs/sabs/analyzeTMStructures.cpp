#include <iostream>
#include <fstream>

#include "CharmmSystemBuilder.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "MslTools.h"
#include "HydrogenBondBuilder.h"


using namespace MSL;
using namespace std;

string programName = "analyzeTMStructures";
string programDescription = "This program counts the number of backbone-backbone/backbone-sidechain interhelical hydrogen bonds in TM segments";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "0.0.1";
string programDate = "23 April 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

time_t startTime, endTime;
double diffTime;

class Segment;
map<string,string> posIdToTmIdMap;
vector<Segment*> segments;
map<string,Segment*> mappedSegments;

struct Options {
	// Required
	string pdbFile;
	vector<string> segment;
	string proteinId;

	//required
	string topFile;
	string parFile;
	string solvFile;
	string hBondFile;

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

Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);

string getSegment(Atom& a);

class Segment {

	public:
	Segment(System& _sys, string _chain, string _tmId, int _start, int _end ) {
		AtomSelection sel(_sys.getAtomPointers());
		atoms = sel.select("seg,chain " + _chain + " and resi " + MslTools::intToString(_start) + "-" + MslTools::intToString(_end));
		caAtoms = sel.select("segCA,chain " + _chain + " and resi " + MslTools::intToString(_start) + "-" + MslTools::intToString(_end) + " and name CA");
		sequence ="";
		for(int i = 0; i < caAtoms.size(); i++) {
			sequence += MslTools::getOneLetterCode(caAtoms[i]->getResidueName());
		}

		tmId = _tmId;
		chain = _chain;
		start = _start;
		end = _end;
	}

	bool exists (Atom& at) {
		if(at.getChainId() == chain) {
			int resNum = at.getResidueNumber();
			if(resNum >= start && resNum <= end) {
				return true;
			}
		}
		return false;
	}

	string getTmId() {
		return tmId;
	}

	string getChain() {
		return chain;
	}
	AtomPointerVector& getAtoms() {
		return atoms;
	}

	AtomPointerVector& getCAAtoms() {
		return caAtoms;
	}

	string getSequence() {
		return sequence;
	}

	void print() {
		cout << "TMSEGMENT " << chain << " " << tmId << " " << start << " " << end << " " << sequence << endl;
	}
	// returns residueNumbers of CAs of closest approach
	void CAofClosestApproach(Segment _other, int& thisSeg, int& otherSeg) {
		AtomPointerVector& otherCA = _other.getCAAtoms();
		int bestI = 0;
		int bestJ = 0;
		double bestDist = MslTools::floatMax;
		for(int i = 0; i < caAtoms.size(); i++) {
			for(int j = 0; j < otherCA.size(); j++) {
				double dist = caAtoms[i]->distance(*otherCA[j]);
				if(dist < bestDist) {
					bestI = i;
					bestJ = j;
					bestDist = dist;
				}
			}
		}
		thisSeg = caAtoms[bestI]->getResidueNumber();
		otherSeg = otherCA[bestJ]->getResidueNumber();
	}
	
	private:
	AtomPointerVector atoms,caAtoms;
	string chain;
	string tmId;
	string sequence;
	int start;
	int end;
};

string getSegment(Atom& a) {
	string posId = a.getPositionId();
	if(posIdToTmIdMap.find(posId) == posIdToTmIdMap.end()) {
		for(int i = 0; i < segments.size(); i++) {
			if(segments[i]->exists(a)) {
				posIdToTmIdMap[posId] =  segments[i]->getTmId();
				return posIdToTmIdMap[posId];
			}
		}
	} else {
		return posIdToTmIdMap[posId];
	}
	//cerr << "Segment not found for " << posId << endl;
	return "NONE";
}

map<string, unsigned int> interfaceResidueCheck(AtomPointerVector & _segA, AtomPointerVector & _segB) {
	map<string, unsigned int> atomsWithIn4AOfPosition;
	for (uint i=0; i < _segA.size(); i++) {
		string segAName = _segA[i]->getName();
		if (segAName != "CA" && segAName != "C" && segAName != "N" && segAName != "HN" && segAName != "O") {
			for (uint j=0; j < _segB.size(); j++) {
				string segBName = _segB[j]->getName();
				if (segBName != "CA" && segBName != "C" && segBName != "N" && segBName != "HN" && segBName != "O") {
					if (_segA[i]->distance(*_segB[j]) < 4.0) {
						if (atomsWithIn4AOfPosition.find(_segA[i]->getIdentityId()) != atomsWithIn4AOfPosition.end()) {
							atomsWithIn4AOfPosition[_segA[i]->getIdentityId()]++;
						}
						else {
							atomsWithIn4AOfPosition[_segA[i]->getIdentityId()] = 1;
						}
						if (atomsWithIn4AOfPosition.find(_segB[j]->getIdentityId()) != atomsWithIn4AOfPosition.end()) {
							atomsWithIn4AOfPosition[_segB[j]->getIdentityId()]++;
						}
						else {
							atomsWithIn4AOfPosition[_segB[j]->getIdentityId()] = 1;
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

struct Hbond {

	Hbond(string _donorTmId, Atom* _donor, string _acceptorTmId, Atom* _acceptor) {
		donorTmId = _donorTmId;
		acceptorTmId = _acceptorTmId;
		donor = _donor;
		acceptor = _acceptor;
		distance = donor->distance(*acceptor);
	}

	string getFormattedLine() {
		// hbonds 1 2 B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832
	
		char bond[1000];
		sprintf(bond,"%s %s %s:%s=%4.3f",donorTmId.c_str(),acceptorTmId.c_str(),donor->getAtomOfIdentityId().c_str(),acceptor->getAtomOfIdentityId().c_str(),distance);
		return string(bond);

	}

	Atom* getDonor() { return donor;}
	string getDonorTmId() {return donorTmId;}
	Atom* getAcceptor() {return acceptor;}
	string getAcceptorTmId() {return acceptorTmId;}
	double getDistance() {return distance;}

	private:
	string donorTmId;
	string acceptorTmId; 
	Atom* donor;
	Atom* acceptor;
	double distance;
};

void getInterHelicalHbonds(EnergySet* & _ESet,vector<Hbond*>& _hbonds) {
	vector<Interaction*> hbondInteractions = (*(_ESet->getEnergyTerms()))["SCWRL4_HBOND"];
	////cout << "Num interactions " << hbondInteractions.size() << endl;

	for(int i = 0; i < hbondInteractions.size(); i++) {
		vector<Atom*> atoms =  hbondInteractions[i]->getAtomPointers();
		string seg1 = getSegment(*atoms[0]);
		string seg2 = getSegment(*atoms[2]);
		//cout << seg1 << " " << seg2 << " " << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
		if(seg1 == seg2 || seg1 == "NONE" || seg2 == "NONE") {
			continue;
		}
		double e = hbondInteractions[i]->getEnergy();
		//cout << seg1 << " " << seg2 << " " << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << " " <<  e << endl;
		if(e < 0) {
			//cout << "HBOND " << atoms[0]->getAtomId() << " " << atoms[2]->getAtomId() << endl;
			_hbonds.push_back(new Hbond(seg1,atoms[0],seg2,atoms[2]));
		}
	}
}

/*
void makePse(char* fileName, string pdbFile, vector<Hbond*> interhelicalHbonds) {
	// make the .inp file

	// PyMOL script
	
	bg_color white
	load MODEL_000.pdb
	rotate X, 90
	color green, chain B and element C
	color cyan, chain B and element C
	show cartoon
	show stick, resi 28+31+32+35+36+39+40
	distance hb000, ///A/35/HA1, ///B/32/O
	color black, hb000
	distance hb001, ///B/35/HA1, ///A/32/O
	color black, hb001
	distance hb002, ///A/36/HA, ///B/35/O
	color black, hb002
	distance hb003, ///B/36/HA, ///A/35/O
	color black, hb003
	distance hb004, ///A/39/HG1, ///B/36/O
	color black, hb004
	distance hb005, ///B/39/HG1, ///A/36/O
	color black, hb005
	zoom all, 3
	set label_color, black
	save MODEL_000.pse


	// make a pse file for each dimer interaction
	ofstream pout;
	pout.open((string(fileName) + ".inp").c_str());
	if(!pout.is_open()) {
		cerr << "Unable to open " << fileName << ".inp" << endl;
		exit(0);
	}
	
	pout << "bg_color white" << endl;
	pout << "load " << pdbFile << endl;
	pout << "rotate X,90,all, 0" << endl;
	pout << "show cartoon" << endl;
	for(int i = 0; i < interhelicalHbonds.size(); i++) {
		char line[1000];
		Atom* donor = interhelicalHbonds[i]->getDonor();
		Atom* acceptor = interhelicalHbonds[i]->getAcceptor();
		sprintf(line,"distance hb%03d, ///%s/%d/%s, ///%s/%d/%s\ncolor black, hb%03d\n",i,donor->getChainId().c_str(),donor->getResidueNumber(),donor->getName().c_str(),acceptor->getChainId().c_str(),acceptor->getResidueNumber(),acceptor->getName().c_str(),i);
		pout << line;
	}

	pout << "zoom all,3" << endl;
	pout << "set label_color, black" << endl;
	pout << "save " << fileName << ".pse" << endl;
	string cmd = "pymol -cqd \"@" + string(fileName) + ".inp\"";
	system(cmd.c_str());
}
*/

string getHbondPseLine(int& hbondNum, vector<Hbond*> interhelicalHbonds) {
	stringstream ss;
	for(int i = 0; i < interhelicalHbonds.size(); i++) {
		char line[1000];
		Atom* donor = interhelicalHbonds[i]->getDonor();
		Atom* acceptor = interhelicalHbonds[i]->getAcceptor();
		sprintf(line,"distance hb%03d, ///%s/%d/%s, ///%s/%d/%s\ncolor black, hb%03d\n",hbondNum,donor->getChainId().c_str(),donor->getResidueNumber(),donor->getName().c_str(),acceptor->getChainId().c_str(),acceptor->getResidueNumber(),acceptor->getName().c_str(),hbondNum);
		hbondNum++;
		ss << line;
	}
	return ss.str();
}

void setCharge(System& _sys, CharmmTopologyReader& _topRead) {
	vector<Position*> positions = _sys.getPositions();
	for(int i =0 ; i < positions.size(); i++) {
		string resName = positions[i]->getResidueName();
		if(_topRead.residueExists(resName)) {
			CharmmTopologyResidue & topRes = _topRead.getLastFoundResidue();
			unsigned int numTopAtoms = topRes.atomSize();
			string name;
			string type;
			double charge;
			string element;
			int group;
			for(unsigned int j = 0; j < numTopAtoms; j++) {
				topRes.getTopolAtom(j,name,type,charge,element,group);
				if(positions[i]->atomExists(name)) {
					Atom& a = positions[i]->getLastFoundAtom();
					a.setCharge(charge);
				}
			}
		}
	}
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
	fout.open((opt.proteinId + ".out" ).c_str());

	if(!fout.is_open()) {
		cerr << "Unable to open " << opt.proteinId << endl;
		exit(0);
	}

	System sys;
	if(!sys.readPdb(opt.pdbFile)) {
		cerr << "Unable to read from " << opt.pdbFile << endl;
		exit(0);
	}

	// Read topology file to get Charge for each atom
	CharmmTopologyReader topRead(opt.topFile);
	if(!topRead.read()) {
		cerr << "Unable to read topology file " << opt.topFile << endl;
		exit(0);
	}

	setCharge(sys,topRead);
/*
	CharmmSystemBuilder CSB(sys,opt.topFile,opt.parFile,opt.solvFile);
	CSB.setSolvent("CHEX");
	CSB.setBuildNonBondedInteractions(false);

	if(!CSB.buildSystemFromPDB(opt.pdbFile)) {
		cerr << "Unable to build system from " << opt.pdbFile << endl;
		exit(0);
	}

	sys.buildAllAtoms();
*/
	// Add hydrogen bond term
	HydrogenBondBuilder hb(sys, opt.hBondFile);
	if(!hb.buildInteractions(10)) {
		cerr << "Unable to build hbond interactions for " << opt.pdbFile << endl;
		exit(0);
	}
	// print format for the logfile
	/*
	tmsegment A,1 233 238
	tmsegment B,2 233 238
	tmsegment C,3 233 238

	sequence 1 GKDIAFIVHGYGGLVFMDLLVRR
	sequence 2 GRDIGFIVHGYGGLVFMDLLVRR
	sequence 3 GRDITFIYHGYGGLVFMDLLVRR

	interfacial 1 00000110110010011001001
	interfacial 2 00100110110010011001001
	interfacial 3 01000110110010011001001

	allhbonds B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832;B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832;B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832

	hbonds 1 2 B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832
	hbonds 1 3 B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832
	hbonds 2 3 B,235,HA2:A,234,O=2.84814;B,238,HA1:A,234,O=3.40832;A,235,HA2:B,234,O=2.84814A,238,HA1:B,234,O=3.40832
	*/

	// create segments;
	for(int i = 0; i < opt.segment.size(); i++) {
		vector<string> toks = MslTools::tokenize(opt.segment[i]);
		segments.push_back(new Segment(sys,toks[0],toks[1],MslTools::toInt(toks[2]),MslTools::toInt(toks[3])));
		fout << "tmsegment " << toks[0] << "," << toks[1] << " " << toks[2] << " " << toks[3] << " " << segments.back()->getSequence() << endl;
		mappedSegments[toks[0] + "," + toks[1]] = segments.back();
	}

	EnergySet* Eset = sys.getEnergySet();
	vector<Hbond*> interhelicalHbonds;
	getInterHelicalHbonds(Eset,interhelicalHbonds);

	//fout << "allhbonds "; 
	//for(int i = 0; i < interhelicalHbonds.size(); i++) {
	//	fout << interhelicalHbonds[i]->getFormattedLine() << ";";
	//}
	//fout << endl;

	// map hbonds based on the helices
	map<string,vector<Hbond*> > hbondsMapped;
	for(int i = 0; i < interhelicalHbonds.size(); i++) {
		string seg1 =  interhelicalHbonds[i]->getDonor()->getChainId() + "," + interhelicalHbonds[i]->getDonorTmId() ;
		string seg2 = interhelicalHbonds[i]->getAcceptor()->getChainId() + "," + interhelicalHbonds[i]->getAcceptorTmId();
		if(hbondsMapped.find(seg1 + " " + seg2) != hbondsMapped.end()) {
			hbondsMapped[seg1 + " " + seg2].push_back(interhelicalHbonds[i]);
		} else if (hbondsMapped.find(seg2 + " " + seg1) != hbondsMapped.end()) {
			hbondsMapped[seg2 + " " + seg1].push_back(interhelicalHbonds[i]);
		} else {
			hbondsMapped[seg1 + " " + seg2].push_back(interhelicalHbonds[i]);
		}
		//cout << "HERE " << seg1 << " " << seg2 << endl;
	}


	stringstream ss;

	ss << "bg_color white" << endl;
	ss << "load " << opt.pdbFile << endl;
	ss << "rotate X,90,all, 0" << endl;
	ss << "show cartoon" << endl;



	int interfaceNum = 0;
	int hbondNum = 0;
	for(map<string,vector<Hbond*> >::iterator it = hbondsMapped.begin(); it != hbondsMapped.end(); it++) {
		//if(it->second.size() < 4) {
		//	continue;
		//}
		ss << getHbondPseLine(hbondNum,it->second );
		interfaceNum++;

		// find interface
		vector<string> toks = MslTools::tokenize(it->first);
		Segment* seg1 = NULL;
		Segment* seg2 = NULL;
		if(mappedSegments.find(toks[0]) != mappedSegments.end()) {
			seg1 = mappedSegments[toks[0]];
		} else {
			cerr << "ERROR cant find segment " << opt.pdbFile << " " << toks[0] << endl;
		}

		if(mappedSegments.find(toks[1]) != mappedSegments.end()) {
			seg2 = mappedSegments[toks[1]];
		} else {
			cerr << "ERROR cant find segment " << opt.pdbFile << " " << toks[1] << endl;
		}

		AtomPointerVector& seg1Atoms = seg1->getAtoms();
		AtomPointerVector& seg2Atoms = seg2->getAtoms();

		map<string,unsigned int> interface = interfaceResidueCheck(seg1Atoms,seg2Atoms); 
		map<string,string> perChainInterface;
		
		for(map<string,unsigned int>::iterator it1 = interface.begin(); it1 != interface.end(); it1++) {
			// split A,23,GLY A,24,GLY and map [A] = 23+24
			vector<string> toks = MslTools::tokenize(it1->first,","); 
			if(perChainInterface.find(toks[0]) == perChainInterface.end()) {
				perChainInterface[toks[0]] = toks[1];
			} else {
				perChainInterface[toks[0]] += "+" + toks[1];
			}
		}

		string interfaceId = "TM_" + seg1->getChain() + "_" + seg1->getTmId() + "_" + seg2->getChain() + "_" + seg2->getTmId() ;
		
		ss << "select " << interfaceId << " , "; 

		for(map<string,string>::iterator it1 = perChainInterface.begin(); it1 != perChainInterface.end(); it1++) {
			if(it1 == perChainInterface.begin()) {
				ss << " chain " << it1->first << " and resi " << it1->second ; 
			} else {
				ss << " or chain " << it1->first << " and resi " << it1->second ; 
			}
		}

		ss << endl;
		ss << "show stick, " << interfaceId << endl;
		
		fout << "hbonds " << it->first << "  ";
		for(vector<Hbond*>::iterator it1 = it->second.begin(); it1 != it->second.end(); it1++) {
			Atom* donor = (*it1)->getDonor();
			Atom* acceptor = (*it1)->getAcceptor();
			double distance = (*it1)->getDistance();
			char tmp[10000];
			sprintf(tmp,"%s:%s=%4.3f;",donor->getAtomOfIdentityId().c_str(),acceptor->getAtomOfIdentityId().c_str(),distance);
			fout << tmp;
		}
		fout << endl;
	}

	ss << "zoom all,3" << endl;
	ss << "set label_color, black" << endl;
	ss << "save " << opt.proteinId << ".pse" << endl;

	if(interfaceNum > 0) {
		ofstream pout;
		pout.open((opt.proteinId + ".inp").c_str());
		if(!pout.is_open()) {
			cerr << "Unable to open " << opt.proteinId << ".inp" << endl;
			exit(0);
		}

		pout << ss.str() ;

		string cmd = "pymol -cqd \"@" + opt.proteinId + ".inp\"";
		system(cmd.c_str());
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

	opt.required.push_back("pdbFile");
	opt.required.push_back("segment");
	opt.required.push_back("proteinId");

	opt.required.push_back("topFile");
	opt.required.push_back("parFile");
	opt.required.push_back("solvFile");
	opt.required.push_back("hBondFile");

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

	opt.pdbFile = OP.getString("pdbFile");
	if (OP.fail()) {
		opt.errorMessages += "pdbFile not specified\n";
		opt.errorFlag = true;
	}
	opt.segment = OP.getMultiString("segment");
	if (OP.fail()) {
		opt.errorMessages += "segment not specified\n";
		opt.errorFlag = true;
	}
	opt.proteinId = OP.getString("proteinId");
	if (OP.fail()) {
		opt.errorMessages += "proteinId not specified\n";
		opt.errorFlag = true;
	}

	opt.topFile = OP.getString("topFile");
	if (OP.fail()) {
		opt.errorMessages += "topFile not specified\n";
		opt.errorFlag = true;
	}
	opt.parFile = OP.getString("parFile");
	if (OP.fail()) {
		opt.errorMessages += "parFile not specified\n";
		opt.errorFlag = true;
	}
	opt.solvFile = OP.getString("solvFile");
	if (OP.fail()) {
		opt.errorMessages += "solvFile not specified\n";
		opt.errorFlag = true;
	}
	opt.hBondFile = OP.getString("hBondFile");
	if (OP.fail()) {
		opt.errorMessages += "hBondFile not specified\n";
		opt.errorFlag = true;
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
	cout << " % analyzeTMStructures " << endl;
	cout << "   --pdbFile <pdbFile> --segment \"A 1 22 45\" --segment \"A 2 78 90\"  --proteinId <proteinId>" << endl;
	cout << "   Optional Parameters " << endl;
	cout << "   --topFile <file> --parFile <file> --solvFile <file> --hBondFile <file>" << endl;
	cout << "   --configfile <file> " << endl;
	cout << "   The segment option takes \"chainId tmId startResNum ednResNum\"" << endl;
	cout << endl;
}

