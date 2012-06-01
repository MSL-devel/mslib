#include "MslTools.h"
#include "SystemRotamerLoader.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"
#include "OptionParser.h"
#include <set>
#include <vector>

using namespace MSL;
using namespace std;
string programName = "makeHonigLibrary";
string programDescription = "This program creates a PDB file from a rotamerlibrary file by building the conformations on a GLY-X-GLY backbone";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "OCT 5 2011";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct Options {
	string commandName; // name of this program
	string rotlibFile; // file containing the atom coordinates in PDB format
	string charmmTopFile;
	string charmmParFile;

	string outputfile; // logFile
	string configfile; // configFile

	map<string,int> numRots; // number of rotamers for each residue type 

	/***** MANAGEMENT VARIABLES ******/
	bool version; // ask for the program version
	bool help; // ask for the program help


	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> required; //list of required options
	vector<string> allowed; //list of allowed options

	string OPerrors; //the errors from the option parser
};
void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " --help" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
       	cout << "Run  as:" << endl;
	cout << " % makeHonigLibrary \n --rotlibfile <rotlibfile> \n --charmmtopfile <charmmTopFile> \n --charmmparfile <charmmParFile>  --outputpdbfile <outputpdbfile>" << endl;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << " --configfile <configfile> \n"<< endl;
	cout << endl;
	cout << "Optional - Num of rotamers per residue type defaults to all the conformations in the library" << endl;
	cout << " --ALA <nALA>\n --ARG <nARG> \n .........\n --HSD <nHSD>\n  ....... \n --VAL <nVAL> " << endl;
	cout << endl;
}
Options parseOptions(int _argc, char * _argv[]) {

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
	opt.required.push_back("rotlibfile"); // rotamerLibrary 
	opt.required.push_back("charmmtopfile"); // rotamerLibrary 
	opt.required.push_back("charmmparfile"); // rotamerLibrary 
	opt.required.push_back("outputfile"); // repacked structure will be written to this file

	opt.allowed.push_back("outputfile"); // all output will be redirected to this logFile
	opt.allowed.push_back("configfile");

	// To specify the number of rotamers for each amino acid type
	opt.allowed.push_back("ALA"); 
	opt.allowed.push_back("ARG");
	opt.allowed.push_back("ASN");
	opt.allowed.push_back("ASP");
	opt.allowed.push_back("CYS");
	opt.allowed.push_back("GLN");
	opt.allowed.push_back("GLU");
	opt.allowed.push_back("GLY");
	opt.allowed.push_back("HSD");
	opt.allowed.push_back("ILE");
	opt.allowed.push_back("LEU");
	opt.allowed.push_back("LYS");
	opt.allowed.push_back("MET");
	opt.allowed.push_back("PHE");
	opt.allowed.push_back("PRO");
	opt.allowed.push_back("SER");
	opt.allowed.push_back("THR");
	opt.allowed.push_back("TRP");
	opt.allowed.push_back("TYR");
	opt.allowed.push_back("VAL");

	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("help"); // --help

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.readArgv(_argc, _argv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version prints the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");

	if (opt.help) {
		help();
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	opt.errorFlag = false;
	opt.warningFlag = false;
	opt.errorMessages = "";
	opt.warningMessages = "";

	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			return opt;
		}
	}
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}

	opt.outputfile = OP.getString("outputfile");
	if (!OP.fail()) {
		if(!freopen(opt.outputfile.c_str(),"w",stdout)) {
			opt.errorFlag = true;
			opt.errorMessages +=  "Unable to redirect output to " + opt.outputfile + "\n";
			return opt;
		}
	}

	// Print Configuration File / Commmand Line Options
	//OP.printConfFile();

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	
	opt.rotlibFile = OP.getString("rotlibfile");
	if (OP.fail()) {
		opt.errorMessages = "rotlibfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.charmmTopFile = OP.getString("charmmtopfile");
	if (OP.fail()) {
		opt.errorMessages = "charmmtopfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.charmmParFile = OP.getString("charmmparfile");
	if (OP.fail()) {
		opt.errorMessages = "charmmparfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}

	int num = OP.getInt("ALA"); 
	if(!OP.fail()) {
		opt.numRots["ALA"] = num;
	}
	num = OP.getInt("ARG");
	if(!OP.fail()) {
		opt.numRots["ARG"] = num;
	}
	num = OP.getInt("ASN");
	if(!OP.fail()) {
		opt.numRots["ASN"] = num;
	}
	num = OP.getInt("ASP");
	if(!OP.fail()) {
		opt.numRots["ASP"] = num;
	}
	num = OP.getInt("CYS");
	if(!OP.fail()) {
		opt.numRots["CYS"] = num;
	}
	num = OP.getInt("GLN");
	if(!OP.fail()) {
		opt.numRots["GLN"] = num;
	}
	num = OP.getInt("GLU");
	if(!OP.fail()) {
		opt.numRots["GLU"] = num;
	}
	num = OP.getInt("GLY");
	if(!OP.fail()) {
		opt.numRots["GLY"] = num;
	}
	num = OP.getInt("HSD");
	if(!OP.fail()) {
		opt.numRots["HSD"] = num;
	}
	num = OP.getInt("ILE");
	if(!OP.fail()) {
		opt.numRots["ILE"] = num;
	}
	num = OP.getInt("LEU");
	if(!OP.fail()) {
		opt.numRots["LEU"] = num;
	}
	num = OP.getInt("LYS");
	if(!OP.fail()) {
		opt.numRots["LYS"] = num;
	}
	num = OP.getInt("MET");
	if(!OP.fail()) {
		opt.numRots["MET"] = num;
	}
	num = OP.getInt("PHE");
	if(!OP.fail()) {
		opt.numRots["PHE"] = num;
	}
	num = OP.getInt("PRO");
	if(!OP.fail()) {
		opt.numRots["PRO"] = num;
	}
	num = OP.getInt("SER");
	if(!OP.fail()) {
		opt.numRots["SER"] = num;
	}
	num = OP.getInt("THR");
	if(!OP.fail()) {
		opt.numRots["THR"] = num;
	}
	num = OP.getInt("TRP");
	if(!OP.fail()) {
		opt.numRots["TRP"] = num;
	}
	num = OP.getInt("TYR");
	if(!OP.fail()) {
		opt.numRots["TYR"] = num;
	}
	num = OP.getInt("VAL");
	if(!OP.fail()) {
		opt.numRots["VAL"] = num;
	}
	
	// return the Options structure
	return opt;

}


void printAtomLine (string _atomName, string _chainId, double x, double y, double z) {
	char line[1000];
	sprintf(line,"%4s %1s %8.3f %8.3f %8.3f",_atomName.c_str(),_chainId.c_str(),x,y,z);
	cout << line << endl;
}


int main(int argc, char* argv[]) {

	Options opt = parseOptions(argc, argv);
	if(opt.errorFlag) {
		cout << opt.OPerrors ;
		cout << opt.errorMessages ;
		exit(0);
	}

	if(opt.warningFlag) {
		cout << opt.warningMessages;
	}
	

	map<string,string> oneLetterCode;
	oneLetterCode["ALA"] = "A"; 
	oneLetterCode["ARG"] = "R";
	oneLetterCode["ASN"] = "N";
	oneLetterCode["ASP"] = "D";
	oneLetterCode["CYS"] = "C";
	oneLetterCode["GLN"] = "Q";
	oneLetterCode["GLU"] = "E";
	oneLetterCode["GLY"] = "G";
	oneLetterCode["HSD"] = "H";
	oneLetterCode["ILE"] = "I";
	oneLetterCode["LEU"] = "L";
	oneLetterCode["LYS"] = "K";
	oneLetterCode["MET"] = "M";
	oneLetterCode["PHE"] = "F";
	oneLetterCode["PRO"] = "P";
	oneLetterCode["SER"] = "S";
	oneLetterCode["THR"] = "T";
	oneLetterCode["TRP"] = "W";
	oneLetterCode["TYR"] = "Y";
	oneLetterCode["VAL"] = "V";




	System sys;
	SystemRotamerLoader rotLoader;
	rotLoader.setSystem(sys);

	if(!rotLoader.readRotamerLibraryFile(opt.rotlibFile)) {
		cerr << "Unable to read " << opt.rotlibFile << endl;
		exit(0);
	}

	RotamerLibrary *rotlib = rotLoader.getRotamerLibrary();

	string seq = "" ;
	set<string> resList = rotlib->getAllResList();

	for(set<string>::iterator it = resList.begin(); it != resList.end(); it++) {
		if(oneLetterCode.find(*it) != oneLetterCode.end() && opt.numRots.find(*it) !=  opt.numRots.end() ) {
			seq += oneLetterCode[*it] + ":" + "GLY " + *it + " GLY \n";
		}
	}
	
	//cout << seq << endl;

	CharmmSystemBuilder csb(sys,opt.charmmTopFile,opt.charmmParFile);
	if(!csb.buildSystem(seq)) {
		cerr << "Unable to build System" << endl;
		exit(0);
	}
	//cout << sys << endl;
	vector<Chain*> chains = sys.getChains();
	for(int i = 0; i < chains.size(); i++) {
		string chainId = chains[i]->getChainId();
		if(!sys.seed(chainId + ",1,C", chainId + ",1,CA",chainId + ",1,N")) {
			cerr << "Unable to seed " << endl;
			exit(0);
		}
	}
	sys.buildAllAtoms();
	//cout << sys.getAtomPointers() << endl;
		
	for(int i = 0; i < chains.size(); i++) {	
		string chainId = chains[i]->getChainId();
		Position& pos = chains[i]->getPosition(1);
		string resName = pos.getCurrentIdentity().getResidueName();
		int tempSize;
		if(opt.numRots.find(resName) == opt.numRots.end()) {
			continue;
		}
		tempSize = opt.numRots[resName];
		if(rotlib->size("",resName) < opt.numRots[resName]) {
			tempSize = rotlib->size("",resName);
		}
		if(tempSize == 0) {
			cerr << "No conformers in the library for " << resName << endl;
			continue;
		}
		if(!rotLoader.loadRotamers(&pos,resName,tempSize)) {
			cerr << "Unable to load rotamers " << resName << endl;
			exit(0);
		}
		if(!chains[i]->atomExists("1,CA")) {
			continue;
		}
		Atom& CA = chains[i]->getLastFoundAtom();
		if(!chains[i]->atomExists("1,C")) {
			continue;
		}
		Atom& C = chains[i]->getLastFoundAtom();
		if(!chains[i]->atomExists("1,O")) {
			continue;
		}
		Atom& O = chains[i]->getLastFoundAtom();
		//cout << resName << " libSize" << rotlib->size("",resName) <<  " loaded " << tempSize << endl;
		for(int j = 0; j < tempSize; j++) {
			(pos.getCurrentIdentity()).setActiveConformation(j);
			AtomPointerVector& atoms = pos.getAtomPointers();
			for(int k = 0; k < atoms.size(); k++) {
				printAtomLine(atoms[k]->getName(),chainId,atoms[k]->getX(),atoms[k]->getY(),atoms[k]->getZ());
			}
			printAtomLine("HEAD","1",CA.getX(),CA.getY(),CA.getZ());
			printAtomLine("HEAD","2",C.getX(),C.getY(),C.getZ());
			printAtomLine("HEAD","3",O.getX(),O.getY(),O.getZ());
			//cout << atoms << endl;
		}
	}
}
