/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include "getChiRecovery.h"

using namespace std;
using namespace MSL;

string programName = "getChiRecovery";
string programDescription = "This programs calculates the chi recovery for a list of repacked pdbs among buried residues as well as exposed residues";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "13 February 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

class ResidueChiStatistics {
	public:
	ResidueChiStatistics (){
	}

	void add(string group,vector<int> count) {
		if(tot.find(group) != tot.end()) {
			tot[group] += 1;
		} else {
			tot[group] = 1;
		}

		if(recovered.find(group) == recovered.end()) {
			recovered[group] = vector<int>(count.size(),0);
		}

		for(int i = 0; i < recovered[group].size(); i++) {
			recovered[group][i] += count[i];
		}

	}

	void printStats(string resName) {
		for(map<string,int>::iterator it = tot.begin(); it != tot.end(); it++) {
			cout << resName << " " << it->first << " " << it->second;
			vector<int> cnts = recovered[it->first];
			for(int i = 0; i < cnts.size(); i++) {
				cout << " " << cnts[i];
			}
			cout << endl;
		}
	}

	private:
	map<string,int> tot; // tot["BURIED"] = 10; tot["EXPOSED"] = 15;
	map<string,vector<int> > recovered; // recovered["BURIED"][0] = 8; recovered["EXPOSED"][0] = 4
};

class ChiRecoveryStatistics {
	public:
	ChiRecoveryStatistics(string _chiDefFile,double _tol) {
		chiStat.read(_chiDefFile);
		tolerance = _tol;
	}

	~ChiRecoveryStatistics() {
		for(map<string,ResidueChiStatistics*>::iterator it = counts.begin(); it != counts.end(); it++) {
			delete it->second;
		}
		counts.clear();
	}

	bool check360(double orig, double rec) {
		// Assume orig and rec are 0-360
		double ul = orig + tolerance; 
		double ll = orig - tolerance;
		if(ul <= 360 && ll >= 0 && rec <= ul && rec >= ll) { // works for the normal case where angles are not too close to 0 and 360
			return true;
		}

		// Cases where the angles are close to 0 or 360
		if(ul > 360 ) {
			if(rec <= (ul - 360) || rec >= ll) {
				return true;
			}
		}

		if(ll < 0) {
			if(rec >= (360 + ll) || rec <= ul) {
				return true;
			}
		}
		return false;
	}

	bool check(double orig, double rec, string resName, int chiNum) {
		// All angles are eith -180 - +180 or 0-360 . 
		// Make all the angles 0-360 
		if(orig < 0) {
			orig += 360;
		}
		if(rec < 0) {
			rec += 360;
		}
		if(check360(orig,rec)) {
			return true;
		}

		
		// There are some special cases to handle for ASP CHI2, GLU CHI3, PHE CHI2, TYR CHI2
		// ASP - OD1 and OD2 could be interchangeably named so check both dihedrals
		// GLU - OE1 and OE2 could be interchangeably named so check both dihedrals
		// PHE - CD1 and CD2 could be interchangeably named so check both dihedrals
		// TYR - CD1 and CD2 could be interchangeably named so check both dihedrals
		// The two dihedrals are ~ 180 degrees apart.
		// Remember chiNum is 0-based

		if((resName == "ASP" && chiNum == 1) || (resName == "GLU" && chiNum == 2) || (resName == "PHE" && chiNum == 1) || (resName == "TYR" && chiNum == 1)) {
			if(orig < 180) {
				orig += 180;
			} else if (orig > 180) {
				orig -= 180;
			}
			if(check360(orig,rec)) {
				return true;
			}
		}
		return false;
	}

	void analyze(Residue& orig, Residue& recovered, string group) {
		string resName = orig.getResidueName();
		vector<double> origChis = chiStat.getChis(orig);
		vector<double> recoveredChis = chiStat.getChis(recovered);
		if(origChis.size() == recoveredChis.size() && origChis.size() > 0) {
			vector<int> count(origChis.size(),0);
			if(check(recoveredChis[0],origChis[0],resName,0) ) {
				count[0] = 1;
			
				for(int i = 1; i < origChis.size(); i++) {
					if(check(recoveredChis[i],origChis[i],resName,i) ) {
						count[i] = 1;
					} else {
						break;
					}
				}
			}
			if(counts.find(resName) == counts.end()) {
				counts[resName] = new ResidueChiStatistics;
			}
			counts[resName]->add(group,count);
		}
	}

	void printStats() {
		cout << "#Format: " << endl;
		cout << "#RES BURIED/EXPOSED/UNKNOWN TOTAL CHI1 CHI1+2 CHI1+2+3 CHI1+2+3+4" << endl;
		for(map<string,ResidueChiStatistics*>::iterator it = counts.begin(); it != counts.end(); it++) {
			it->second->printStats(it->first);
		}
	}

	private:
	ChiStatistics chiStat; // stores chi definitions
	map<string,ResidueChiStatistics*> counts;
	double tolerance;
};


void computeRecovery(ChiRecoveryStatistics& cRS, string orig, string repacked,map<string,double>& refSasa ) {
	System sys1;
	if(!sys1.readPdb(orig)) {
		cerr << "Unable to read " << orig << endl;
		return;
	}
	System sys2;
	if(!sys2.readPdb(repacked)) {
		cerr << "Unable to read " << repacked << endl;
		return;
	}

	SasaCalculator sc(sys1.getAtomPointers());
	sc.calcSasa();

	for(int i = 0; i < sys1.positionSize(); i++) {
		Residue & o = sys1.getResidue(i);
		if(sys2.identityExists(o.getIdentityId())) {
			Residue & r = sys2.getLastFoundIdentity();
			string resName = r.getResidueName();
			double sasa = o.getSasa();
			string group = "UNKNOWN";
			if(refSasa.find(resName) != refSasa.end()) {
				if(sasa < refSasa[resName]) {
					group = "BURIED";
				} else {
					group = "EXPOSED";
				}
			}
			cRS.analyze(o,r,group);
		} else {
			cerr << "Missing " << o.getIdentityId() << " in " << repacked << endl;
		}
	}
}

void readListFile(string filename, vector<string>& pdbs) {
	ifstream f;
	f.open(filename.c_str());
	if(f.is_open()) {
		string line;
		while(getline(f,line)) {
			line = MslTools::trim(MslTools::uncomment(line));
			if(line.length() > 0) {
				pdbs.push_back(line);
			}
		}
	}
}

int main(int argc, char* argv[]) {
	Options opt = parseOptions(argc,argv);
	if(opt.errorFlag) {
		cerr << opt.errorMessages;
		exit(0);
	}

	if(opt.warningFlag) {
		cerr << opt.warningMessages;
	}

	// Read list of original files
	vector<string> orig;
	readListFile(opt.origList,orig);

	// Read list of repacked files
	vector<string> repacked;
	readListFile(opt.repackedList,repacked);

	int numPdbs = orig.size();

	if(orig.size() > repacked.size()  ) {
		numPdbs = repacked.size();
	}

	// compute and store counts
	ChiRecoveryStatistics cRS(opt.chiDefFile,opt.tolerance);


	cout << "#Computing chi recovery over " << numPdbs << " PDBs." << endl;
	for(int i = 0; i < numPdbs; i++ ) {
		computeRecovery(cRS,orig[i],repacked[i],opt.refSasa);
	}

	// print stats
	cRS.printStats();
}


void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help() {
	cout << "Program: " << programName << endl;
	cout << "Program description: " << programDescription << endl;
	cout << "Program author: " << programAuthor << endl;
	cout << "Program version: " << programVersion << " " << programDate << endl;
	cout << "MSL version: " << mslVersion << " " << mslDate << endl;
	cout << endl;
	usage();
}

void usage() {
	cout << "Usage:" << endl;
	cout << endl;
	cout << "getChiRecovery --origpdbs <filename> --repackedpdbs <filename> --chideffile <filename> [--tolerance <20.0> --sasathreshfile <filename]\n";
	cout << endl;
	cout << "The sasa threshold file can be used to specify the sasa threshold that defines buried/exposed residues" << endl;
}

Options parseOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	opt.required.push_back("origpdbs");
	opt.required.push_back("repackedpdbs");
	opt.required.push_back("chideffile");

	opt.optional.push_back("sasathreshfile");
	opt.optional.push_back("tolerance");
	
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.setShortOptionEquivalent(opt.equivalent);

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
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
	
	if (OP.countOptions() == 0){
		usage();
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.errorMessages = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;
	opt.errorMessages = "";
	opt.warningMessages = "";

	

	opt.origList = OP.getString("origpdbs");
	if (OP.fail()){
		cerr << "ERROR 1111 origpdbs not specified.\n";
		opt.errorFlag = true;
		opt.errorMessages += "ERROR 1111 origpdbs not specified\n";
	}
	opt.repackedList = OP.getString("repackedpdbs");
	if (OP.fail()){
		opt.errorFlag = true;
		opt.errorMessages += "ERROR 1112 repackedpdbs not specified\n";
	}

	opt.chiDefFile = OP.getString("chideffile");
	if (OP.fail()){
		opt.errorFlag = true;
		opt.errorMessages += "ERROR 1113 chideffile not specified\n";
	}

	/*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti

	Using < 25% of the sasa as the cutoff for buried residues
      */

	opt.refSasa["ALA"] = 0.25 * 116.40;
	opt.refSasa["ARG"] = 0.25 * 249.26;
	opt.refSasa["ASN"] = 0.25 * 168.87;
	opt.refSasa["ASP"] = 0.25 * 155.37;
	opt.refSasa["CYS"] = 0.25 * 141.48;
	opt.refSasa["GLU"] = 0.25 * 187.16;
	opt.refSasa["GLN"] = 0.25 * 189.17;
	opt.refSasa["GLY"] = 0.25 * 83.91;
	opt.refSasa["HSD"] = 0.25 * 198.51;
	opt.refSasa["HSP"] = 0.25 * 198.51;
	opt.refSasa["HSE"] = 0.25 * 198.51;
	opt.refSasa["HIS"] = 0.25 * 198.51;
	opt.refSasa["ILE"] = 0.25 * 189.95;
	opt.refSasa["LEU"] = 0.25 * 197.99;
	opt.refSasa["LYS"] = 0.25 * 207.49;
	opt.refSasa["MET"] = 0.25 * 210.55;
	opt.refSasa["PHE"] = 0.25 * 223.29;
	opt.refSasa["PRO"] = 0.25 * 144.80;
	opt.refSasa["SER"] = 0.25 * 125.68;
	opt.refSasa["THR"] = 0.25 * 148.06;
	opt.refSasa["TRP"] = 0.25 * 265.42;
	opt.refSasa["TYR"] = 0.25 * 238.30;
	opt.refSasa["VAL"] = 0.25 * 162.24;

	
	opt.sasaThreshFile = OP.getString("sasaThreshFile");
	if (OP.fail()){
		opt.warningFlag = true;
		opt.warningMessages += "WARNING no sasaThreshFile set, using internal thresholds.\n";
	} else {
		ifstream f;
		f.open(opt.sasaThreshFile.c_str());
		if(f.is_open()) {
			string line;
			while(getline(f,line)) {
				MslTools::uncomment(line);
				vector<string> toks = MslTools::tokenizeAndTrim(line);
				if(toks.size() >= 2) {
					opt.refSasa[toks[0]] = MslTools::toDouble(toks[1]);
				}
			}
		} else {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read file " + opt.sasaThreshFile + "\n";
		}
	}

	opt.tolerance = OP.getDouble("tolerance");
	if (OP.fail()){
		opt.warningFlag = true;
		opt.warningMessages += "WARNING no tolerance set, using 20.\n";
		opt.tolerance = 20.0;
	} 

	return opt;
}



