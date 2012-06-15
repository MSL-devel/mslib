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
#include "MslTools.h"
#include "RotamerLibrary.h"
#include "OptionParser.h"

#include<fstream>
#include<map>

using namespace MSL;
using namespace std;


string programName = "createEBL";
string programDescription = "This program reads a rotamer library and a raw energy table for one residue type.\n It then sorts the conformers based on their ability to fit natural environements.\n A conformer fits an environment if its energy is less than a specified offset from the best energy for that environment [lower of i) energy of the best rotamer or ii) crystal energy]";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "April 03 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


struct Options {
	string resName; // The residuetype in the rotlibFile and energyTable
	string energyTableFile; // File containing the raw energy table
	string rotlibFile; // The rotlibFile that produced the energyTable

	string outputRotlibFile; // Sorted rotlib will be written to this file

	int numRotsToSort; // optional arguement to say how many rotamers to sort - can ask the program to sort just the first n rotamers
	double offset; // optional argument to specify the offset from best energy to use as threshold

	string configfile; // configFile


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
	cout << " % createEBL \n --resname <resName> \n --energytablefile <energyTableFile> \n --rotlibfile <rotlibfile> \n --outputrotlibfile <outputRotlibFile> \n" << endl ;
	cout << endl;
	cout << "Optional Parameters " << endl;
	cout << "--numrotstosort <numRotsToSort> --configfile <configfile> --offset <double - offset from best energy to use as threshold>" << endl;
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
	opt.required.push_back("resname"); 
	opt.required.push_back("energytablefile"); 
	opt.required.push_back("rotlibfile"); // rotamerLibrary 
	opt.required.push_back("outputrotlibfile"); 

	opt.allowed.push_back("numrotstosort"); 
	opt.allowed.push_back("offset"); 
	opt.allowed.push_back("configfile"); 

	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("help"); // --help
	opt.allowed.push_back("v"); // -v is equivalent to --version
	opt.allowed.push_back("h"); // -h is equivalent to --help

	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
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

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	
	opt.resName= OP.getString("resname");
	if (OP.fail()) {
		opt.errorMessages = "resname not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.energyTableFile= OP.getString("energytablefile");
	if (OP.fail()) {
		opt.errorMessages = "energytablefile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.rotlibFile = OP.getString("rotlibfile");
	if (OP.fail()) {
		opt.errorMessages = "rotlibfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.outputRotlibFile = OP.getString("outputrotlibfile");
	if (OP.fail()) {
		opt.errorMessages = "outputrotlibfile not specified\n";
		opt.errorFlag = true;
		return opt;
	}
	opt.numRotsToSort= OP.getInt("numrotstosort");
	if (OP.fail()) {
		opt.warningMessages = "numrotstosort not specified - sorting all\n";
		opt.warningFlag = true;
		opt.numRotsToSort = -1;
	}

	map<string,double> offset;
	
/*
	If no offsets are specified we can use this
	All bonded terms  +  vdw + hb
	/data01/sabs/EnergyTables/tables/bonded_vdw_hb/cutoffs.txt
*/
	offset["ALA"] =  1.0;  
	offset["ARG"] =  8.0; 
	offset["ASN"] =  1.5;
	offset["ASP"] =  1.5;
	offset["CYS"] =  1.0;
	offset["GLN"] =  4.0;
	offset["GLU"] =  4.5;
	offset["HSD"] =  2.5;
	offset["HSE"] =  2.5;
	offset["HSP"] =  6.5;
	offset["ILE"] =  1.0;
	offset["LEU"] =  1.0;
	offset["LYS"] =  5.0;
	offset["MET"] =  3.0;
	offset["PHE"] =  1.5;
	offset["SER"] =  2.0;
	offset["THR"] =  1.0;
	offset["TRP"] =  3.5;
	offset["TYR"] =  2.0;
	offset["VAL"] =  1.0;

	opt.offset = OP.getDouble("offset");
	if(OP.fail()) {
		if(offset.find(opt.resName) != offset.end()) {
			opt.warningMessages = "offset not specified - using " + MslTools::doubleToString(offset[opt.resName]) + "\n";
			opt.warningFlag = true;
			opt.offset = offset[opt.resName];
		} else {
			opt.errorMessages = "Unknown resName " + opt.resName + "\n";
			opt.errorFlag = true;
			return opt; 
		}
	}

		
	// Print Configuration File / Commmand Line Options
	OP.printConfFile();

	
	// return the Options structure
	return opt;

}




class Table {
	
	public:

	void build(string file, double _offset) {
		ifstream in;
		in.open(file.c_str());
		if(!in.is_open()) {
			cerr << "unable to open " << file << endl;
			exit(0);
		}
		offset = _offset; 
		upperLimit = offset;
		string confline;
		while(getline(in,confline)) {
			vector<string> val = MslTools::tokenizeAndTrim(confline);
			if(val.size() < 3) {
				continue;
			}
			t.push_back(vector<double> ());
			// File format: 1A2J A,23,LEU xtalE rot1E rot2E ....
			// xtalE = crystal Energy 
			// rotNE = energy for Nth rotamer

			// The best energy will be either x-tal energy or the best rotamer energy
			double best = MslTools::toDouble(val[2]);

			// skip x-tal energy
			for(int i = 3; i < val.size(); i++) {
				double value = MslTools::toDouble(val[i]);
				t.back().push_back(value);
				if(best > value) {
					best = value;
				}
			}
			bestEnergies.push_back(best);
			deletedCavity.push_back(false);
			oldEnergy.push_back(MslTools::doubleMax); // will check this only if the cavity is covered. Hopefully, it will contain the old cover energy if the cavity was indeed covered by some rotamer.
		}
		for(int i = 0; i < t[0].size(); i++) {
			deletedRotamer.push_back(false);
		}
		//cout << "Built " << endl;
	}

	inline void deleteCavity (int r) {
		if(r < deletedCavity.size()) {
			deletedCavity[r] = true;
		}
	}
	inline void recoverCavity (int r) {
		if(r < deletedCavity.size()) {
			deletedCavity[r] = false;
		}
	}

	inline void deleteRotamer (int c) {
		if(c < deletedRotamer.size()) {
			deletedRotamer[c] = true;
		}
	}

	int getBest () {
		vector<int> cover(t[0].size(),0);
		
		for(int i = 0; i < t.size(); i++) {
			if(deletedCavity[i] == true) {
				continue;
			}
			for(int j = 0; j < t[i].size(); j++) {
				if(deletedRotamer[j] == true) {
					continue;
				}
				if(t[i][j] < offset + bestEnergies[i]) {
					cover[j]++;
				}
			}
		}
		int best  = 0;
		int best_covered = 0;

		// find the best rotamer
		for(int i = 0; i < t[0].size(); i++) {
			if(deletedRotamer[i] == true) {
				continue;
			}
			if(cover[i] >= best_covered) {
				best = i;
				best_covered = cover[i];	
			}
		}
		 return best;
	}

	void recoverCavities() {
		for(int i = 0; i < bestEnergies.size(); i++) {
			// if cavity is deleted
			if(deletedCavity[i] == true) {
				// if cavity not covered by old rotamer recover it.
				if(oldEnergy[i] >= offset + bestEnergies[i]) {
					deletedCavity[i] = false;
					covered--;
				}
			}
		}

	}

	void updateOldEnergy(int id) {
		// check if this conformer has better energies for any deleted cavity
		for(int i = 0; i < t.size(); i++) {
			if(deletedCavity[i] == true) {
				if(oldEnergy[i] > t[i][id]) {
					oldEnergy[i] = t[i][id];	
				}
			}
		}

	}

	vector<int> getMaxCover(int n) {
		vector<int> res;
		covered = 0; // number of cavities covered at a particular point in time
		int count = 0;

		// until n conformers have been ranked or all conformers in the library have been ranked or we cannot decrease the offset any further
		while(count < n && count < t[0].size() && offset > 0) {
			int id = getBest();
			count++; // no of conformers ranked so far
			//cout << count << " Covered " << covered << endl;
			//cout << "Best: " << id << endl;
			for(int i = 0; i < t.size();i++) {
				//cout << t[i][id] << " " ;
				if(deletedCavity[i] == true) {
					continue;
				}
				if(t[i][id] < offset + bestEnergies[i] ) {
					covered++;
					deleteCavity(i);

					updateOldEnergy(id); // check if this conformer has better energies for any deleted cavity
					//oldEnergy[i] = t[i][id];
				//	i--;
				}
			}
			//cout << endl;
			// if all cavities have been covered or no new cavities were covered
			offset = offset - upperLimit/n;
			//cout << "Covered " << covered << endl;
			recoverCavities();
			deleteRotamer(id);
			res.push_back(id);
		}
		//cout << "At exit: # rotamers " << count << " Requested " << n << " size " << t[0].size()  << " Threshold " << offset << endl;

		return res;
	}

	int getNumberCovered () {
		return covered;
	}

	void printUncovered() {
		for(int i = 0; i < t.size(); i++) {
			if(deletedCavity[i] == true) {
				continue;
			}
			cout << "(" << i << ") " << endl;
		}
	}
	

	void print() {
		for(int i = 0; i < t.size(); i++) {
			for(int j = 0; j < t[i].size(); j++) {
				cout << t[i][j] << " " ;
			}
			cout << endl;
		}
	}
	int getNumberOfRotamers() {
		return t[0].size();
	}

	private:
	int covered;
	vector<bool> deletedCavity;
	vector<bool> deletedRotamer;
	vector<vector <double> > t; // contains the energy Table
	vector<double> oldEnergy; // contains the old cover energy for each cavity
	vector<double> bestEnergies; // contains the bestEnergies for each cavity
	double offset;
	double upperLimit;

};



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


	Table t;
	t.build(opt.energyTableFile,opt.offset);
	if(opt.numRotsToSort == -1) {
		opt.numRotsToSort = t.getNumberOfRotamers();
	}
	

	vector<int> idx = t.getMaxCover(opt.numRotsToSort);

	//cout << "Number Covered " << idx.size() << " " << t.getNumberCovered() << endl;

	//t.printUncovered();

	//cout << "Number Covered " << idx.size() << " " << t.getNumberCovered() << endl;
	RotamerLibrary rotlib;
	if(!rotlib.readFile(opt.rotlibFile)) {
		cerr << "Unable to read " << opt.rotlibFile << endl;
		exit(0);
	}

	string library = rotlib.getDefaultLibrary();
	
	RotamerLibrary rotlib_out (rotlib);
	rotlib_out.removeAllConformations();

	for(int i = 0; i < idx.size(); i++) {
		//cout << "Rotamer # : " << idx[i] << endl;
		string line = rotlib.getInternalCoorLine(library,opt.resName,idx[i]);
		vector<string> val = MslTools::tokenizeAndTrim(line);
		vector<double> v;

		for(int j = 0; j < val.size(); j++) {
			v.push_back(MslTools::toDouble(val[j]));
		}

		if(!rotlib_out.addConformation(library,opt.resName,v)) {
			cout << "Unable to add " << opt.resName << " to library " << library << endl;
		}
	}

	//cout << "Number of rotamers " << idx.size() << endl;
	if(!rotlib_out.writeFile(opt.outputRotlibFile)) {
		cerr << "Unable to write to " << opt.outputRotlibFile << endl;
	}
}
