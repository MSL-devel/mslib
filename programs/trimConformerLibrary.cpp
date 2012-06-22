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

#include "RotamerLibrary.h"
#include "OptionParser.h"
#include "MslTools.h"

using namespace std;
using namespace MSL;

string programName = "trimConformerLibrary";
string programDescription = "The conformer library distributed with MSL is large and takes significant time to read into memory. This program reads such a library with defined levels and creates a smaller library at the input \"level\" identifier.";
string programAuthor = "Sabareesh Subramaniam";
string programVersion = "1.0.0";
string programDate = "June 22 2012";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;


// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("conflib");
		required.push_back("level");
		allowed.push_back("outfile");
		allowed.push_back("version"); // --version
		allowed.push_back("help"); // --help
		allowed.push_back("v"); // -v is equivalent to --version
		allowed.push_back("h"); // -h is equivalent to --help
		equivalent.push_back(vector<string>());
		equivalent.back().push_back("v");
		equivalent.back().push_back("version");
		equivalent.push_back(vector<string>());
		equivalent.back().push_back("h");
		equivalent.back().push_back("help");

	}

	// Storage for the vales of each optional
	std::string confLib;
	std::string level;
	std::string outfile;

	bool help;
	bool version;

	// Storage for different types of options
	std::vector<std::string> required;
	std::vector<std::string> allowed;
	std::vector<std::vector<std::string> > equivalent;

};

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --conflib libfile --level level_id [--outfile filename]" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"

	if (OP.countOptions() == 0){
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
		usage();
		exit(0);
	}

	opt.confLib = OP.getString("conflib");
	if (OP.fail()){
		cerr << "ERROR 1111 conflib not specified.\n";
		exit(1111);
	}
	opt.level = OP.getString("level");
	if (OP.fail()){
		cerr << "ERROR 1111 level not specified.\n";
		exit(1111);
	}

	// This is not implemented yet, but is a good idea (dwkulp 3/28/10)
	opt.outfile = OP.getString("outfile");
	if (OP.fail()){
		opt.outfile = opt.level + "_" + MslTools::pathTail(opt.confLib);
		cerr << "WARNING 1111 outfile not specified writing to " << opt.outfile << ".\n";
	}
	
	return opt;
}

int main (int argc, char* argv[]) {
	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Read the library file
	RotamerLibrary rotlib;
	if(!rotlib.readFile(opt.confLib)) {
		cerr << "Unable to read " << opt.confLib << endl;
		exit(0);
	}

	// trim all libraries in the file
	vector<string> libraries = rotlib.getLibraryNames();

	for(vector<string>::iterator it = libraries.begin(); it != libraries.end(); it++) {
		if(!rotlib.trimToLevel(*it,opt.level)) {
			cerr << "Unable to trim library " << *it << " in " << opt.confLib << endl;
			exit(0);
		}
	}

	if(!rotlib.writeFile(opt.outfile)) {
		cerr << "ERROR Unable to write " << opt.outfile << endl;
		exit(0);
	} else {
		cerr << "SUCCESS Trimmed library written to " << opt.outfile << endl;
		exit(0);
	}



}
