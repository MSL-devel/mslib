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

#include <string>
#include "grepSequence.h"
#include "OptionParser.h"
#include "System.h"
#include "Transforms.h"
#include "RegEx.h"
#include "MslTools.h"

using namespace std;

using namespace MSL;


int main(int argc, char *argv[]){

	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Read-in list of PDBS
	vector<string> lines;
	MslTools::readTextFile(lines,opt.pdbs);

	// Regular expression object
	RegEx re;

	// For each PDB in list..
	bool firstPDB = true;
	AtomPointerVector ref;
	for (uint i = 0; i < lines.size();i++){

		string sysFileName = MslTools::getFileName(lines[i]);
		System sys;
		sys.readPdb(lines[i]);
		AtomPointerVector &sysAts = sys.getAtomPointers();
		sysAts.saveCoor("pre");

		// For each chain
		for (uint c = 0; c < sys.chainSize();c++){
			Chain &ch = sys.getChain(c);

			
			vector<pair<int,int> > matches = re.getResidueRanges(ch,opt.regex);

			for (uint m = 0; m < matches.size();m++){
				
				AtomPointerVector tmp;
				for (uint r = matches[m].first; r < matches[m].second;r++){

					if (firstPDB){
						ref.push_back(new Atom(ch.getResidue(r)("CA")));
					} else {
						tmp.push_back(new Atom(ch.getResidue(r)("CA")));
					}
				}

				if (firstPDB){
					firstPDB = false;
					continue;

				}
				tmp.saveCoor("pre");

				// Now align tmp to ref, apply it to system.
				Transforms t;
				t.rmsdAlignment(tmp,ref,sysAts);

				// Align again for RMSD purposes
				tmp.applySavedCoor("pre");
				t.rmsdAlignment(tmp,ref);

				double rmsd = tmp.rmsd(ref);



				fprintf(stdout, "%25s %1s %3d-%3d %8.3f\n", sysFileName.c_str(), ch.getChainId().c_str(), matches[m].first, matches[m].second,rmsd);
				
				char fname[200];
				sprintf(fname,"%s/%s-%1s-%03d.pdb",opt.outdir.c_str(),sysFileName.c_str(),ch.getChainId().c_str(),matches[m].first);
				sys.writePdb(fname);

				sysAts.applySavedCoor("pre");

				tmp.deletePointers();
				
			}
			
		}
		
	}

}

Options setupOptions(int theArgc, char * theArgv[]){

	// Create the options
	Options opt;
	

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: grepSequence conf" << endl;
		cout << endl;
		cout << "\n";
		cout << "pdblist LIST\n";
		cout << "regex   G...G\n";
		cout << "#outdir .\n";
		cout << endl;
		exit(0);
	}

	opt.configFile = OP.getString("config");
	
	if (opt.configFile != "") {
		OP.readFile(opt.configFile);
		if (OP.fail()) {
			string errorMessages = "Cannot read configuration file " + opt.configFile + "\n";
			cerr << "ERROR 1111 "<<errorMessages<<endl;
		}
	}

	opt.pdbs = OP.getString("pdblist");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.regex = OP.getString("regex");
	if (OP.fail()){
		cerr << "ERROR 1111 no regex specifed."<<endl;
		exit(1111);
	}
	opt.outdir = OP.getString("outdir");
	if (OP.fail()){
		cerr << "WARNING 1111 no outdir specified using current directory"<<endl;
		opt.outdir = ".";
	}

	return opt;
}
