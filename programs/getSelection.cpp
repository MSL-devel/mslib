/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan
 Sabareesh Subramaniam, Ben Mueller

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
#include <vector>
#include <map>
#include "getSelection.h"
#include "OptionParser.h"
#include "System.h"
#include "ResidueSelection.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include "CharmmTopologyReader.h"
#include "AtomContainer.h"

using namespace std;

using namespace MSL;


int main(int argc, char *argv[]){


	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Read-in list of PDBS
	System sys;
	sys.readPdb(opt.pdb);



	if (opt.resSel != ""){
		ResidueSelection sel(sys);
		vector<Residue *> &res = sel.select(opt.resSel);
		
		if (opt.sequence){
			cout << "SEQ: ";
			for (uint i = 0; i < res.size();i++){
				cout << MslTools::getOneLetterCode(res[i]->getResidueName());
			}
			cout <<endl;
		} else {
		        if (opt.length){
			  cout << res.size()<<endl;
			} else {
			  for (uint i = 0; i < res.size();i++){
			    cout <<res[i]->toString()<<endl;
			  }
			}
		}

	}


	if (opt.atomSel != ""){
		AtomSelection sel(sys.getAtomPointers());
		AtomPointerVector a = sel.select(opt.atomSel);

		if (opt.sequence){

			cout << "SEQ: ";
			for (uint i = 0; i < a.size();i++){
				cout << MslTools::getOneLetterCode(a[i]->getResidueName());
			}
			cout <<endl;
		} else {


			if (opt.outPdb != ""){
			  cout << "Writing "<<opt.outPdb<<endl;
				PDBWriter pout;
				pout.open(opt.outPdb);
				pout.write(a);
				pout.close();
			} else {

				for (uint i = 0; i < a.size();i++){
					cout <<a[i]->toString()<<endl;
				}
			}
			
		}
	}


	if (opt.charmmTop != ""){

	        CharmmTopologyReader ctr(opt.charmmTop);
		ctr.read();		
		vector<CharmmTopologyResidue*> residues = ctr.getResidues();
		map<string,bool> charmmResidueNames;
		for (uint i = 0; i < residues.size();i++){
		    charmmResidueNames[residues[i]->getName()] = true;
		}

		if (opt.addCharmmHis){
		    charmmResidueNames["HIS"] = true;
		}
		
		AtomContainer storeValidResidues;
		for (uint i = 0; i < sys.positionSize();i++){
		    if (charmmResidueNames[sys.getPosition(i).getResidueName()]){
			    storeValidResidues.addAtoms(sys.getPosition(i).getCurrentIdentity().getAtomPointers());
			    if (opt.outPdb == ""){
				    cout << "Residue: "<<sys.getPosition(i).getCurrentIdentity().toString()<< " is found in Charmm Toplogy File: "<<opt.charmmTop<<endl;
			    }

		    } else {
			    cout << "Residue: "<<sys.getPosition(i).getCurrentIdentity().toString()<< " is NOT found in Charmm Toplogy File: "<<opt.charmmTop<<endl;
		    }
	
		}

		if (opt.outPdb != ""){
		    	PDBWriter pout;
	    		pout.open(opt.outPdb);
			pout.write(storeValidResidues.getAtomPointers());
			pout.close();
		} 


			
	}
}

Options setupOptions(int theArgc, char * theArgv[]){

	// Create the options
	Options opt;
	

	// Parse the options
	OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.readArgv(theArgc, theArgv);
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option

	if (OP.countOptions() == 0){
		cout << "Usage: grepSequence conf" << endl;
		cout << endl;
		cout << "pdb PDB\n";
		cout << "\n#One of these types of selections\n";
		cout << "resSel    SELE_STATEMENT\n";
		cout << "atomSel   SELE_STATEMENT\n";
		cout << "charmmTop CHARMM_TOPOLOGY_FILE\n";
		cout << "\n#Optionally ask for only a sequence string back\n";
		cout << "sequence\n";
		cout << endl;
		cout << "Examples:\n\n";
		cout << "getSelection --pdb foo.pdb --atomSel \"resn ALA+ARG+ASN+ASP+CYS+GLN+GLU+GLY+HIS+HSD+HSE+HSP+ILE+LEU+LYS+MET+PRO+PHE+SER+THR+TRP+TYR+VAL\""<<endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.sequence = OP.getBool("sequence");
	if (OP.fail()){
		opt.sequence = false;
	}
	opt.length   = OP.getBool("length");
	if (OP.fail()){
	  opt.length = false;
	}
	opt.resSel = OP.getString("resSel");
	if (OP.fail()){
		opt.resSel = "";
	}
	opt.atomSel = OP.getString("atomSel");
	if (OP.fail()){
		opt.atomSel = "";
	}

	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  cout << "No output pdb !"<<endl;
		opt.outPdb ="";
	}

	opt.charmmTop = OP.getString("existsInCharmm");
	if (OP.fail()){
	  opt.charmmTop = "";
	}

	if (opt.resSel == "" && opt.atomSel == "" && opt.charmmTop == ""){
		cerr << "ERROR 1111 either resSel or atomSel or charmmTop has to be specified.\n";
		exit(1111);
	}

	if (  (opt.resSel != "" && opt.atomSel != "" && opt.charmmTop != "") ||
	      (opt.resSel != "" && opt.atomSel != "")  ||
	      (opt.resSel != "" && opt.charmmTop != "")  ||
	      (opt.atomSel != "" && opt.charmmTop != "")) {
		cerr << "ERROR 1111 either resSel OR atomSel OR charmmTop has to be specified, but all three or two/three.\n";
		exit(1111);
	}


	opt.addCharmmHis = OP.getBool("addCharmmHis");
	return opt;
}
