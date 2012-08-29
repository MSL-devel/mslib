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

// MSL Includes
#include "System.h"
#include "OptionParser.h"
#include "MslTools.h"
#include "CharmmTopologyResidue.h"
#include "CharmmTopologyReader.h"
#include "CharmmSystemBuilder.h"
#include "SystemRotamerLoader.h"
#include "PolymerSequence.h"
#include "Position.h"
#include "Residue.h"
#include "MslTools.h"
#include "findClashes.h"

// STL Includes
#include <iostream>
#include <string>
#include <signal.h>

using namespace MSL;
using namespace std;



int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	cout << "Read in PDB structure"<<endl;
	System sys;
	sys.readPdb(opt.pdb);


	for (uint i = 0; i < sys.positionSize();i++){
	  Residue &r1 = sys.getResidue(i);

	  for (uint j = i+2; j < sys.positionSize();j++){
	    Residue &r2 = sys.getResidue(j);

	    if (r1.getResidueName() == "HOH" || r2.getResidueName() == "HOH") continue;
	    for (uint a1 = 0; a1 < r1.size();a1++){
	      Atom &at1 = r1.getAtom(a1);
	      for (uint a2 = 0; a2 < r2.size();a2++){

		Atom &at2 = r2.getAtom(a2);

		if (at1.getElement() == "H" || at2.getElement() == "H") continue;

		if (opt.interfaceOnly && at1.getChainId() == at2.getChainId()) continue;

		if (opt.atomPair != ""){
		  if (   (MslTools::stringf("%s-%s",at1.getName().c_str(),at2.getName().c_str()) != opt.atomPair) &&
			 (MslTools::stringf("%s-%s",at2.getName().c_str(),at1.getName().c_str()) != opt.atomPair) ){
		    continue;
		  }

		}

		if (opt.elementPair != ""){
		  if (   (MslTools::stringf("%s-%s",at1.getElement().c_str(),at2.getElement().c_str()) != opt.elementPair) &&
			 (MslTools::stringf("%s-%s",at2.getElement().c_str(),at1.getElement().c_str()) != opt.elementPair) ){
		    continue;
		  }

		}

		double dist = at1.distance(at2);
		if (dist > opt.tol){
		  continue;
		}

		cout << "CLASH "<<dist<<" "<<at1<<" and "<<at2<<endl;
		if (opt.printAll){
		  cout << MslTools::stringf("%12s %12s %8.3f\n", MslTools::getAtomId(at1.getChainId(),at1.getResidueNumber(),at1.getResidueIcode(),at1.getName()).c_str(),
					                           MslTools::getAtomId(at2.getChainId(),at2.getResidueNumber(),at2.getResidueIcode(),at2.getName()).c_str(),
					    dist);
		}

	      }
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
    OP.setAllowed(opt.optional);
    //    OP.setShortOptionEquivalent(opt.equivalent);
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"


    if (OP.countOptions() == 0){
        cout << "Usage:" << endl;
        cout << endl;
        cout << "pdb PDB"<<endl;
        exit(0);
    }

    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	    cerr << "ERROR 1111 no pdb file"<<endl;
    }
    
    opt.tol = OP.getDouble("tol");
    if (OP.fail()){
      opt.tol = 1.88;
      cerr << "WARNING tol set to default: "<<opt.tol<<endl;
    }

    opt.interfaceOnly = OP.getBool("interfaceOnly");
    if (OP.fail()){
      opt.interfaceOnly = false;
    }

    opt.atomPair = OP.getString("atomPair");
    if (OP.fail()){
      opt.atomPair = "";
    }
    opt.elementPair = OP.getString("elementPair");
    if (OP.fail()){
      opt.elementPair = "";
    }
    opt.printAll = OP.getBool("printAll");
    if (OP.fail()){
      opt.printAll = false;
    }
    return opt;
}
