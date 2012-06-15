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


#include "OptionParser.h"
#include "System.h"
#include "MslTools.h"
#include "BackRub.h"
#include "MslOut.h"
#include "backrubPdb.h"

using namespace std;

using namespace MSL;

static MslOut MSLOUT("backrubPdb");

int main(int argc, char *argv[]){


	// Option Parser
	Options opt = setupOptions(argc,argv);


	MSLOUT.stream() << "Read the PDB: "<<opt.pdb<<endl;
	// Read-in list of PDBS

	
	System sys;
	sys.readPdb(opt.pdb);

	System orig;
	orig.readPdb(opt.pdb);

	int index = 1;
	do { 

	  // Re-sample 3 times
	  for (uint re = 0; re < 3; re++){
	  
	    MSLOUT.stream()<<"Resample "<<re<<endl;

	    // For each conformation do it..
	    AtomPointerVector &activeAtoms = sys.getAtomPointers();
	    int minAltConf = activeAtoms.getMinAltConf();

	    System newSystem;
	    for (uint i = 0; i < minAltConf;i++){

	      for (uint a = 0; a< activeAtoms.size();a++){
		activeAtoms[a]->setActiveConformation(i);
	      }


	      // A BackRub object
	      BackRub br;

	      string chainId    = sys.getPosition(opt.startRes).getChainId();
	      int startResIndex = sys.getChain(chainId).getPositionIndex(&sys.getPosition(opt.startRes));
	      int endResIndex   = sys.getChain(chainId).getPositionIndex(&sys.getPosition(opt.endRes));

	      //MSLOUT.stream() << "Backrub like this: "<<chainId<<" "<<opt.startRes<<"("<<startResIndex<<") "<<opt.endRes<<"("<<endResIndex<<")"<<endl;

	    
	      // Do local sampling inside BackRub object
	      br.localSample(sys.getChain(chainId),startResIndex,endResIndex,opt.numSamples);


	      //MSLOUT.stream()<<"Number of conformationsB: "<<br.getAtomPointers().getMinAltConf()<<" "<<br.getAtomPointers().getMaxAltConf()<<endl;

	      newSystem.addAtoms(br.getAtomPointers());

	      //MSLOUT.stream()<<"Number of conformationsT: "<<newSystem.getAtomPointers().getMinAltConf()<<" "<<newSystem.getAtomPointers().getMaxAltConf()<<endl;
	    
	    }


	    //MSLOUT.stream()<<"Number of conformations1: "<<sys.getAtomPointers().getMinAltConf()<<" "<<sys.getAtomPointers().getMaxAltConf()<<endl;
	    sys.reset();
	    sys.addAtoms(newSystem.getAtomPointers());
	    newSystem.reset();

	    MSLOUT.stream()<<"Number of conformations: "<<sys.getAtomPointers().getMinAltConf()<<" "<<sys.getAtomPointers().getMaxAltConf()<<endl;


	  }
	
	  // Write out mulitple pdb files...
	  string brChain = orig.getPosition(opt.startRes).getChainId();
	  AtomPointerVector origChain = orig.getChain(brChain).getAtomPointers();
	  AtomPointerVector origOtherChains;
	  for (uint c = 0; c < orig.chainSize();c++){
	    if (orig.getChain(c).getChainId() != brChain){
	      origOtherChains += orig.getChain(c).getAtomPointers();
	    }
	  }

	  int index = 1;
	  int minAltConf = sys.getAtomPointers().getMinAltConf();
	  for (uint i = 0; i < minAltConf;i++){



	    // Set this alt conformation
	    for (uint a  = 0; a < sys.getAtomPointers().size();a++){
	      sys.getAtom(a).setActiveConformation(i);
	    }

	    if (sys.getAtomPointers().rmsd(origChain) < opt.rmsd){


	      char tmp[100];
	      sprintf(tmp,"%s_%06d.pdb",opt.outPdb.c_str(),index);
	      PDBWriter pout;
	      pout.open((string)tmp);
	      pout.write(sys.getAtomPointers());
	      pout.write(origOtherChains);
	      pout.close();

	      index++;

	      if (index > opt.numModels){
		MSLOUT.stream() << "Done."<<endl;
		exit(9);
	      }
	    
	    }
	  
	  }
	} while (index < opt.numModels);
	  
	//sys.writeMultiplePdbs(opt.outPdb,0.5);
			

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
		cout << "Usage: backrub" << endl;
		cout << endl;
		cout << "pdb PDB\n";
		cout << "startRes Chain,Position\n";
		cout << "endRes Chain,Position\n";
		cout << "numSamples X\n";
		cout << "outPdb PDB\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	opt.startRes = OP.getString("startRes");
	if (OP.fail()){
		cerr << "ERROR 1111 no startRes specified."<<endl;
		exit(1111);
	}
	opt.endRes = OP.getString("endRes");
	if (OP.fail()){
		cerr << "ERROR 1111 no endRes specified."<<endl;
		exit(1111);
	}

	opt.numModels = OP.getInt("numModels");
	if (OP.fail()){
	  opt.numModels = 200;
	  cerr << "WARNING 1111 numModels is defaulting to"<<opt.numModels<<"\n";

	}

	opt.numSamples = OP.getInt("numSamples");
	if (OP.fail()){
	  opt.numSamples = 10;
	  cerr << "WARNING 1111 numSamples is defaulting to"<<opt.numSamples<<"\n";

	}
	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  opt.rmsd = 0.5;
	  cerr << "WARNING rmsd set to "<<opt.rmsd<<endl;
	}

	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  opt.outPdb = MslTools::getFileName(opt.pdb)+".br.pdb";
	  cerr << "WARNING 1111 outPdb is defaulting to "<<opt.outPdb<<endl;

	}

	return opt;
}
