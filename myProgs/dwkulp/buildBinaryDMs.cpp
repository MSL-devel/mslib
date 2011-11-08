//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes

#include "System.h"
#include "MatrixWindow.h"
#include "DistanceMatrix.h"
#include "DistanceMatrixDatabase.h"
#include "OptionParser.h"

#include "MslTools.h"
#include "buildBinaryDMs.h"

using namespace std;
using namespace MSL;
using namespace MslTools;


int main(int argc, char *argv[]){
    
    // Option Parser
    Options opt = setupOptions(argc,argv);

    //read in list of PDBs to compare to first PDB
    vector<string> list;
    ifstream fs;

    fs.open(opt.pdbList.c_str());
    if (fs.fail()){
	cerr<<"Cannot open file "<<opt.pdbList<<endl;
	exit(1);
    }

    while(true){
	string line;
	getline(fs, line);

	if(fs.fail()){
	    //no more lines to read, quite the while.
	    break;
	}

	if(line==""){
	    continue;
	}
	list.push_back(line);
    }

    fs.close();

    
    // List of distance matrices, one for each PDB
    DistanceMatrixDatabase dmd;
    dmd.setName("distance matrix database");

    // Create DistanceMatrix and System Objects for list of PDBs
    for(int i=0; i<list.size(); i++){

      cout<<i<<" create sys and dm. "<<list[i]<<endl;

	//PDBReader r1;
	// r1.setSingleAltLocationFlag(true);
	//r1.open(list[i]);
	//r1.read();
	//r1.close();

        System sys;
        sys.readPdb(list[i]);

	//sysVec[i] =new System();
	//sysVec[i]->readPdb(list[i]);
	    
	DistanceMatrix *dm = new DistanceMatrix();

	map<int, bool> nearALoopMap;

	//add CA atoms to the atom vectors
	for (int j=0; j<sys.positionSize(); j++){
	    Residue &tempRes=sys.getResidue(j);
	    if (tempRes.atomExists("CA")){

		//only add CA if it is on a helix
		string segID = tempRes("CA").getSegID();

		bool isALoop = false;
		bool isNearALoop = false;
		bool isChainTerminal = false;
		
		if (j > 1 && j < sys.positionSize()-1 && (sys.getResidue(j).getChainId() != sys.getResidue(j-1).getChainId() || sys.getResidue(j).getChainId() != sys.getResidue(j+1).getChainId())){ isChainTerminal = true; }


		if (segID == "LLLL" || segID == "TTTT" || segID == "SSSS"){
		    isALoop = true;
		    nearALoopMap[j] = true;
		  }

		  // Look ahead 2 residues for loop residues
		  for (uint a2 = j+1; a2 < j+2; a2++){
		    if (a2 < sys.positionSize()){
		      if (sys.getResidue(j).getChainId() == sys.getResidue(a2).getChainId() && sys.getResidue(a2).atomExists("CA")){
			if (sys.getResidue(a2)("CA").getSegID()  == "LLLL" || sys.getResidue(a2)("CA").getSegID() == "TTTT" || sys.getResidue(a2)("CA").getSegID() == "SSSS"){
			  nearALoopMap[a2] = true;
			  isNearALoop = true;
			}
		      } else {
			break;
		      }
		    }
		  }

		  // If near an N-term loop
		  if (!isChainTerminal && (j < 2 || nearALoopMap[j-2] ||nearALoopMap[j-1])){
		    isNearALoop = true;
		  }


		  if (isALoop || isNearALoop){ 
		    dm->addAtom(tempRes("CA"));
		  }
	    }
	  }//end for on j


	cout << MslTools::getFileName(list[i])<< " "<<dm->getAtomVector().size()<<" residues were added out of "<<sys.positionSize()<<" total residues in the file.\n";

	//fill the DistanceMatrix and set window size
	dm->setGeneralWinSize(opt.windowSize);
	dm->createDistanceMatrix();
	dm->setIntraChain(opt.intraChainCompare);
	dm->setPDBid(list[i]);
	dm->setDiagnolMatrixWindowsOnlyFlag(opt.diagnolMatrixWindowsOnly);


	//create matrix windows
	dm->createMatrixWindows();
	cout << "Number of matrix windows: "<<dm->getDiagnolMWs().size()<<endl;
	dmd.addDistanceMatrix(dm);
	dm = NULL;

	cout << "DONE LOOP : "<<i<<endl;

    }//end for on i


    //dmd.save_checkpoint("dmd.helixPairs.bin");
    dmd.save_checkpoint(opt.outName);
    
}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "buildBinaryDMs --pdbList LIST --windowSize SIZE --allowIntraChainCompare [ --outName OUTNAME ]"<<endl;
		exit(0);
	}
	opt.pdbList = OP.getString("pdbList");
	if (OP.fail()){
		cerr << "ERROR 1111 pdbList not specified.\n";
		exit(1111);
	}


	opt.windowSize = OP.getInt("windowSize");
	if (OP.fail()){
		cerr << "ERROR 1111 windowSize not specified."<<endl;
		exit(1111);
	}

	opt.intraChainCompare = OP.getBool("allowIntraChainCompare");
	if (OP.fail()){
		opt.intraChainCompare = false;
	}

	opt.diagnolMatrixWindowsOnly = OP.getBool("diagnolMWonly");
	if (OP.fail()){
	  opt.diagnolMatrixWindowsOnly = false;
	}
	opt.outName = OP.getString("outName");
	if (OP.fail()){
		opt.outName = "dmd.bin";
	}
	opt.debug = OP.getBool("debug");
	if (OP.fail()){
		opt.debug = false;
	}


	return opt;
}
