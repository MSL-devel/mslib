//STL Includes
#include<fstream>
#include<string>
#include<vector>
#include<iostream>

//MSL Includes

#include "PDBReader.h"
#include "PDBWriter.h"
#include "System.h"
#include "Matrix.h"
#include "MatrixWindow.h"
#include "DistanceMatrix.h"
#include "OptionParser.h"

#include "Transforms.h"
#include "MslTools.h"
#include "AtomSelection.h"
#include "ManageDistanceMatrixResults.h"
#include "multiSearchDM.h"

using namespace std;
using namespace MslTools;


int main(int argc, char *argv[]){
    
    // Option Parser
    Options opt = setupOptions(argc,argv);

    ifstream fs2;

    //create system and dm for first PDB
    PDBReader reader;
    reader.open(opt.inputPDB);
    reader.read();
    reader.close();

    System *constSys = new System(reader.getAtoms());
    DistanceMatrix constDM;
    
    //add CA atoms to the atom vectors
    for (int j=0; j<constSys->residueSize(); j++){
	Residue &tempRes=constSys->getResidue(j);
	if (tempRes.exists("CA")){
	    constDM.addAtom(tempRes("CA"));
	}
    }//end for on j

    
    //fill the DistanceMatrix and set window size
    constDM.setGeneralWinSize(opt.windowSize);
    constDM.createDistanceMatrix();
    constDM.setIntraChain(opt.intraChainCompare);
    constDM.setPDBid(opt.inputPDB);
    constDM.setDebug(opt.debug);

    //create matrix windows
    constDM.createMatrixWindows();

    delete(constSys);


    if (constDM.getMatrixWindows().size()==0){
	    cout<<"Uh-oh.All the windows got filtered in the PDB you wanted to compare against."<<endl;
	    exit(111);
    }

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
    vector<DistanceMatrix> DMVec(list.size());

    // A system object for each PDB
    vector<System*> sysVec(list.size(), NULL);

    // Create DistanceMatrix and System Objects for list of PDBs
    for(int i=0; i<list.size(); i++){

	cout<<i<<"create sys and dm."<<endl;

	PDBReader rAv;
	rAv.open(list[i]);
	rAv.read();
	rAv.close();

	sysVec[i] =new System(rAv.getAtoms());
    
	//add CA atoms to the atom vectors
	for (int j=0; j<sysVec[i]->residueSize(); j++){
	    Residue &tempRes=sysVec[i]->getResidue(j);
	    if (tempRes.exists("CA")){
		//only add CA if it is on a helix
		string segID = tempRes("CA").getSegID();
		//if(segID == "" || segID.at(0) == 'H'){
		    DMVec[i].addAtom(tempRes("CA"));
		    //	}
	    }
	}//end for on j
	//fill the DistanceMatrix and set window size
	DMVec[i].setGeneralWinSize(opt.windowSize);
	DMVec[i].createDistanceMatrix();
	DMVec[i].setIntraChain(opt.intraChainCompare);
	DMVec[i].setPDBid(list[i]);

	//create matrix windows
	DMVec[i].createMatrixWindows();
	
	delete(sysVec[i]);

    }//end for on i


    //ManageResults to take care of printing/sorting at end
    ManageDistanceMatrixResults resultManager;

	    
    for(int i=0; i<DMVec.size(); i++){

	    cout<< "Trying "<<DMVec[i].getPDBid()<<" ("<<i<<") # Residues: "<<DMVec.size()<<" Number of MatrixWindows to compare: "<<DMVec[i].getMatrixWindows().size();

	    //don't compare if all of the windows got filtered out
	    if (DMVec[i].getMatrixWindows().size() == 0){
		    cout << " Sorry Zero Matrix Windows !"<<endl;
		    continue;
	    }
	    
	    cout <<endl;

	    vector<DistanceMatrixResult> resultsToAdd;

	    if(opt.searchCriteria=="standard"){
		resultsToAdd = constDM.multiCompareAllWindows(DMVec[i], DistanceMatrix::standard, opt.numberOfIterations);
	    }//end if
	    if(opt.searchCriteria=="diagonal"){
		resultsToAdd = constDM.multiCompareAllWindows(DMVec[i], DistanceMatrix::diag, opt.numberOfIterations);
       	    }
	    if(opt.searchCriteria=="doubleDiagonal"){
		resultsToAdd = constDM.multiCompareAllWindows(DMVec[i], DistanceMatrix::doubleDiag, opt.numberOfIterations);
	    }
	    if(opt.searchCriteria=="minDistance"){
		resultsToAdd = constDM.multiCompareAllWindows(DMVec[i], DistanceMatrix::minDist, opt.numberOfIterations);
	    }
	    if(opt.searchCriteria=="minDistanceRow"){
		resultsToAdd = constDM.multiCompareAllWindows(DMVec[i], DistanceMatrix::minDistRow, opt.numberOfIterations);
	    }

	    bool addFlag = false;
 	    for (uint j = 0; j< resultsToAdd.size();j++){
 		    if (opt.likenessTolerance == MslTools::doubleMax || resultsToAdd[j].getLikeness() <= opt.likenessTolerance){
			    addFlag = true;
			    break;
 		    }
	    }

	    if (addFlag &&  resultsToAdd.size() > 0){

		    resultManager.addResults(resultsToAdd);	    		    
	    }


    }//end for on i
    
    cout << "Printing"<<endl;
    resultManager.setAlignPdbs(opt.alignPdbs);
    resultManager.setRmsdTol(opt.rmsdTol);
    resultManager.printResults();

    cout << "Done."<<endl;
    return 0;
}




void getRMSD(MatrixWindow *_win1, MatrixWindow *_win2, System *_sys2){
    Transforms t;
    AtomVector ca1 = (*_win1).getSmallAVec();
    AtomVector ca2 = (*_win2).getSmallAVec();
    AtomVector a2 = (*_sys2).getAtoms();

    bool result = t.align(ca2, ca1, a2);
    if(!result){
	cout<<"Alignment has failed!"<<endl;
	exit(1211);
    }

    double r=ca1.rmsd(ca2);
    fprintf(stdout, "RMSD: %8.3f\n", r);

    //write out aligned pdb

    /*char a[80];
      sprintf(a, "/snap/cluster/jdegrado/pizza/%03d.aligned.%03d.pdb",j,i);

      PDBWriter w(a);
      w.write(a2);
      w.close();*/

}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;

	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "multiSearchDM --inputPDB PDB --pdbList LIST --searchCriteria SEARCH_STRING --windowSize NUM [ --numberOfIterations NUM --allowIntraChainCompare  --likenessTolernace TOL --alignPdbs --debug ]\n";
		cout << endl<<"\tsearchCriteria can be:\n";
		cout << "\tstandard       - sum of the diff. for every item between two MatrixWindows.\n";
		cout << "\tdiagnol        - sum of the diff. for the diagnol (top left-bottom right) between two MatrixWindows.\n";
		cout << "\tdoubleDiagnol  - sum of the diff. for both diagnols between two MatrixWindows.\n";
		cout << "\tminDistance    - sum of the diff. between the minimal distance for each row,col of the given MatrixWindows.\n";
		cout << "\tminDistanceRow - sum of the diff. between the minimal distance for each row of the given MatrixWindows.\n";
		exit(0);
	}

	opt.inputPDB = OP.getString("inputPDB");
	if (OP.fail()){
		cerr << "ERROR 1111 inputPDB not specified.\n";
		exit(1111);
	}

	opt.pdbList = OP.getString("pdbList");
	if (OP.fail()){
		cerr << "ERROR 1111 pdbList not specified.\n";
		exit(1111);
	}

	opt.searchCriteria = OP.getString("searchCriteria");
	if (OP.fail()){
		cerr <<"ERROR 1111 searchCriteria not specified."<<endl;
		exit(1111);
	}

	opt.windowSize = OP.getInt("windowSize");
	if (OP.fail()){
		cerr << "ERROR 1111 windowSize not specified."<<endl;
		exit(1111);
	}

	opt.numberOfIterations = OP.getInt("numberOfIterations");
	if (OP.fail()){
		opt.numberOfIterations = 1;
	}

	opt.intraChainCompare = OP.getBool("allowIntraChainCompare");
	if (OP.fail()){
		opt.intraChainCompare = false;
	}


	opt.likenessTolerance = OP.getDouble("likenessTolerance");
	if (OP.fail()){
		opt.likenessTolerance = MslTools::doubleMax;
	}

	opt.alignPdbs         = OP.getBool("alignPdbs");
	if (OP.fail()){
		opt.alignPdbs = false;
	}

	opt.rmsdTol   = OP.getDouble("rmsdTol");
	if (OP.fail()) {
		opt.rmsdTol = 2.0;
	}


	opt.debug = OP.getBool("debug");
	if (OP.fail()){
		opt.debug = false;
	}


	return opt;
}
