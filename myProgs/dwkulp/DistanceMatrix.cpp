#include "DistanceMatrix.h"
#include "MatrixWindow.h"
#include "DistanceMatrixResult.h"
#include <algorithm>
using namespace std;
using namespace MSL;

//constructor of DistanceMatrix
DistanceMatrix::DistanceMatrix(){
    generalWinSize = 10;
    intraChain = false;
    debug = false;
    archiveType = "binary";
    diagnolMatrixWindowsOnly = false;
}

DistanceMatrix::DistanceMatrix(const DistanceMatrix &_dm){
    copy(_dm);
}

//destructor of DistanceMatrix
DistanceMatrix::~DistanceMatrix(){

    for (int i=0; i<listMW.size(); i++){
	delete(listMW[i]);
    }

    mapMW.clear();
}

void DistanceMatrix::copy(const DistanceMatrix &_dm){
    //copy list of matrixwindows.
    for (int i=0; i<_dm.listMW.size(); i++){
	listMW.push_back(new MatrixWindow(*(_dm.listMW[i])));
    }
    
    //copy everything else
    atomVec = _dm.atomVec;
    dimension = _dm.dimension;
    generalWinSize = _dm.generalWinSize;
    intraChain = _dm.intraChain;
    mapMW = _dm.mapMW;
    PDBid = _dm.PDBid;
    weight = _dm.weight;
    debug = _dm.debug;
    archiveType = _dm.archiveType;
    diagnolMatrixWindowsOnly = _dm.diagnolMatrixWindowsOnly;
}

void DistanceMatrix::addAtom(const Atom & _atm){
           
           Atom *a = new Atom(_atm.getAtomId(),_atm.getX(),_atm.getY(),_atm.getZ(),_atm.getElement());
            atomVec.push_back(a);
}

void DistanceMatrix::createDistanceMatrix(){
    dimension = atomVec.size();
    initialize(dimension, dimension);
    for(int i=0; i<dimension; i++){
	for(int j=0; j<dimension; j++){
	    //fills the distance matrix
	    (*this)[i][j] = (*atomVec[i]).distance(*atomVec[j]);
	}
    }
    
}
  
pair<vector<int>, double>  DistanceMatrix::compareAllWindows(DistanceMatrix &_distMat, int choice){
   vector<MatrixWindow*> listOther = _distMat.getMatrixWindows();
    int length = listMW.size();
    int lengthOther = listOther.size();
    vector<int> minCoord(2, 0.0);
    MatrixWindow *minWindow1 = NULL;
    MatrixWindow *minWindow2 = NULL;
    double minLikeness= MslTools::doubleMax;

    for (int i =0; i<length; i++){//loops through listMW to get comparee

	MatrixWindow *win1 = listMW[i];
	for (int j=0; j<lengthOther; j++){//loops through listOther to get comparor

	    MatrixWindow *win2 = listOther[j];
	    
	    double likeness;
	    //decides which compare method from MatrixWindow to call
	    switch(choice){
		    case standard:
			likeness = (*win1).compare((*win2));
			break;
		    case diag:
			likeness = (*win1).compareDiagonal((*win2));
			break;
		    case doubleDiag:
			likeness = (*win1).compareDoubleDiagonal((*win2));
			break;
		    case minDist:
			likeness = (*win1).compareMinDist((*win2));
			break;
		    case minDistRow:
			likeness=(*win1).compareMinRow((*win2));
			break;
		    case minDistCol:
			likeness=(*win1).compareMinCol((*win2));
			break;
		    default:
			cout<<"Invalid argument in function DistanceMatrix::compareAll()."<<endl;
			exit(333);
	    }//end switch
	    
	    if (likeness<minLikeness && abs(likeness - minLikeness) > 0.001){
		    cout << "New likeness: "<<likeness<<endl;
		minLikeness=likeness;
		minWindow1 = win1;
		minWindow2 = win2;
		minCoord[0]= i;
		minCoord[1]=j;

	    }//endif
	    win2=NULL;
	}//end for on j
	win1= NULL;
    }//end for on i

    
    pair<vector<int>, double> coordsAndLikeness(minCoord, minLikeness);

    minWindow1=NULL;
    minWindow2=NULL;

    return coordsAndLikeness;

}

vector<DistanceMatrixResult> DistanceMatrix::multiCompareAllWindows(DistanceMatrix &_distMat, int choice, int _numCompare){
/*
    //which chains to skip:
    map<string, bool> forbiddenIDMat1;
    map<string, bool> forbiddenIDMat2;

    // Maintain the proper spacing between residues of different matrix window pairs. 
    map<string,int> properRegisterRow;
    map<string,int> properRegisterCol;
    map<string,int>::iterator findRegistry;
*/
    //get list of MWs to compare
    vector<MatrixWindow*> listOther = _distMat.getMatrixWindows();
    int length = listMW.size();
    int lengthOther = listOther.size();

    //vector of objects to return
    vector<DistanceMatrixResult> returnVec;

    //loop over number of times to compare all
    for(int k=0; k<_numCompare; k++){

	//minimum windows and indices
	MatrixWindow *minWindow1 = NULL;
	MatrixWindow *minWindow2 = NULL;
	double minLikeness = 1000000;
	int minIndex1=0;
	int minIndex2=0;

	//IDs to skip in the future
	string i1IDWin;
	string j1IDWin;
	string i2IDWin;
	string j2IDWin;


    //which chains to skip:
    map<string, bool> forbiddenIDMat1;
    map<string, bool> forbiddenIDMat2;

    // Maintain the proper spacing between residues of different matrix window pairs. 
    map<string,int> properRegisterRow;
    map<string,int> properRegisterCol;
    map<string,int>::iterator findRegistry;


	for (int i =0; i<length; i++){//loops through listMW to get compare

	    MatrixWindow *win1 = listMW[i];
	    for (int j=0; j<lengthOther; j++){//loops through listOther to get comparor
		
		MatrixWindow *win2 = listOther[j];

		//get the ID (seg or chain) so we can filter ones we want to skip
		int i1 = win1->getLeftR();
		int j1 = win1->getLeftC();
		int i2 = win2->getLeftR();
		int j2 = win2->getLeftC();
     
		string i1ID = atomVec[i1]->getSegID();
		string j1ID = atomVec[j1]->getSegID();
		string i2ID = _distMat.getAtomVector()[i2]->getSegID();
		string j2ID = _distMat.getAtomVector()[j2]->getSegID();
    
		if(i1ID=="" || j1ID=="" || i2ID=="" ||j2ID==""){
   
		    i1ID = atomVec[i1]->getChainId();
		    j1ID = atomVec[j1]->getChainId();
		    i2ID = _distMat.getAtomVector()[i2]->getChainId();
		    j2ID = _distMat.getAtomVector()[j2]->getChainId();
		}//end if
		
		
		// Skip if both chains are forbidden within a matrix
		if (forbiddenIDMat1.find(i1ID+":"+j1ID)!=forbiddenIDMat1.end() || forbiddenIDMat2.find(i2ID+":"+j2ID)!=forbiddenIDMat2.end()) continue;
		

		// Skip if both inter-matrix segids found and if difference in residue number is not the same.
		findRegistry = properRegisterRow.find(i1ID+":"+i2ID);
		double diffInResidueNumber = atomVec[i1]->getResidueNumber() - _distMat.getAtomVector()[i2]->getResidueNumber();
		if (findRegistry != properRegisterRow.end() && findRegistry->second != diffInResidueNumber) continue;

		findRegistry        = properRegisterCol.find(j1ID+":"+j2ID);
		diffInResidueNumber = atomVec[j1]->getResidueNumber() - _distMat.getAtomVector()[j2]->getResidueNumber();
		if (findRegistry != properRegisterCol.end() && findRegistry->second != diffInResidueNumber) continue;

		
		double likeness;
		//decides which compare method from MatrixWindow to call
		switch(choice){
		    case standard:
			likeness = (*win1).compare((*win2));
			break;
		    case diag:
			likeness = (*win1).compareDiagonal((*win2));
			break;
		    case doubleDiag:
			likeness = (*win1).compareDoubleDiagonal((*win2));
			break;
		    case minDist:
			likeness = (*win1).compareMinDist((*win2));
			break;
		    case minDistRow:
			likeness=(*win1).compareMinRow((*win2));
			break;
		    case minDistCol:
			likeness=(*win1).compareMinCol((*win2));
			break;
		    default:
			cout<<"Invalid argument in function DistanceMatrix::compareAll()."<<endl;
			exit(333);
		}//end switch

		if (likeness<minLikeness && abs(likeness - minLikeness) > 0.001){
		    minLikeness = likeness;
		    minWindow1 = win1;
		    minWindow2 = win2;
		    minIndex1 = i;
		    minIndex2 = j;
		    
		    
		    //set the ID to avoid in the future
		    i1IDWin = i1ID;
		    j1IDWin = j1ID;
		    i2IDWin = i2ID;
		    j2IDWin = j2ID;
		    
		 
		}//endif
		win2=NULL;
	    }//end for on j
	    win1= NULL;
	}//end for on i
	

	// Allowed   inter matrix identifier[SEGID] difference in residue numbers
	vector<int> residueNumbers1 = minWindow1->getUpLeftResidueNumbers();
	vector<int> residueNumbers2 = minWindow2->getUpLeftResidueNumbers();

	properRegisterRow[i1IDWin+":"+i2IDWin] =   (residueNumbers1[0] - residueNumbers2[0]);
	properRegisterRow[i2IDWin+":"+i1IDWin] = - (residueNumbers1[0] - residueNumbers2[0]);
	properRegisterCol[j1IDWin+":"+j2IDWin] =   (residueNumbers1[1] - residueNumbers2[1]);
	properRegisterCol[j2IDWin+":"+j1IDWin] = - (residueNumbers1[1] - residueNumbers2[1]);

	// Forbidden intra matrix identifier[SEGID] pairs
	forbiddenIDMat1[i1IDWin+":"+j1IDWin] = false;
	forbiddenIDMat1[j1IDWin+":"+i1IDWin] = false;	
	forbiddenIDMat2[i2IDWin+":"+j2IDWin] = false;	
	forbiddenIDMat2[j2IDWin+":"+i2IDWin] = false;	


	DistanceMatrixResult currentResult(*this, *minWindow1, _distMat, *minWindow2, minLikeness);

	returnVec.push_back(currentResult);

	minWindow1=NULL;
	minWindow2=NULL;

    }//end for on k

    return returnVec;

}



pair<int, double> DistanceMatrix::compareAllWindows(vector<double> minDistVec1, vector<double> minDistVec2, vector<string> allowedAminoAcids1,vector<string> allowedAminoAcids2, int choice){
    int length = listMW.size();
    int bestIndex;
    MatrixWindow *bestWindow = NULL;
    double minLikeness = 1000000;//for choice == minRowWeighted or choice==minRowAndColWeighted
    double maxLikeness = -1.0;//for choice == correlation

    double deviationFromMaxLikeness = 0.5;

    vector<pair<double,int> > results;
    
    for (int i=0; i<length; i++){
	MatrixWindow *win = listMW[i];


        //cout << "Testing window: "<<endl<<win->toString()<<endl;
	double likeness;
	switch(choice){
		case minRowWeighted:
		    likeness = win->compareMinRowWeighted(minDistVec1);
		    if(likeness<(minLikeness+deviationFromMaxLikeness) && abs(likeness - minLikeness) > 0.001){
			minLikeness=likeness;
			bestWindow=win;
			bestIndex = i;

			results.push_back(pair<double,int>(minLikeness,bestIndex));
			}//end if
		    break;
		case minRowAndColWeighted:
		    likeness = win-> compareMinRCWeighted(minDistVec1);
		    if(likeness<(minLikeness+deviationFromMaxLikeness) && abs(likeness - minLikeness) > 0.001){
			minLikeness=likeness;
			bestWindow=win;
			bestIndex = i;
			results.push_back(pair<double,int>(minLikeness,bestIndex));
			}//end if
		    break;
		case correlation:
		    likeness = win->compareCorrelation(minDistVec1);
		    if(likeness>(maxLikeness-deviationFromMaxLikeness) && abs(likeness-maxLikeness) > 0.001){
			maxLikeness = likeness;
			bestWindow = win;
			bestIndex = i;
			results.push_back(pair<double,int>(maxLikeness,bestIndex));
		    }
		    break;
		case correlationRowCol:
			if (win->getUpLeftResidueNumbers()[0] == win->getUpLeftResidueNumbers()[1]){
				likeness = win -> compareCorrelationRowColAverage(minDistVec1);
				if(likeness>(maxLikeness-deviationFromMaxLikeness) && abs(likeness-maxLikeness) > 0.001){
					maxLikeness = likeness;
					bestWindow = win;
					bestIndex = i;

					results.push_back(pair<double,int>(maxLikeness,bestIndex));

				}
			}
		    break;
	        case correlationHeteroRowCol:
 		    likeness = win->compareCorrelationHeteroRowCol(minDistVec1, minDistVec2);
		    if(likeness>(maxLikeness-deviationFromMaxLikeness) && abs(likeness-maxLikeness) > 0.001){
			maxLikeness = likeness;
			bestWindow = win;
			bestIndex = i;

			results.push_back(pair<double,int>(maxLikeness,bestIndex));
		    }
		    break;
			
		default:
		    cout<< "Invalid int (choice) in function compareAllMinRowWeighted in DistanceMatrix."<<endl;
		    exit(336);
	}//end switch
	
	//fprintf(stdout, "\t%8.3f\n", likeness);
	win=NULL;
    }//end for on i

    pair<int, double> indexAndLikeness (0, MslTools::doubleMax);
    if (results.size() == 0){
	    return indexAndLikeness;
    }
    // Sort the results
    if(choice==minRowWeighted || choice==minRowAndColWeighted){
	    sort(results.begin(), results.end());
    } else{
	    sort(results.begin(), results.end(),greater<pair<double,int> >());
    }

		    
    // Apply Sequence Filtering to results.
    // Instead of results.size(), do results[0].first - 0.1 or something.
    int bestResultIndex = 0;
    if (allowedAminoAcids1.size() == listMW[0]->getWinSize() && allowedAminoAcids2.size() == listMW[0]->getWinSize()) {
	    for (uint i  = 0; i < results.size();i++){
		    if (listMW[results[i].second]->doesSequenceMatchRow(allowedAminoAcids1) && listMW[results[i].second]->doesSequenceMatchCol(allowedAminoAcids2)){
			    bestResultIndex = i;
			    break;
		    }
	    }
    } else  if (allowedAminoAcids1.size() == listMW[0]->getWinSize()){
	    
	    for (uint i  = 0; i < results.size();i++){
		    if (listMW[results[i].second]->doesSequenceMatchRow(allowedAminoAcids1)){
			    bestResultIndex = i;
			    break;
		    }
	    }
    } else  if (allowedAminoAcids2.size() == listMW[0]->getWinSize()){
	    for (uint i  = 0; i < results.size();i++){
		    if (listMW[results[i].second]->doesSequenceMatchCol(allowedAminoAcids2)){
			    bestResultIndex = i;
			    break;
		    }
	    }
    }



    indexAndLikeness.first  = results[bestResultIndex].second;
    indexAndLikeness.second = results[bestResultIndex].first;


    bestWindow=NULL;
    return indexAndLikeness;


}



void DistanceMatrix::printCompareInfo(vector<double> _inputVec1, vector<double> _inputVec2, pair<int, double> _result, int choice){
    //retrieve values from the pair
    int index = _result.first;
    double minLikeness = _result.second;

   
    //retrieve the winning Matrix Window from listMW
    MatrixWindow *mw = listMW[index];

    //retrieve and print relevent info
    int i = mw->getLeftR();
    int j = mw->getLeftC();
   
    string iID = atomVec[i]->getSegID().c_str();
    string jID = atomVec[j]->getSegID().c_str();
    
    if(atomVec[i]->getSegID()==""){
	iID = atomVec[i]->getChainId().c_str();
	jID = atomVec[j]->getChainId().c_str();
    }
    
    int ires = atomVec[i]->getResidueNumber();
    int jres = atomVec[j]->getResidueNumber();

    string PDBname= getFileName(PDBid);
    string PDBnameShort = PDBname.substr(0,17);

	
    cout<<"Comparing with PDB "<<PDBnameShort<<endl;

    switch(choice){
	    case minRowWeighted: 
		fprintf(stdout, "MinRowWeighted compare: \tWindow %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i, j, iID.c_str(), ires, jID.c_str(), jres, minLikeness);
		break;
	    case minRowAndColWeighted:
		fprintf(stdout, "MinRowAndColWeighted compare: \tWindow %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i, j, iID.c_str(), ires, jID.c_str(), jres, minLikeness);
		break;
	    case correlation: 
		fprintf(stdout, "Correlation compare: \tWindow %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i, j, iID.c_str(), ires, jID.c_str(), jres, minLikeness);
		break;
	    case correlationRowCol:
		fprintf(stdout, "CorrelationRowCol compare: \tWindow %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i, j, iID.c_str(), ires, jID.c_str(), jres, minLikeness);
		break;
	    case correlationHeteroRowCol:
		fprintf(stdout, "CorrelationHeteroRowCol compare: \tWindow %3d,%3d (Residues: %1s %3d, %1s %3d)\t%8.3f\n", i, j, iID.c_str(), ires, jID.c_str(), jres, minLikeness);
		break;
	    default:
		cout<<"you suck! You fucked up printCompareInfo in DistanceMatrix--the one that takes the vector. Shithead."<<endl;
		exit(3334);
    }//end switch	
	
    
    //print out the vectors we compared

    string label = "Weights:  ";
    if(weight.size() != 0){
        fprintf(stdout, "%-20s", label.c_str());
	for (uint i = 0; i < weight.size();i++){
	fprintf(stdout, "%8.3f, ", weight[i]);
    }
    fprintf(stdout, "\n");

    }
    label = "DistVectRow: ";
    fprintf(stdout, "%-20s", label.c_str());
    for (uint i = 0; i < _inputVec1.size();i++){
	fprintf(stdout, "%8.3f, ", _inputVec1[i]);
    }
    fprintf(stdout, "\n");

    //print row vector and correlation
    label = "RowsVect:";
    fprintf(stdout, "%-20s", label.c_str());
    vector<double> minRowValues = mw->getMinRowValues();

    for (int i=0; i<minRowValues.size(); i++){
	fprintf(stdout, "%8.3f, ", minRowValues[i]);
    }
    fprintf(stdout,"\n");
    if(choice != correlation){//if choice is correlation, we already printed out the correlation...
	double correlateR = MslTools::correlate(minRowValues, _inputVec1);
	fprintf(stdout, "RowCorrelation: \t %8.3f, ", correlateR);
    }
    cout<<endl;

    //print col vector and correlation if choice is MinRowAndColWeighted
    if(choice==minRowAndColWeighted || choice == correlationRowCol){
	label = "ColsVect: ";
	fprintf(stdout, "%-20s", label.c_str());

	vector<double> minColValues = mw->getMinColValues();
	for (int i=0; i<minColValues.size(); i++){
	    fprintf(stdout, "%8.3f, ", minColValues[i]);
	}
	fprintf(stdout,"\n");
	double correlateC = MslTools::correlate(minColValues, _inputVec1);
	fprintf(stdout, "ColCorrelation: \t %8.3f, ", correlateC);
	cout<<endl;
    }

    if (choice == correlationHeteroRowCol){

	    label = "DistVectCol: ";
	    fprintf(stdout, "%-20s", label.c_str());
	    for (uint i = 0; i < _inputVec2.size();i++){
		    fprintf(stdout, "%8.3f, ", _inputVec2[i]);
	    }
	    fprintf(stdout, "\n");

	    label = "ColsVect: ";
	    fprintf(stdout, "%-20s", label.c_str());
	    vector<double> minColValues = mw->getMinColValues();
	    for (int i=0; i<minColValues.size(); i++){
		    fprintf(stdout, "%8.3f, ", minColValues[i]);
	    }
	    fprintf(stdout,"\n");
	    double correlateC = MslTools::correlate(minColValues, _inputVec2);
	    fprintf(stdout, "ColCorrelation: \t %8.3f, ", correlateC);
	    cout<<endl;
    }

   

}


//must add segID
void DistanceMatrix::printCompareInfo(DistanceMatrix &_distMat, pair<vector<int>, double> _result, int choice){
    
    //retrieve values from the pair
    vector<int> mwIndex(2, 0.0);
    mwIndex= _result.first;
    double minLikeness = _result.second;

    //retrieve the winning Matrix Windows
    vector<MatrixWindow*> listMW2 = _distMat.getMatrixWindows();
    MatrixWindow *minWindow1 = listMW[mwIndex[0]];
    MatrixWindow *minWindow2 = listMW2[mwIndex[1]];

    //print information
    int i1 = (*minWindow1).getLeftR();
    int j1 = (*minWindow1).getLeftC();
    int i2 = (*minWindow2).getLeftR();
    int j2 = (*minWindow2).getLeftC();
     
    string i1ID = atomVec[i1]->getSegID().c_str();
    string j1ID = atomVec[j1]->getSegID().c_str();
    string i2ID = _distMat.getAtomVector()[i2]->getSegID().c_str();
    string j2ID = _distMat.getAtomVector()[j2]->getSegID().c_str();
    
    if(i1ID=="" || j1ID=="" || i2ID=="" ||j2ID==""){
   
	i1ID = atomVec[i1]->getChainId().c_str();
	j1ID = atomVec[j1]->getChainId().c_str();
	i2ID = _distMat.getAtomVector()[i2]->getChainId().c_str();
	j2ID = _distMat.getAtomVector()[j2]->getChainId().c_str();
    }

    int i1res = atomVec[i1]->getResidueNumber();
    int j1res = atomVec[j1]->getResidueNumber();
    int i2res = _distMat.getAtomVector()[i2]->getResidueNumber();
    int j2res = _distMat.getAtomVector()[j2]->getResidueNumber();

    string PDBname= getFileName(PDBid);
    string PDBnameShort = PDBname.substr(0,17);

    string PDBname2 = getFileName(_distMat.getPDBid());
    string PDBnameShort2 = PDBname2.substr(0,17);

    cout<<"Comparing PDBs "<<PDBnameShort<<", "<<PDBnameShort2<<endl;


    switch(choice){
	    case standard:
		fprintf(stdout, "Standard compare:\t\tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
		break;
	    case diag:
		fprintf(stdout, "Diagonal compare: \t\tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
	     	break;
	    case doubleDiag:
		fprintf(stdout, "Double Diagonal compare: \tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
		break;
	    case minDist:
		fprintf(stdout, "Minimum Distance compare: \tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
		break;
	    case minDistRow:
		fprintf(stdout, "Minimum Distance Row compare: \tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
		break;
	    case minDistCol:
		fprintf(stdout, "Minimum Distance Column compare: \tWindow1 %3d,%3d (Residues: %1s%3d, %1s%3d)\tWindow2 %3d,%3d (Residues: %1s%3d, %1s%3d)\t%8.3f\n", i1, j1, i1ID.c_str(), i1res, j1ID.c_str(), j1res, i2, j2, i2ID.c_str(), i2res, j2ID.c_str(), j2res, minLikeness);
		break;
	    default:
		cout<<"Error. Incorrect int value (choice) in DistanceMatrix::printCompareInfo(...)"<<endl;
		exit(334);
    }//end switch
}

void DistanceMatrix::createMatrixWindows(){
    int length = dimension - generalWinSize+1;

    //cout << "Dimension,generateWindowSize: "<<dimension<<","<<generalWinSize<<endl;
    for(int i=0; i<length; i++){
	for(int j=0; j<length; j++){//i and j specify the upper left corner of our window

	    // Store diagnol matrix windows
	    if (i == j){

	      
	      // Make sure diagnol window doesn't cross chain borders...and that there is no skip in residue numbering..
	      if ((atomVec[i]->getChainId() == atomVec[i+generalWinSize-1]->getChainId())){
		  //	  (atomVec[i]->getResidueNumber()+generalWinSize-1 == atomVec[i+generalWinSize-1]->getResidueNumber())){ 

		if (atomVec[i]->getResidueNumber()+generalWinSize-1 == atomVec[i+generalWinSize-1]->getResidueNumber()){
		  //		  cout << "Adding: "<<atomVec[i]->getResidueNumber()<<" and "<<atomVec[i+generalWinSize-1]->getResidueNumber()<<endl;

		  MatrixWindow *mw =  new MatrixWindow(i, j, *this, generalWinSize);
		  diagMW.push_back(mw);
		  mw = NULL;
		}
	      }
	    }

	    // Don't add any other matrix windows if this flag is set..
	    if (diagnolMatrixWindowsOnly) continue;

	    bool skip = filterWindow(i, j, generalWinSize);

	    if (!skip){

		// Create a MW.
		MatrixWindow *mw = new MatrixWindow(i, j, *this, generalWinSize);

		// Add to list of MWs
		listMW.push_back(mw);

		// Make a key for a map of MWs
		char a[80];
		sprintf(a,"%d:%d",i,j);

		// Insert MW pointer into map, using key
		mapMW[(string)a] = mw;

		// Bye-bye mw
		mw = NULL; 
		
	    }//end else

	}//end for on j
    }//end for on i
}


bool DistanceMatrix::filterWindow(int _i, int _j, int _dim){

    int i1 = _i;
    int i2 = i1 +(_dim-1);
    int j1 = _j;
    int j2 = j1 + (_dim-1);

    
    
    string i1ID =(*(atomVec[i1])).getSegID();
    string i2ID =(*(atomVec[i2])).getSegID();
    string j1ID =(*(atomVec[j1])).getSegID();
    string j2ID =(*(atomVec[j2])).getSegID();

    if(i1ID==""||i2ID==""||j1ID==""||j2ID==""){
	i1ID = (*(atomVec[i1])).getChainId();
	i2ID = (*(atomVec[i2])).getChainId();
	j1ID = (*(atomVec[j1])).getChainId();
	j2ID = (*(atomVec[j2])).getChainId();
    }

    int test1 = i2-j1;//these variables will help us skip the main diagonal
    int test2 = i1-j2;
    if(test1 != 0)
	test1= test1/fabs(test1);
    if(test2 !=0)  
	test2=test2/fabs(test2);

        //if the bool returns "false" for intraChain, this will filter out all A-A or B-B comparisons.
    if(intraChain==false && i1ID==j1ID){

	    if (debug){
		    cout << "Filtered due to intraChain matrixWindow\n";
	    }
	    return true;
    }
	
     //this section kills windows that cross either the diagnal or the vertical/horizontal lines differentiating A-A from A-B, eg.
     if((i1ID!=i2ID) || (j1ID != j2ID)){
	    if (debug){
		    cout << "Filtered due to diagnol crossing matrixWindow\n";
	    }
	    return true;
     }
/*     else if((test1 != test2) || (test1==0) || (test2==0)){//filter out main diagonal

	     if (debug){
		    cout << "Filtered due to main diagnol crossing matrixWindow\n";
	     }
	     return true;
     }
*/     else{//filters whitespace
	 int upR = _i;
         int upL = _j;
         double sum = 0;
         for(int i=0; i<_dim; i++){
	     for(int j=0; j<_dim; j++){
		 sum += (*this)[upR + i][upL + j];
	     }
	 }
	 double average = sum/(_dim*_dim);
	 if(average > 15){
	    if (debug){
		    cout << "Filtered due to 'noise' matrixWindow (15 Angstrom average)\n";
	    }
	     return true;
	 }
	 else{
	     return false;
	 }
     }//end else

}


MatrixWindow * DistanceMatrix::operator()(string _key){

    map<string,MatrixWindow *>::iterator lookup;

    lookup = mapMW.find(_key);

    if (lookup != mapMW.end()){
	return lookup->second;
    } 

    return NULL;
}


vector<MatrixWindow*> & DistanceMatrix::getDiagnolMWs(){
  return diagMW;
}


