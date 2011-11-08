#include "MatrixWindow.h"

using namespace MSL;

//constructor of MatrixWindow
MatrixWindow::MatrixWindow(){
	upLeftR = 0;
	upLeftC = 0;
	bigMatrix = NULL;
	winSize = 0;
}
MatrixWindow::MatrixWindow(int initLeftR, int initLeftC, DistanceMatrix &m, int _winSize){
	bigMatrix = &m;
	upLeftR = initLeftR;
	upLeftC = initLeftC;
	winSize = _winSize;

	//sets the small atom vector
	setSmallAVec();

	//sets min distance
	setMinDistanceRow();
	setMinDistanceCol();


	archiveType = "binary";
}

//destructor of MatrixWindow
MatrixWindow::~MatrixWindow(){
	bigMatrix = NULL;
}

MatrixWindow::MatrixWindow(const MatrixWindow&_mw){
	copy(_mw);
}

void MatrixWindow::copy(const MatrixWindow &_mw){
	upLeftR = _mw.upLeftR;
	upLeftC = _mw.upLeftC;
	winSize = _mw.winSize;
	bigMatrix = _mw.bigMatrix;
	smallAtomVec = _mw.smallAtomVec;
	minDistanceRow = _mw.minDistanceRow;
	minDistanceCol = _mw.minDistanceCol;
	minRowValues = _mw.minRowValues;
	minColValues = _mw.minColValues;
	archiveType = _mw.archiveType;

}


//makes the small atomVector
void MatrixWindow::setSmallAVec(){
	for (int i=0; i<winSize; i++){
		smallAtomVec.push_back((bigMatrix->getAtomVector())[i+upLeftR]);
	}//end for on i
	for (int j=0; j<winSize; j++){
		smallAtomVec.push_back((bigMatrix->getAtomVector())[j+upLeftC]);
	}//end for on j
}

void MatrixWindow::setMinDistanceRow(){

	for(int i=0; i<winSize; i++){
		double minValue=(*bigMatrix)[upLeftR+i][upLeftC];
		int minCol = upLeftC;
		for(int j=0; j<winSize; j++){
	    
			if((*bigMatrix)[upLeftR+i][upLeftC+j] < minValue){
				minValue = (*bigMatrix)[upLeftR+i][upLeftC+j];
				minCol = upLeftC+j;
			}//end if

		}//end for on j
		minDistanceRow.push_back(minCol);
		minRowValues.push_back(minValue);
	}//end for on i
}


void MatrixWindow::setMinDistanceCol(){
	for(int i=0; i<winSize; i++){
		double minValue = (*bigMatrix)[upLeftR][upLeftC+i];
		int minRow = upLeftR;
		for(int j=0; j<winSize; j++){
			if((*bigMatrix)[upLeftR+j][upLeftC+i] < minValue){
				minValue = (*bigMatrix)[upLeftR+j][upLeftC+i];
				minRow=upLeftR+j;
			}
		}//end for on j
		minDistanceCol.push_back(minRow);
		minColValues.push_back(minValue);
	}//end for on i
}

double MatrixWindow::compare(MatrixWindow & mw){
	DistanceMatrix *thisMatrix= bigMatrix;
	DistanceMatrix *otherMatrix = &mw.getMatrix();
	double likeness=0; 
	for(int i=0; i<winSize; i++){
		for(int j=0; j<winSize; j++){
			double diff = (*thisMatrix)[getLeftR()+i][getLeftC()+j]-(*otherMatrix)[mw.getLeftR()+i][mw.getLeftC()+j];
			likeness = likeness + diff*diff;
		}//end for on j
	}//end for on i
	return likeness;
}

double MatrixWindow::compareDiagonal(MatrixWindow & _mw){
	DistanceMatrix *thisMatrix=bigMatrix;
	DistanceMatrix *otherMatrix=&_mw.getMatrix();
	double likeness=0;
	for(int i=0; i<winSize; i++){
		double diff = (*thisMatrix)[getLeftR()+i][getLeftC()+i]-(*otherMatrix)[_mw.getLeftR()+i][_mw.getLeftC()+i];
		likeness =likeness + diff*diff;
	}//end for on i
	return likeness;
}

double MatrixWindow::compareDoubleDiagonal(MatrixWindow &_mw){
	DistanceMatrix *thisMatrix=bigMatrix;
	DistanceMatrix *otherMatrix=&_mw.getMatrix();
	double likeness=0;
	for(int i=0; i<winSize; i++){
		double diff1 = (*thisMatrix)[getLeftR()+i][getLeftC()+i]-(*otherMatrix)[_mw.getLeftR()+i][_mw.getLeftC()+i];
		double diff2 = (*thisMatrix)[getLeftR()+i][getLeftC()+(winSize-1)-i]-(*otherMatrix)[_mw.getLeftR()+i][_mw.getLeftC()+(winSize-1)-i];
		likeness =likeness + diff1*diff1 +diff2*diff2;
	}//end for on i
	return likeness;
}

double MatrixWindow::compareMinDist(MatrixWindow &_mw){
	DistanceMatrix *thisMatrix=bigMatrix;
	DistanceMatrix *otherMatrix=&_mw.getMatrix();
	double likeness=0;
	for(int i=0; i<winSize; i++){//compares min distance along each row
		int row1 = upLeftR+i;
		int row2 = _mw.getLeftR()+i;
		int col1 = minDistanceRow[i];
		int col2 = _mw.getMinDistanceRow()[i];

		double diff = (*thisMatrix)[row1][col1]-(*otherMatrix)[row2][col2];
		likeness += diff*diff;
	}//end for on i

	double likeness2=0;
	for(int i=0; i<winSize; i++){//compares min distance along each column
		int row1 = minDistanceCol[i];
		int row2 = _mw.getMinDistanceCol()[i];
		int col1 = upLeftC+i;
		int col2 = _mw.getLeftC()+i;

		double diff = (*thisMatrix)[row1][col1]-(*otherMatrix)[row2][col2];
		likeness2 += diff*diff;
	}//end for on i
	double returnValue=likeness+likeness2;
	return returnValue;
}

double MatrixWindow::compareMinRow(MatrixWindow &_mw){
	DistanceMatrix *thisMatrix=bigMatrix;
	DistanceMatrix *otherMatrix=&_mw.getMatrix();
	double likeness = 0;
	for(int i=0; i<winSize; i++){
		int row1 = upLeftR+i;
		int row2 = _mw.getLeftR()+i;
		int col1 = minDistanceRow[i];
		int col2 = _mw.getMinDistanceRow()[i];

		double diff = (*thisMatrix)[row1][col1]-(*otherMatrix)[row2][col2];
		likeness += diff*diff;
	}//end for on i
	return likeness;
}

double MatrixWindow::compareMinCol(MatrixWindow &_mw){
	DistanceMatrix *thisMatrix=bigMatrix;
	DistanceMatrix *otherMatrix=&_mw.getMatrix();
	double likeness =0;
	for(int i=0; i<winSize; i++){
		int row1 = minDistanceCol[i];
		int row2 = _mw.getMinDistanceCol()[i];
		int col1 = upLeftC+i;
		int col2 = _mw.getLeftC()+i;

		double diff = (*thisMatrix)[row1][col1]-(*otherMatrix)[row2][col2];
		likeness += diff*diff;
	}//end for on i
	return likeness;
}

double MatrixWindow::compareMinRowWeighted(vector<double> _minDistVec){
	DistanceMatrix *thisMatrix=bigMatrix;
	double likeness=0;
	double sumDm=0;

	for (int i=0; i<_minDistVec.size(); i++){
		int row = upLeftR +i;
		int col = minDistanceRow[i];

		double dm = _minDistVec[i];
		double dn = (*thisMatrix)[row][col];

		double factor = 1/dm; //default
		vector<double> weightVec = bigMatrix->getWeight();
		if(weightVec.size() != 0){//if we have specified weights
			factor = weightVec[i];
		}
	

		likeness += (dm-dn)*(dm-dn)*factor;

		sumDm += factor;
	}//end for on i
	double normalizedLikeness = likeness/sumDm;

	return normalizedLikeness;
}

double MatrixWindow::compareMinRCWeighted(vector<double> _minDistVec){
	DistanceMatrix *thisMatrix=bigMatrix;
	double likeness = 0;
	double sumDm=0;

	//min row calculations
	for (int i=0; i<_minDistVec.size(); i++){
		int row = upLeftR +i;
		int col = minDistanceRow[i];

		double dm = _minDistVec[i];
		double dn = (*thisMatrix)[row][col];
	
		double factor = 1/dm; //default
		vector<double> weightVec = bigMatrix->getWeight();
		if(weightVec.size() != 0){//if we have specified weights
			factor = weightVec[i];
		}

		likeness += (dm-dn)*(dm-dn)*factor;
		sumDm += factor;
	}//end for on i
   
	//min col calculations
	for (int i=0; i<_minDistVec.size(); i++){
		int row = minDistanceCol[i];
		int col = upLeftC + i;

		double dm = _minDistVec[i];
		double dn = (*thisMatrix)[row][col];

		double factor = 1/dm; //default
		vector<double> weightVec = bigMatrix->getWeight();
		if(weightVec.size() != 0){//if we have specified weights
			factor = weightVec[i];
		}

		likeness += (dm-dn)*(dm-dn)*factor;
		sumDm += factor;
	}//end for on i

	double normalizedLikeness = likeness/sumDm;
	return normalizedLikeness;
}

double MatrixWindow::compareCorrelation(vector<double> _minDistVec){
	vector<double> weightVec = bigMatrix->getWeight();

	vector<double> filtedMinDistVec;
	vector<double> filtedMinRowValues;
	for (uint i = 0; i < _minDistVec.size();i++){
		if  (weightVec[i] != 0.0){
			filtedMinDistVec.push_back(_minDistVec[i]);
			filtedMinRowValues.push_back(minRowValues[i]);
		}
	}
	
	//double likeness = MslTools::correlate(minRowValues, _minDistVec);
	double likeness = MslTools::correlate(filtedMinRowValues,filtedMinDistVec);
	return likeness;
}

double MatrixWindow::compareCorrelationRowCol(vector<double> _minDistVec){
	vector<double> bigMinDistVec(2*_minDistVec.size(), 0.0);
	vector<double> bigRCVec(2*_minDistVec.size(), 0.0);

	//note: minRowValues.size()=minColValues.size()=_minDistVec.size() or else something is seriously wrong.
	for(int i=0; i<_minDistVec.size(); i++){
		bigMinDistVec[i]                     = _minDistVec[i];
		bigMinDistVec[_minDistVec.size() +i] = _minDistVec[i];
		bigRCVec[i]                          = minRowValues[i];  
		bigRCVec[_minDistVec.size() + i]     = minColValues[i];
	}
	double likeness = MslTools::correlate(bigMinDistVec, bigRCVec);
	return likeness;
}

double MatrixWindow::compareCorrelationRowColAverage(vector<double> _minDistVec){
	vector<double> bigRCVec(_minDistVec.size(), 0.0);

	//note: minRowValues.size()=minColValues.size()=_minDistVec.size() or else something is seriously wrong.
	for(int i=0; i<_minDistVec.size(); i++){
		bigRCVec[i]  = (minRowValues[i] + minColValues[i]) / 2;  

	}
	double likeness = MslTools::correlate(_minDistVec, bigRCVec);
	return likeness;
}

double MatrixWindow::compareCorrelationHeteroRowCol(vector<double> _minDistVec1, vector<double> _minDistVec2){
	vector<double> bigMinDistVec(_minDistVec1.size()+_minDistVec2.size(), 0.0);
	vector<double> bigRCVec(_minDistVec1.size()+_minDistVec2.size(), 0.0);

	//note: minRowValues.size()=minColValues.size()=_minDistVec.size() or else something is seriously wrong.
	stringstream rowVals;
	stringstream colVals;
	stringstream vec1Vals;
	stringstream vec2Vals;
	for(int i=0; i<_minDistVec1.size(); i++){
		bigMinDistVec[i]                     = _minDistVec1[i];
		bigRCVec[i]                          = minRowValues[i];  
		rowVals << minRowValues[i] <<", ";
		vec1Vals << _minDistVec1[i] <<", ";
	}

	int offset = _minDistVec1.size()-1;
	for(int i=0; i<_minDistVec2.size(); i++){
		if (i == 1 || i == 5 || i == 9 || i == 10 || i == 11) continue;
		int index = offset+i;
		bigMinDistVec[index] = _minDistVec2[i];
		bigRCVec[index]      = minColValues[i];  

		colVals << minColValues[i] <<", ";
		vec2Vals << _minDistVec2[i] <<", ";
		
	}


	cout << rowVals.str()<<endl;
	cout << vec1Vals.str()<<endl;
	cout << "--"<<endl;
	cout << colVals.str()<<endl;
	cout << vec2Vals.str()<<endl;
	double likeness = MslTools::correlate(bigMinDistVec, bigRCVec);

	
	return likeness;
}

void MatrixWindow::initialize(int _leftR, int _leftC, DistanceMatrix &m, int winSize){
	setLeftR(_leftR);
	setLeftC(_leftC);
	setMatrix(m);
	setSmallAVec();
	setMinDistanceRow();
	setMinDistanceCol();
}

string MatrixWindow::toString(){

	
	stringstream result;

	Atom *row = bigMatrix->getAtomVector()[upLeftR];
	Atom *col = bigMatrix->getAtomVector()[upLeftC];
	result << row->toString()<<endl<<col->toString()<<endl;
// 	for(int i=0; i<winSize; i++){
// 		for(int j=0; j<winSize; j++){
// 			char a[80];
// 			sprintf(a,"%8.3f ",(*bigMatrix)[getLeftR()+i][getLeftC()+j]);
// 			result << a;
// 		}
// 		result << endl;
// 	}

	return result.str();
}



vector<int> MatrixWindow::getUpLeftResidueNumbers(){
	vector<int> tmp;
    
	tmp.push_back(bigMatrix->getAtomVector()[upLeftR]->getResidueNumber());
	tmp.push_back(bigMatrix->getAtomVector()[upLeftC]->getResidueNumber());

	return tmp;
}


bool MatrixWindow::doesSequenceMatchRow(vector<string> _allowableAminoAcids){
	

	AtomPointerVector &_av = bigMatrix->getAtomVector();

	
	bool found = false;
	for (uint i = 0 ;  i < winSize;i++){
		Atom *a = _av[upLeftR+i];

		found = false;
		vector<string> AAs = MslTools::tokenize(_allowableAminoAcids[i],":");
		for (uint j = 0; j < AAs.size();j++){
			
			//cout << "ResName: "<<a->getResidueName()<<" "<<AAs[j]<<"."<<endl;
			if (a->getResidueName() == AAs[j] || AAs[j] == "X"){
				found = true;
				break;
			}
		}

		if (!found) break;

	}


	return found;
}


bool MatrixWindow::doesSequenceMatchCol(vector<string> _allowableAminoAcids){
	

	AtomPointerVector &_av = bigMatrix->getAtomVector();

	
	bool found = false;
	for (uint i = 0 ;  i < winSize;i++){
		Atom *a = _av[upLeftC+i];

		found = false;
		vector<string> AAs = MslTools::tokenize(_allowableAminoAcids[i],":");
		for (uint j = 0; j < AAs.size();j++){
			
			if (a->getResidueName() == AAs[j] || AAs[j] == "X"){
				found = true;
				break;
			}
		}

		if (!found) break;

	}


	return found;
}
