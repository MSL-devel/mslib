#include "DistanceMatrixResult.h"
#include "MslTools.h"
using namespace std;
using namespace MSL;

//constructor
DistanceMatrixResult::DistanceMatrixResult(){
	dm1 = NULL;
	dm2 = NULL;
	mw1 = NULL;
	mw2 = NULL;
	likeness = MslTools::doubleMax;
}
DistanceMatrixResult::DistanceMatrixResult(DistanceMatrix &_dm1, MatrixWindow &_mw1, DistanceMatrix &_dm2, MatrixWindow &_mw2, double _likeness){
	dm1 = &_dm1;
	dm2 = &_dm2;

//	cout << "DM1: "<<dm1<<" DM2: "<<dm2<<endl;
	mw1 = &_mw1;
	mw2 = &_mw2;
	likeness = _likeness;
}

DistanceMatrixResult::DistanceMatrixResult(const DistanceMatrixResult &_result){
	copy(_result);
}


DistanceMatrixResult::~DistanceMatrixResult(){

}


void DistanceMatrixResult::copy(const DistanceMatrixResult &_result){

	dm1 = _result.dm1;
	dm2 = _result.dm2;
	mw1 = _result.mw1;
	mw2 = _result.mw2;
	likeness = _result.likeness;
}

bool DistanceMatrixResult::operator > (const DistanceMatrixResult &_leftHandSide) const{
	if(likeness < _leftHandSide.getLikeness())
		return true;
	else
		return false;
}

bool DistanceMatrixResult::operator < (const DistanceMatrixResult &_leftHandSide) const{
	if (likeness > _leftHandSide.getLikeness())
		return true;
	else
		return false;
}


string DistanceMatrixResult::toString(){


	char str[500];
	sprintf(str, "%s %s %8.3f\n", MslTools::getFileName(dm1->getPDBid()).c_str(), MslTools::getFileName(dm2->getPDBid()).c_str(), likeness);

	return (string)(str);

}
