#ifndef RESULT_H
#define RESULT_H

#include "DistanceMatrix.h"
#include "MatrixWindow.h"


using namespace std;
namespace MSL {

class DistanceMatrixResult{
  public:
	DistanceMatrixResult();
	DistanceMatrixResult(DistanceMatrix &_dm1, MatrixWindow &_mw1, DistanceMatrix &_dm2, MatrixWindow &_mw2, double _likeness);
	DistanceMatrixResult(const DistanceMatrixResult &_result);
	~DistanceMatrixResult();

	void operator=(const DistanceMatrixResult &_result) { copy(_result); }

	//inline accesor functions
	DistanceMatrix & getDistanceMatrix1() const {return *dm1;}
	DistanceMatrix & getDistanceMatrix2() const {return *dm2;}
	MatrixWindow & getMatrixWindow1()const {return *mw1;}
	MatrixWindow & getMatrixWindow2() const {return *mw2;}
	double getLikeness() const {return likeness;}

	bool operator > (const DistanceMatrixResult &_leftHandSide) const;
	bool operator < (const DistanceMatrixResult &_leftHandSide) const;


	string toString();
 private:
	void copy(const DistanceMatrixResult &_result);
    
	DistanceMatrix *dm1;
	DistanceMatrix *dm2;
	MatrixWindow *mw1;
	MatrixWindow *mw2;
	double likeness;

};
}

#endif
