#ifndef MANAGERESULTS_H
#define MANAGERESULTS_H

#include "DistanceMatrixResult.h"
#include "MatrixWindow.h"

#include<iostream>

using namespace std;
namespace MSL{
class ManageDistanceMatrixResults{
  public:
    ManageDistanceMatrixResults();
    ~ManageDistanceMatrixResults();

    void addResults(vector<DistanceMatrixResult> &_inputVec){allResults.push_back(_inputVec);}

    void printResults();

    void setAlignPdbs(bool _flag) { alignPdbs = _flag; }
    bool getAlignPdbs() { return alignPdbs; }

    void setRmsdTol(double _tol) { rmsdTol = _tol; }
    double getRmsdTol() { return rmsdTol; }
    
  private:
    bool alignPdbs;
    double rmsdTol;
    vector<vector<DistanceMatrixResult> > allResults; //allResults[0] gives a vector of the results from one multicompare between two DM's



};
}

#endif
