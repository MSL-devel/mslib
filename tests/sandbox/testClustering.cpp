/*
 */


#include "Clustering.h"
#include "MslOut.h"
using namespace MSL;
using namespace std;

static MslOut MSLOUT("testClustering");
int main(){

  MSLOUT.turnAllOn();

  vector<vector<double> > matrix;
  vector<string> labels;

  matrix.resize(6);
  labels.resize(6);
  for (uint i = 0; i < 6; i++){
    matrix[i].resize(6);
  }

  labels[0] = "Julie";
  labels[1] = "John";
  labels[2] = "Ryan";
  labels[3] = "Bob";
  labels[4] = "Ted";
  labels[5] = "Kristi";

  cout << "LABELS DONE "<<matrix.size()<<" "<<matrix[0].size()<<endl;
  matrix[0][0] = 0.0;
  matrix[0][1] = 4.0;
  matrix[0][2] = 36.0;
  matrix[0][3] = 81.0;
  matrix[0][4] = 196.0;
  matrix[0][5] = 225.0;
    cout << "MATRIXXXDONE"<<endl;
  matrix[1][0] = 4.0;
  matrix[1][1] = 0.0;
  matrix[1][2] = 16.0;
  matrix[1][3] = 49.0;
  matrix[1][4] = 144.0;
  matrix[1][5] = 169.0;

  matrix[2][0] = 36.0;
  matrix[2][1] = 16.0;
  matrix[2][2] = 0.0;
  matrix[2][3] = 9.0;
  matrix[2][4] = 64.0;
  matrix[2][5] = 81.0;

  matrix[3][0] = 81.0;
  matrix[3][1] = 49.0;
  matrix[3][2] = 9.0;
  matrix[3][3] = 0.0;
  matrix[3][4] = 25.0;
  matrix[3][5] = 36.0;    

  matrix[4][0] = 196.0;
  matrix[4][1] = 144.0;
  matrix[4][2] = 64.0;
  matrix[4][3] = 25.0;
  matrix[4][4] = 0.0;
  matrix[4][5] = 1.0;    
  cout << "MATRIX1 DONE"<<endl;
  matrix[5][0] = 225.0;
  matrix[5][1] = 169.0;
  matrix[5][2] = 81.0;
  matrix[5][3] = 36.0;
  matrix[5][4] = 1.0;
  matrix[5][5] = 0.0;      

  cout << "MATRIX2 DONE"<<endl;
  Clustering clust(&matrix,&labels);
  cout << "START CLUSTERING"<<endl;
  
  clust.SingleLinkage();
}
