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

using namespace std;
using namespace MSL;
using namespace MslTools;


int main(int argc, char *argv[]){

    DistanceMatrixDatabase dmd;
    dmd.load_checkpoint("degrado.bin");

    vector<DistanceMatrix *> &dms = dmd.getDistanceMatrixList();

    for (uint i = 0;i < dms.size();i++){
	    cout << "DM["<<i<<"]: "<<dms[i]->getPDBid()<<endl;
    }
    cout << "Done"<<endl;
}
