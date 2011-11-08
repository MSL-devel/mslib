#include "DistanceMatrixDatabase.h"
#include "MatrixWindow.h"

using namespace MSL;

DistanceMatrixDatabase::DistanceMatrixDatabase(){
	name = "Generic";
	archiveType = "binary";
}
DistanceMatrixDatabase::DistanceMatrixDatabase(const DistanceMatrixDatabase &_dm){
	name = _dm.name;
	archiveType = _dm.archiveType;
}
DistanceMatrixDatabase::~DistanceMatrixDatabase(){
	for (uint i = 0;i < listDMs.size();i++){
		if (listDMs[i] != NULL){
			delete(listDMs[i]);
		}
	}
	listDMs.clear();
}
