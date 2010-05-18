#include "MslOut.h"
#include "Transforms.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include <stdio.h>
using namespace std;
using namespace MSL;

static MslOut MSLOUT("testMSLOut2");
int main(int argc, char *argv[]){

	OptionParser OP;
	OP.readArgv(argc, argv);

	Transforms t;
	AtomPointerVector av;

	MSLOUT.stream(MslOut::GENERAL) << "Hey general!"<<std::endl;
	MSLOUT.stream(MslOut::SPECIFIC) << "Hey specific!"<<std::endl;
}

