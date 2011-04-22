#include "MslOut.h"
#include "Transforms.h"
#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "OptionParser.h"
#include "Timer.h"
#include <stdio.h>
using namespace std;
using namespace MSL;

static MslOut MSLOUT("testMSLOut2");
int main(int argc, char *argv[]){

	OptionParser OP;
	OP.readArgv(argc, argv);

	Transforms t;
	AtomPointerVector av;
	MSLOUT.stream(MslOut::GENERAL) << "Its a general statement, this will always output."<<std::endl;
	MSLOUT.stream(MslOut::SPECIFIC) << "Its a specific statement, this will only output when this object output is on"<<std::endl;
	MSLOUT.stream(MslOut::WARNING) << "Its a WARNING that is specific to this object"<<std::endl;
	MSLOUT.stream(MslOut::ERROR) << "Its an ERROR that is specific to this object"<<std::endl;
	MSLOUT.debug() << "A debug statement like this is only seen when __MSL_MSLOUT_DEBUG_OFF__ is unset!"<<std::endl;

	// Performance test for MSLOUT;
	Timer tme;
	double start = tme.getWallTime();
	AtomContainer a;
	for (uint i = 0; i < 100000;i++){
		Atom *b = new Atom();
		a.addAtom(*b);
	}
	double end = tme.getWallTime();
	fprintf(stdout, "Time: %8.3f\n",(end-start));
	
}

