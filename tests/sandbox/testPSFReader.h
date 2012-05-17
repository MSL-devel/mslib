#include "PSFReader.h"



int main(){

	PSFReader psf;
	psf.open("/home/dwkulp/sofwtare/msl/bpti.psf");
	psf.read();
	psf.close();
	
	return 1;
}
