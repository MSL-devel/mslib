#include "SysEnv.h"
#include <iostream>

using namespace MSL;
using namespace std;

static SysEnv SYSENV;

/*
  Run this program twice to test. First just run with no variables set. Then set MSL_DIR to a random value and see that the output changes.
 */

int main(){

	cout << "MSL_DIR is set to: "<<SYSENV.getEnv("MSL_DIR")<<endl;

	cout << "MSL_CHARMM_PAR is set to: "<<SYSENV.getEnv("MSL_CHARMM_PAR")  <<endl;
	cout << "MSL_CHARMM_TOP is set to: "<<SYSENV.getEnv("MSL_CHARMM_TOP")  <<endl;
                  
	cout << "MSL_HBOND_PAR is set to: "<<SYSENV.getEnv("MSL_HBOND_PAR")   	<<endl;
	cout << "MSL_HBOND_CA_PAR is set to: "<<SYSENV.getEnv("MSL_HBOND_CA_PAR")   	<<endl;
                  
	cout << "MSL_PDB_TOP is set to: "<<SYSENV.getEnv("MSL_PDB_TOP")   	<<endl;
                  
	cout << "MSL_ROTLIB is set to: "<<SYSENV.getEnv("MSL_ROTLIB")     <<endl;

	cout << "MSL_EXAMPLE_FILE_DIR is set to: "<<SYSENV.getEnv("MSL_EXAMPLE_FILE_DIR")     <<endl;

	cout << "FOOBAR is set to: "<<SYSENV.getEnv("FOOBAR") <<endl;
}
