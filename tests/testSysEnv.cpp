#include "SysEnv.h"
using namespace MSL;
using namespace std;
#include <iostream>

static SysEnv SYSENV;

/*
  Run this program twice to test. First just run with no variables set. Then set MSLHOME to a random value and see that the output changes.
 */

int main(){

	cout << "MSLHOME is set to: "<<SYSENV.getEnv("MSLHOME")<<endl;

	cout << "HBONDPARDIR is set to: "<<SYSENV.getEnv("HBONDPARDIR")<<endl;
	cout << "HBONDPAR is set to: "<<SYSENV.getEnv("HBONDPAR")   	<<endl;
                  
	cout << "CHARMMDIR is set to: "<<SYSENV.getEnv("CHARMMDIR")  <<endl;
	cout << "CHARMMPAR is set to: "<<SYSENV.getEnv("CHARMMPAR")  <<endl;
	cout << "CHARMMTOP is set to: "<<SYSENV.getEnv("CHARMMTOP")  <<endl;
                  
	cout << "ROTLIBDIR is set to: "<<SYSENV.getEnv("ROTLIBDIR")  <<endl;
	cout << "ROTLIB is set to: "<<SYSENV.getEnv("ROTLIB")     <<endl;

	cout << "FOOBAR is set to: "<<SYSENV.getEnv("FOOBAR") <<endl;
}
