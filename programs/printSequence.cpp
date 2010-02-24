#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "Residue.h"

using namespace std;
using namespace MSL;
using namespace MslTools;

#include "printSequence.h"

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	System sys;
	sys.readPdb(opt.pdb);

	for (uint c = 0; c < sys.size();c++){
	  for (uint r = 0; r < sys.getChain(c).size();r++){
	      Residue &res = sys.getChain(c).getResidueByIndex(r);
	      fprintf(stdout,"%1s",MslTools::getOneLetterCode(res.getResidueName()).c_str());
	  }
      fprintf(stdout,"\n");
	}

}


Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option


	if (OP.countOptions() == 0){
		cout << "Usage: printSequence " << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}

	return opt;
}
