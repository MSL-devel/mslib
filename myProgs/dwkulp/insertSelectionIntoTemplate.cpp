#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "AtomContainer.h"
#include "insertSelectionIntoTemplate.h"

using namespace MSL;
using namespace std;
using namespace MslTools;

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	// Template is the structure to insert into
	System templatePDB;
	templatePDB.readPdb(opt.templatePDB);

	// Fragment is the structure to take peice of structure from
	System fragmentPDB;
	fragmentPDB.readPdb(opt.fragmentPDB);

	// Discover stems by distance
	int stem1index = -1;
	int stem2index = -1;
	if (opt.templateStem1 == "" || opt.templateStem2 == ""){

	  Residue &frag1 = fragmentPDB.getPosition(0).getCurrentIdentity();
	  for (uint j = 0; j < templatePDB.positionSize();j++){
	    Residue &template1 = templatePDB.getPosition(j).getCurrentIdentity();
	    if (frag1.atomExists("CA") && template1.atomExists("CA") &&
		  frag1.getLastFoundAtom().distance(template1.getLastFoundAtom()) < 0.3){
	            stem1index = j;
		    break;
	      }
	  }

	  Residue &frag2 = fragmentPDB.getPosition(fragmentPDB.positionSize()-1).getCurrentIdentity();
	  for (uint j = 0; j < templatePDB.positionSize();j++){
	    Residue &template1 = templatePDB.getPosition(j).getCurrentIdentity();
	    if (frag2.atomExists("CA") && template1.atomExists("CA") &&
		  frag2.getLastFoundAtom().distance(template1.getLastFoundAtom()) < 0.3){
	            stem2index = j;
		    break;
	      }
	  }
	  
	} else {
	  stem1index = templatePDB.getPositionIndex(opt.templateStem1);
	  stem2index = templatePDB.getPositionIndex(opt.templateStem2);
	}
	Position &stem1 = templatePDB.getPosition(stem1index);
	Position &stem2 = templatePDB.getPosition(stem2index);

	fprintf(stdout, "Stems: %s , %s\n",stem1.toString().c_str(),stem2.toString().c_str());

	cout << "Add upto stem1"<<endl;
	AtomContainer fusedProtein;
	for (uint i = 0; i < stem1.getIndexInSystem();i++){
	  fusedProtein.addAtoms(templatePDB.getPosition(i).getAtomPointers());
	}




	cout << "Inserting fragment"<<endl;
	int newResidueNumber = stem1.getResidueNumber();
	bool stem1added = false;
	for (uint i = 0; i < fragmentPDB.positionSize();i++){
	  Residue &fragRes = fragmentPDB.getPosition(i).getCurrentIdentity();
	  
	  if (fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem1.getCurrentIdentity()("CA")) < 0.3) {
	       fprintf(stdout, "FUSION Template Stem1 to %s\n",fragRes.toString().c_str());

	       AtomContainer newStemRes;
	       newStemRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newStemRes.size();a++){
		 newStemRes[a].setChainId(stem1.getChainId());
		 newStemRes[a].setResidueNumber(newResidueNumber);
		 newStemRes[a].setResidueIcode(stem1.getResidueIcode());
	       }
	       
	       fusedProtein.addAtoms(newStemRes.getAtomPointers());

	       stem1added = true;
	       
	       newResidueNumber++;
	       continue;
	  }





	  if (stem1added && fragRes.atomExists("CA") && fragRes.getLastFoundAtom().distance(stem2.getCurrentIdentity()("CA")) < 0.3) {
	       fprintf(stdout, "FUSION Template Stem2 to %s\n",fragRes.toString().c_str());

	       AtomContainer newStemRes;
	       newStemRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newStemRes.size();a++){
		 newStemRes[a].setChainId(stem2.getChainId());
		 newStemRes[a].setResidueNumber(newResidueNumber);
	       }
	       
	       fusedProtein.addAtoms(newStemRes.getAtomPointers());
	       newResidueNumber++;
	       // We should be done by now..
	       break;
	  }
	  
	  if (stem1added) {

	       AtomContainer newFragRes;
	       newFragRes.addAtoms(fragRes.getAtomPointers());
	       for (uint a = 0; a < newFragRes.size();a++){
		 newFragRes[a].setChainId(stem1.getChainId());
		 newFragRes[a].setResidueNumber(newResidueNumber);
	       }
	       
	       fusedProtein.addAtoms(newFragRes.getAtomPointers());
	       newResidueNumber++;
	  }

	}

	cout << "Add rest of protein"<<endl;

	// Add the rest of the template protein
	for (uint i = stem2.getIndexInSystem()+1; i < templatePDB.positionSize();i++){
	  AtomContainer tmp;
	  tmp.addAtoms(templatePDB.getPosition(i).getAtomPointers());
	  for (uint t = 0; t< tmp.size();t++){
	    tmp[t].setResidueNumber(newResidueNumber);
	  }
	  fusedProtein.addAtoms(tmp.getAtomPointers());
	  newResidueNumber++;
	}
		

	char fname[200];
	sprintf(fname, "%s_%s.pdb",MslTools::getFileName(opt.templatePDB).c_str(),MslTools::getFileName(opt.fragmentPDB).c_str());
	cout << "Fused Protein: "<< fname <<endl;

	fusedProtein.writePdb(fname);

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
		cout << "Usage: insertSelectionIntoTemplate " << endl;
		cout << endl;
		cout << "\n";
		cout << "template PDB\n";
		cout << "templateStem1 positionId\n";
		cout << "templateStem2 positionId\n";
		cout << "fragment PDB\n";
		cout << endl;
		exit(0);
	}

	opt.templatePDB = OP.getString("template");
	if (OP.fail()){
		cerr << "ERROR 1111 no template specified."<<endl;
		exit(1111);
	}


	opt.fragmentPDB = OP.getString("fragment");
	if (OP.fail()){
		cerr << "ERROR 1111 no fragment specified."<<endl;
		exit(1111);
	}


	opt.templateStem1= OP.getString("templateStem1");
	opt.templateStem2 = OP.getString("templateStem2");


	return opt;
}
