#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "Timer.h"
#include "Transforms.h"
#include "AtomContainer.h"
#include "AtomSelection.h"
#include "dualLoopSearch.h"
#include <omp.h>

using namespace MSL;
using namespace std;
using namespace MslTools;

int main(int argc, char *argv[]) {
	Timer t;
	double startTime = t.getWallTime();

	Options opt = setupOptions(argc, argv);


	vector<string> theLines;
	MslTools::readTextFile(theLines,opt.pdblist);
	

	int i = 0;
        #pragma omp parallel private(i) shared(opt,theLines)
	{
	  vector<string> lines;

	  #pragma omp critical
	  lines = theLines;

	 #pragma omp for schedule(dynamic)
         for (i = 0; i < lines.size();i++){
	  string fileName = lines[i];

	   int dualLoopIndex = 0;

	   System sys;
	   sys.readPdb(fileName);	   


	   System stem1;
	   stem1.readPdb(opt.stem1pdb);
	   Residue &stem1res1 = stem1.getResidue(opt.stem1res1);
	   Residue &stem1res2 = stem1.getResidue(opt.stem1res2);

	   if (! (stem1res1.atomExists("CA") && stem1res1.atomExists("N") && stem1res1.atomExists("C"))) {
	     cerr << "ERROR STEM1RES1 does not have N-CA-C atoms\n";
	     exit(1);
	   };

	   if (! (stem1res2.atomExists("CA") && stem1res2.atomExists("N") && stem1res2.atomExists("C"))) {
	     cerr << "ERROR STEM1RES2 does not have N-CA-C atoms\n";
	     exit(1);
	   };

	   System stem2;
	   stem2.readPdb(opt.stem2pdb);
	   Residue &stem2res1 = stem2.getResidue(opt.stem2res1);
	   Residue &stem2res2 = stem2.getResidue(opt.stem2res2);


	   System stem2_orig;
	   stem2_orig.readPdb(opt.stem2pdb);
	   AtomSelection stem2origsel(stem2_orig.getAtomPointers());
	   AtomPointerVector stem2_orig_ca = stem2origsel.select("name CA");
	   
	   if (! (stem2res1.atomExists("CA") && stem2res1.atomExists("N") && stem2res1.atomExists("C"))) {
	     cerr << "ERROR STEM2RES1 does not have N-CA-C atoms\n";
	     exit(1);
	   };
	   if (! (stem2res2.atomExists("CA") && stem2res2.atomExists("N") && stem2res2.atomExists("C"))) {
	     cerr << "ERROR STEM2RES2 does not have N-CA-C atoms\n";
	     exit(1);
	   };
	   cout << MslTools::getFileName(fileName)<<" # chains : "<< sys.chainSize() << " # residues: "<<sys.positionSize() <<endl;
 

	   vector<vector<int> > residueQuadsThatMatch;
	   
	   // Each chain
	   for (uint c = 0; c < sys.chainSize();c++){
	     
	     // Walk through residues
	     Chain &ch = sys.getChain(c);

	     // Skip overly large and small chains to speed up calculations
	     //if (ch.size() > 600 || ch.size() < 100) {
	     //	 continue;
	     //}

	     for (uint r1 = 0; r1 < ch.positionSize();r1++){
	       Residue &loop1Res1  = ch.getResidue(r1);

	       if (! (loop1Res1.atomExists("CA") && loop1Res1.atomExists("N") && loop1Res1.atomExists("C"))) continue;

	       for (uint r2 = r1+opt.loop1min; r2 < r1+opt.loop1max && r2 < ch.positionSize();r2++){
		 Residue &loop1Res2  = ch.getResidue(r2);
		 if (! (loop1Res2.atomExists("CA") && loop1Res2.atomExists("N") && loop1Res2.atomExists("C"))) continue;


		 for (uint r3 = r2+1; r3 < ch.positionSize();r3++){
		   if (r2 == r3) continue;
		   Residue &loop2Res1  = ch.getResidue(r3);
		   if (! (loop2Res1.atomExists("CA") && loop2Res1.atomExists("N") && loop2Res1.atomExists("C"))) continue;


		   // Look for R2-R3 distance match
		   double dist2 = loop1Res2.distance(loop2Res1,"CA",true); // true = distance squared
		   if (abs(dist2 - opt.distanceStem2) > 20) {
		     continue;
		   }


		   for (uint r4 = r3+opt.loop2min; r4 < r3+opt.loop2max && r4 < ch.positionSize();r4++){
		     if (r1 == r4 || r1 == r3) continue;

		     Residue &loop2Res2  = ch.getResidue(r4);
		     if (! (loop2Res2.atomExists("CA") && loop2Res2.atomExists("N") && loop2Res2.atomExists("C"))) continue;



		     
		     // Look for R1-R4 distance match (< 1 Angstroms)
		     double dist1 = loop1Res1.distance(loop2Res2,"CA",true); // true = distance squared
		     if (abs(dist1 - opt.distanceStem1) < 20){ 

		       // RMSD Check to stems
		       AtomContainer move;
		       move.addAtom(loop1Res1("N"));
		       move.addAtom(loop1Res1("CA"));
		       move.addAtom(loop1Res1("C"));
		       move.addAtom(loop2Res2("N"));
		       move.addAtom(loop2Res2("CA"));
		       move.addAtom(loop2Res2("C"));

		       AtomContainer dualLoopAtoms;
		       for (uint m = r1; m <= r2;m++){
			 dualLoopAtoms.addAtoms(ch.getResidue(m).getAtomPointers());
		       }

		       for (uint m = r3; m <= r4;m++){
			 dualLoopAtoms.addAtoms(ch.getResidue(m).getAtomPointers());
		       }

		       AtomContainer ref;
		       ref.addAtom(stem1res1("N"));
		       ref.addAtom(stem1res1("CA"));
		       ref.addAtom(stem1res1("C"));
		       ref.addAtom(stem1res2("N"));
		       ref.addAtom(stem1res2("CA"));
		       ref.addAtom(stem1res2("C"));	   

		       AtomContainer ref2;
		       ref2.addAtom(loop1Res2("N"));
		       ref2.addAtom(loop1Res2("CA"));
		       ref2.addAtom(loop1Res2("C"));
		       ref2.addAtom(loop2Res1("N"));
		       ref2.addAtom(loop2Res1("CA"));
		       ref2.addAtom(loop2Res1("C"));

		       Transforms t;
		       if (!t.rmsdAlignment(move.getAtomPointers(),ref.getAtomPointers(),dualLoopAtoms.getAtomPointers())) {
			 cerr << "Error doing rmsd alignment"<<endl;
			 exit(0);
		       }

		       t.rmsdAlignment(move.getAtomPointers(),ref.getAtomPointers(),ref2.getAtomPointers());
		       t.rmsdAlignment(move.getAtomPointers(),ref.getAtomPointers());

		       double rmsd = ref.getAtomPointers().rmsd(move.getAtomPointers());

		       if (rmsd > opt.rmsdTol){
			 continue;
		       }
		       

		       // Good now do the other stem
		       AtomContainer stem2copy;
		       stem2copy.addAtom(stem2res1("N"));
		       stem2copy.addAtom(stem2res1("CA"));
		       stem2copy.addAtom(stem2res1("C"));
		       stem2copy.addAtom(stem2res2("N"));
		       stem2copy.addAtom(stem2res2("CA"));
		       stem2copy.addAtom(stem2res2("C"));	   

		       



		       t.rmsdAlignment(stem2copy.getAtomPointers(),ref2.getAtomPointers(),stem2.getAtomPointers());
		       t.rmsdAlignment(stem2copy.getAtomPointers(),ref2.getAtomPointers());
		       double rmsd2 = stem2copy.getAtomPointers().rmsd(ref2.getAtomPointers());

		       if (rmsd2 > opt.rmsdTol) continue;



		       dualLoopIndex++;

		       char dualLoopFile[100];
		       sprintf(dualLoopFile,"%s.loop.%05d.pdb",MslTools::getFileName(fileName).c_str(),dualLoopIndex);
		       dualLoopAtoms.writePdb(dualLoopFile);

		       char binderFile[100];
		       sprintf(binderFile,"%s.binder.%05d.pdb",MslTools::getFileName(fileName).c_str(),dualLoopIndex);
		       stem2.writePdb(binderFile);

		       AtomSelection stem2sel(stem2.getAtomPointers());
		       
		       double stem2rmsd = stem2_orig_ca.rmsd(stem2sel.select("name CA"));

		       fprintf(stdout, "%25s Lp1: %1s %4d %3s, %1s %4d %3s  Lp2: %1s %4d %3s, %1s %4d %3s  ; dist1 = %8.3f , dist2 = %8.3f, RMSD = %8.3f, %8.3f , %8.3f -> PDBs: %s %s\n",fileName.c_str(),
			       loop1Res1.getChainId().c_str(), loop1Res1.getResidueNumber(),loop1Res1.getResidueName().c_str(),
			       loop1Res2.getChainId().c_str(), loop1Res2.getResidueNumber(),loop1Res2.getResidueName().c_str(),
			       loop2Res1.getChainId().c_str(), loop2Res1.getResidueNumber(),loop2Res1.getResidueName().c_str(),
			       loop2Res2.getChainId().c_str(), loop2Res2.getResidueNumber(),loop2Res2.getResidueName().c_str(),dist1,dist2,rmsd,rmsd2,stem2rmsd,
			       dualLoopFile,binderFile);			       

		       
		     }


		   } // R4
		 } // R3

	       } // R2
	     } // R1
	   

	   } // FOR CHAINS
	 } // FOR lines

	 } // PRAGMA OMP

	fprintf(stdout, "Done. Time: %8.6f\n", (t.getWallTime() - startTime));
}
	


Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;
	OP.setRequired(opt.required);	
	OP.setAllowed(opt.optional);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
	OP.readArgv(theArgc, theArgv);


	if (OP.countOptions() == 0){
		cout << "Usage: dualLoopSearch " << endl;
		cout << endl;
		cout << "\n";
		cout << "pdblist PDBLIST\n";
		cout << "distanceStem1 DIST\n";
		cout << "distanceStem2 DIST\n";

		cout << "loop1min LEN\n";
		cout << "loop1max LEN\n";
		cout << "loop2min LEN\n";
		cout << "loop2max LEN\n";
		
		cout << endl;
		exit(0);
	}

	opt.pdblist = OP.getString("pdblist");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdblist specified."<<endl;
		exit(1111);
	}


	opt.distanceStem1 = OP.getDouble("distanceStem1");
	if (OP.fail()){
		cerr << "ERROR 1111 no distanceStem1 specified."<<endl;
		exit(1111);
	}

	opt.distanceStem2 = OP.getDouble("distanceStem2");
	if (OP.fail()){
		cerr << "ERROR 1111 no distanceStem2 specified."<<endl;
		exit(1111);
	}


	opt.loop1min = OP.getInt("loop1min");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop1min specified."<<endl;
		exit(1111);
	}

	opt.loop1max = OP.getInt("loop1max");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop1max specified."<<endl;
		exit(1111);
	}

	opt.loop2min = OP.getInt("loop2min");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop2min specified."<<endl;
		exit(1111);
	}

	opt.loop2max = OP.getInt("loop2max");
	if (OP.fail()){
		cerr << "ERROR 1111 no loop2max specified."<<endl;
		exit(1111);
	}

	opt.stem1pdb = OP.getString("stem1pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem1pdb specified."<<endl;
		exit(1111);
	}
	opt.stem1res1 = OP.getString("stem1res1");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem1res1 specified."<<endl;
		exit(1111);
	}
	opt.stem1res2 = OP.getString("stem1res2");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem1res2 specified."<<endl;
		exit(1111);
	}

	opt.stem2pdb = OP.getString("stem2pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem2pdb specified."<<endl;
		exit(1111);
	}
	opt.stem2res1 = OP.getString("stem2res1");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem2res1 specified."<<endl;
		exit(1111);
	}
	opt.stem2res2 = OP.getString("stem2res2");
	if (OP.fail()){
		cerr << "ERROR 1111 no stem2res2 specified."<<endl;
		exit(1111);
	}

	opt.rmsdTol = OP.getDouble("rmsdTol");
	if (OP.fail()){
	  opt.rmsdTol = 0.3;
	}
	return opt;
}
