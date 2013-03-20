#include "Transforms.h"
#include "CharmmSystemBuilder.h"
#include "AtomSelection.h"

using namespace std;
using namespace MSL;

int main(int argc, char* argv[]) {
	// Goal is to stitch the TMD and CCD of FtsB using extracted fragments
	// Align the fragments against the CCD and TM

	if(argc != 7) {
		cerr << "Usage: connectWithFragments <TMD> <CCD> <fragFileList> <outDir> <NresStart> <CresEnd>" << endl;
		exit(0);
	}

	string outDir = string(argv[4]);

	int nStart = MslTools::toInt(string(argv[5]));
	int cEnd = MslTools::toInt(string(argv[6]));
	// lets assume we align only the two residues at each end
	
	string resNRange = MslTools::intToString(nStart) + "-" + MslTools::intToString(nStart+1);
	string resCRange = MslTools::intToString(cEnd-1) + "-" + MslTools::intToString(cEnd);
	string rigidTM = "1-" + MslTools::intToString(nStart-1);
	string rigidCCD = MslTools::intToString(cEnd+1) + "-60";

	//cout << resNRange << endl;
	//cout << resCRange << endl;
	//cout << rigidTM << endl;
	//cout << rigidCCD << endl;

	System tmd;
	if(!tmd.readPdb(string(argv[1]))) {
		cerr << "Unable to open " << argv[1] << endl;
		exit(0);
	}

	AtomSelection selTMD(tmd.getAtomPointers());
	//AtomPointerVector& tmdC = selTMD.select("tmdC,resi 20-21 and name N+CA+C+O");
	AtomPointerVector& tmdC = selTMD.select("tmdC,resi " + resNRange + " and name N+CA+C+O");

	AtomPointerVector& tmdAtoms = tmd.getAtomPointers();

	System ccd;
	if(!ccd.readPdb(string(argv[2]))) {
		cerr << "Unable to open " << argv[2] << endl;
		exit(0);
	}

	AtomSelection selCCD(ccd.getAtomPointers());
	//AtomPointerVector& ccdAN = selCCD.select("ccdAN,chain A and resi 26-27 and name N+CA+C+O");
	//AtomPointerVector& ccdBN = selCCD.select("ccdBN,chain B and resi 26-27 and name N+CA+C+O");
	AtomPointerVector& ccdAN = selCCD.select("ccdAN,chain A and resi " + resCRange + " and name N+CA+C+O");
	AtomPointerVector& ccdBN = selCCD.select("ccdBN,chain B and resi " + resCRange + " and name N+CA+C+O");

	ifstream fragListFile;
	fragListFile.open(argv[3]);
	if(!fragListFile.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		exit(0);
	}

	vector<string> fragList;
	string line;
	while(getline(fragListFile,line)) {
		line = MslTools::uncomment(line);
		line = MslTools::trim(line);

		if(line.length() > 1) {
			fragList.push_back(line);
			//cout << fragList.back() << endl;
		}
	}

	Transforms trans;


	for(int i = 0; i < fragList.size(); i++) {
		System frag;
		if(!frag.readPdb(fragList[i])) {
			cerr << "Unable to open " << fragList[i] << endl;
			exit(0);
		}

		Chain& chain1 = frag.getChain(0);
		chain1.setChainId("A");
		//chain1.renumberChain(20);
		chain1.renumberChain(nStart);

		if(!frag.duplicateChain("A","B")) {
			cerr << "Unable to duplicate chain " << fragList[i] << endl;
			exit(0);
		}

		Chain& chain2 = frag.getChain(1);

		AtomPointerVector& fragAtoms = frag.getAtomPointers();
		AtomPointerVector& chain1Atoms = chain1.getAtomPointers();
		AtomPointerVector& chain2Atoms = chain2.getAtomPointers();

		AtomSelection selFrag(fragAtoms);
		//AtomPointerVector& fragN = selFrag.select("fragN, resi 20-21 and name N+CA+C+O");
		//AtomPointerVector& fragCA = selFrag.select("fragCA, chain A and resi 26-27 and name N+CA+C+O");
		//AtomPointerVector& fragCB = selFrag.select("fragCB, chain B and resi 26-27 and name N+CA+C+O");

		AtomPointerVector& fragN = selFrag.select("fragN, resi " + resNRange + " and name N+CA+C+O");
		AtomPointerVector& fragCA = selFrag.select("fragCA, chain A and resi " + resCRange + " and name N+CA+C+O");
		AtomPointerVector& fragCB = selFrag.select("fragCB, chain B and resi " + resCRange + " and name N+CA+C+O");

		double score = 0.0;
		if(!trans.rmsdAlignment(fragCA,ccdAN,chain1Atoms)) {
			cerr << "Unable to align chainA " << fragList[i] << endl;
			exit(0);
		}

		score += trans.getLastRMSD();
		double score1 = trans.getLastRMSD();

		if(!trans.rmsdAlignment(fragCB,ccdBN,chain2Atoms)) {
			cerr << "Unable to align chainB " << fragList[i] << endl;
			exit(0);
		}
		score += trans.getLastRMSD();
		double score2 = trans.getLastRMSD();
		
		if(!trans.rmsdAlignment(tmdC,fragN,tmdAtoms)) {
			cerr << "Unable to align tmd " << fragList[i] << endl;
			exit(0);
		}
		score += trans.getLastRMSD();
		double score3 = trans.getLastRMSD();

	
		if(score3 < 5.0) {
			
			/*
			if(!frag.writePdb(outDir + "/aligned_" +fragList[i] )) {
				cerr << "unable to write output" << endl;
				exit(0);
			}

			//frag.addAtoms(selTMD.select("rigidTM, resi 1-18"));	
			//frag.addAtoms(selCCD.select("rigidCCD,resi 29-60"));	
			*/

			frag.addAtoms(selTMD.select("rigidTM, resi " + rigidTM ));	
			frag.addAtoms(selCCD.select("rigidCCD,resi " + rigidCCD));	

			if(!frag.writePdb(outDir + "/stitched_" +fragList[i] )) {
				cerr << "unable to write output" << endl;
				exit(0);
			}

			/*
			if(!tmd.writePdb(outDir + "/tmd_" +fragList[i] )) {
				cerr << "unable to write output" << endl;
				exit(0);
			}
			*/
			
		}

		cout << fragList[i] << " " << score1 << " " << score2 << " " << score3 << " " << score << endl;
	}

}



/*
int main(int argc, char* argv[]) {
	// Goal is to stitch the TMD and CCD of FtsB using extracted fragments
	// Align the fragments against the CCD so that the CCD is extended using the fragment
	// Then score each such extension against the TM

	if(argc != 4) {
		cerr << "Usage: connectWithFragments <TMD> <CCD> <fragFileList>" << endl;
		exit(0);
	}
	
	System tmd;
	if(!tmd.readPdb(string(argv[1]))) {
		cerr << "Unable to open " << argv[1] << endl;
		exit(0);
	}

	AtomSelection selTMD(tmd.getAtomPointers());
	AtomPointerVector& tmdBB = selTMD.select("tmdBB,resi 20-21 and name N+CA+C+O");

	vector<vector<double> > interHelixDist;
	// compute the interhelical distances
	for(int i  = 0; i < 8; i++) {
		interHelixDist.push_back(vector<double>());
		for(int j = 8; j < 16; j++) {
			interHelixDist.back().push_back(tmdBB[i]->distance(*tmdBB[j]));
		}
	}

	//for(vector<vector<double> >::iterator it1 = interHelixDist.begin(); it1 != interHelixDist.end(); it1++) {
	//	for(vector<double>::iterator it2 = it1->begin(); it2 != it1->end(); it2++) {
	//		cout << *it2 << " ";
	//	}
	//	cout << endl;
	//}

	System ccd;
	if(!ccd.readPdb(string(argv[2]))) {
		cerr << "Unable to open " << argv[2] << endl;
		exit(0);
	}

	AtomSelection selCCD(ccd.getAtomPointers());
	AtomPointerVector& ccdAEnd = selCCD.select("ccdAEnd,chain A and resi 25-27 and name N+CA+C+O");
	AtomPointerVector& ccdBEnd = selCCD.select("ccdBEnd,chain B and resi 25-27 and name N+CA+C+O");

	ifstream fragListFile;
	fragListFile.open(argv[3]);
	if(!fragListFile.is_open()) {
		cerr << "Unable to open " << argv[3] << endl;
		exit(0);
	}

	vector<string> fragList;
	string line;
	while(getline(fragListFile,line)) {
		line = MslTools::uncomment(line);
		line = MslTools::trim(line);

		if(line.length() > 1) {
			fragList.push_back(line);
			//cout << fragList.back() << endl;
		}
	}

	vector<double> scores;

	Transforms trans;

	for(int i = 0; i < fragList.size(); i++) {
		System frag;
		if(!frag.readPdb(fragList[i])) {
			cerr << "Unable to open " << fragList[i] << endl;
			exit(0);
		}

		Chain& chain1 = frag.getChain(0);
		AtomPointerVector& chain1Atoms = chain1.getAtomPointers();
		chain1.setChainId("A");
		chain1.renumberChain(0);

		if(!frag.duplicateChain("A","B")) {
			cerr << "Unable to duplicate chain " << fragList[i] << endl;
			exit(0);
		}
		Chain& chain2 = frag.getChain(1);
		AtomPointerVector& chain2Atoms = chain2.getAtomPointers();

		AtomSelection selFrag(frag.getAtomPointers());
		AtomPointerVector& fragA = selFrag.select("fragA, chain A and resi 5-7 and name N+CA+C+O");
		AtomPointerVector& fragB = selFrag.select("fragB, chain B and resi 5-7 and name N+CA+C+O");

		AtomPointerVector& distAtoms = selFrag.select("distAtoms, resi 0-1 and name N+CA+C+O");

		if(!trans.rmsdAlignment(fragA,ccdAEnd,chain1Atoms)) {
			cout << fragA.size() << " " << ccdAEnd.size() << endl;
			cerr << "Unable to align chainA " << fragList[i] << endl;
			exit(0);
		}
		
		if(!trans.rmsdAlignment(fragB,ccdBEnd,chain2Atoms)) {
			cout << fragB.size() << " " << ccdBEnd.size() << endl;
			cerr << "Unable to align chainB " << fragList[i] << endl;
			exit(0);
		}
		// Good - extended both ccd chains by aligning the fragment
		// now measure the distances to compute the score

		vector<vector<double> > interHelixFragDist;
		// compute the interhelical distances
		for(int r  = 0; r < 8; r++) {
			interHelixFragDist.push_back(vector<double>());
			for(int c = 8; c < 16; c++) {
				interHelixFragDist.back().push_back(distAtoms[r]->distance(*distAtoms[c]));
			}
		}

		double score = 0.0;
		//for(vector<vector<double> >::iterator it1 = interHelixFragDist.begin(); it1 != interHelixFragDist.end(); it1++) {
		for(int r = 0; r < 8; r++) {
			//for(vector<double>::iterator it2 = it1->begin(); it2 != it1->end(); it2++) {
			for(int c = 0; c < 8; c++) {
				//cout << *it2 << " ";
				//cout << interHelixFragDist[r][c] << " ";
				score += (interHelixFragDist[r][c]-interHelixDist[r][c]) * (interHelixFragDist[r][c]-interHelixDist[r][c]);
			}
			//cout << endl;
		}
		if(!frag.writePdb("aligned_" +fragList[i] )) {
			cerr << "unable to write output" << endl;
			exit(0);
		}

		cout << fragList[i] << " " << score << endl;
	}






}
*/
