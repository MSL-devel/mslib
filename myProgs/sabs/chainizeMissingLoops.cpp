#include "System.h"
#include "AtomBondBuilder.h"
#include "PolymerSequence.h"

using namespace std;
using namespace MSL;


// Read a pdb file identify missing loops, the ends opf the missing loops will be assigned to different chains.
// Make sure the inp.pdb file has charmm names and PDB format.
// Identify missing loops using the atom bond builder and assign different chain Ids to each broken segment.
// patch the termini of these new chains with uncharged patches.
// Also convert HIS to correct HSD/HSE/HSP
// if HD1 does not exist --> HSE
// if HE2 does not exist --> HSD
// else HSP

int main(int argc, char* argv[]) {
	if(argc != 4) {
		cerr << "chainizeMissingLoops <inp.pdb> <out.pdb> <out.seq>" << endl;
		exit(0);
	}

	System sysIn, sysOut;
	if(!sysIn.readPdb(string(argv[1]))) {
		cerr  << "Unable to read" << argv[1] << endl;
		exit(0);
	}

	AtomBondBuilder abb;
	abb.buildConnections(sysIn.getAtomPointers());

	int chainNum = 0;

	vector<Chain*> chains = sysIn.getChains();
	for(int i = 0; i < chains.size(); i++) {
		vector<Position*> pos = chains[i]->getPositions();
		// check all positions except first 
		for(int j = 0; j < pos.size(); j++) {
			string resName = pos[j]->getResidueName();
			if(resName == "HIS") {
				if(!pos[j]->atomExists("HD1")) {
					pos[j]->setResidueName("HSE");
				} else if (!pos[j]->atomExists("HE2")) {
					pos[j]->setResidueName("HSD");
				} else {
					pos[j]->setResidueName("HSP");
				}
			}
			if(j == 0 ) {
				continue;
			}
			if(pos[j]->atomExists("N")) {
				Atom & N = pos[j]->getLastFoundAtom();
				if(pos[j-1]->atomExists("C")) {
					Atom& prevC = pos[j-1]->getLastFoundAtom();
					if(N.isBoundTo(&prevC)) {
						if(pos[j]->getChainId() != prevC.getChainId()) {
							cout << argv[1] << " " << pos[j]->getChainId() << " " << pos[j]->getResidueNumber() << pos[j]->getResidueIcode() <<  " " << prevC.getChainId() << endl;
						}
						pos[j]->setParentChain(NULL);
						pos[j]->setChainId(prevC.getChainId());
					} else {
						cout << argv[1] << " " << pos[j]->getChainId() << " " << pos[j]->getResidueNumber() << pos[j]->getResidueIcode() <<  " " << chainNum << endl;
						pos[j]->setParentChain(NULL);
						pos[j]->setChainId(MslTools::intToString(chainNum));
						pos[j-1]->setResidueName(pos[j-1]->getResidueName() + "-CT2");
						if(pos[j]->getResidueName() == "PRO") {
							pos[j]->setResidueName(pos[j]->getResidueName() + "-ACP");
						} else {
							pos[j]->setResidueName(pos[j]->getResidueName() + "-ACE");
						}
						chainNum++;
					}
				} else {
					cerr << argv[1] << " " << pos[j-1]->getPositionId() <<  " Missing C." << endl;
				}
			} else {
				cerr << argv[1] << " " << pos[j]->getPositionId() <<  " Missing N." << endl;
			}
		}
	}
	/*
	sysOut.addAtoms(sysIn.getAtomPointers());
	if(!sysOut.writePdb(string(argv[2]))) {
		cerr << "Unable to write " << argv[2] << endl;
		exit(0);
	}

	PolymerSequence seq(sysOut);
	ofstream out(argv[3]);
	if(out.is_open()) {
		out << seq << endl;
	} else {
		cerr << "Unable to open " << argv[3] << endl;
	}
	*/
}
