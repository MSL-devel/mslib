#include"System.h"
using namespace std;
using namespace MSL;


bool uniq(map<string,string>& seqMap, string sequence) {
	for(map<string,string>::iterator it = seqMap.begin(); it != seqMap.end(); it++) {
		string orig = it->first;
		// depending on which one is longer align sequence with orig or orig with seq and find the longest match
		if(orig.length() > sequence.length()) {
			if(string::npos != orig.find(sequence)) {
				return false;
			}
		} else {
			if(string::npos != sequence.find(orig)) {
				return false;
			}
		}

	}
	return true;
}

int main(int argc, char* argv[]) {
	if(argc != 2) {
		cerr << "Usage: getUniqueChains <pdbFile>" << endl;
		exit(0);
	}

	System sys;
	if(!sys.readPdb(string(argv[1]))) {
		cerr << "Unable to read " << argv[1] << endl;
		exit(0);
	}

	vector<Chain*>& chains = sys.getChains();
	map<string,string> seq; // map from seq to chainId
	vector<string> uniqIds;

	for(int i = 0; i < chains.size(); i++) {
		string chainId = chains[i]->getChainId();
		string sequence = "";
		for(int j = 0; j < chains[i]->positionSize(); j++) {
			string threeLetterName = (chains[i]->getPosition(j)).getResidueName();
			string oneLetterName = MslTools::getOneLetterCode(threeLetterName);
			if(threeLetterName == "TIP") {
				continue;
			}
			sequence += oneLetterName;
		}
		cout << argv[1] << " " << chainId << " " << sequence << endl;
		if(uniq(seq,sequence)) {
			seq[sequence] = chainId;
			uniqIds.push_back(chainId);
		}
	}

	cout << "Unique " << argv[1];
	for(int i = 0; i < uniqIds.size(); i++) {
		cout << " " << uniqIds[i];
	}
	cout << endl;

}

