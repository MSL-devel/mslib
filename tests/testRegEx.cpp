
#include "testData.h"
#include "RegEx.h"
#include "System.h"


int main(){



	writePdbFile();
	System sys;
	sys.readPdb("/tmp/xtalLattice.pdb");
	             
	cout << "\n****\nSearch chain A of /tmp/xtalLattice for A..K\n"<<endl;
	if (sys.exists("A")){
		RegEx re;
		vector<pair<int,int> > matches = re.getResidueRanges(sys("A"),"A..K");

		cout <<endl;
		for (uint i = 0; i < matches.size();i++){

			for (uint j = matches[i].first; j <= matches[i].second;j++){
				Residue &r = sys("A").getResidueByIndex(j);

				cout << "R: "<<r.toString()<<endl;
			}
			cout << " -- "<<endl;
		}

	}



}
