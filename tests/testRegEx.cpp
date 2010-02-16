
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




	// Generic string matching
	string regExpression = "^CLUSTAL W (\\S+)\\s+([\\S\\s]+)$";
        string lineToMatch="CLUSTAL W 2.1 multiple sequence alignment";
        string lineNotToMatch="LUTAL W 2.1 multiple sequence alignment";

	vector<string> results;

	// Line which should match RegEx
	if (MslTools::regex(lineToMatch,regExpression,results)){

	  cout << "PASS test 1\n";

	  // Print out the matched '()' substrings
	  for (uint i = 0; i < results.size();i++){
	    cout << "Match["<<i<<"]: "<<results[i]<<endl;
	  }

	} else {
	  cerr << "ERROR did not match line properly.\n";
	}


	// Line which should NOT match RegEx
	if (MslTools::regex(lineNotToMatch,regExpression,results)){
	  cerr << "ERROR matched line when it was not suppose to\n";
	} else {
	  cout << "PASS test 2 , did not match a string that it should not match\n";
	}


}
