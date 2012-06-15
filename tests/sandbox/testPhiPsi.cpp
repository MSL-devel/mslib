/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/


#include "PhiPsiReader.h"
#include "PhiPsiStatistics.h"
#include "PDBReader.h"
#include "System.h"
#include "testData.h"


using namespace std;

using namespace MSL;


#include "MslOut.h"
static MslOut MSLOUT("testPhiPsi");
int main(){

	// Mute all objects and turn this test output on
	MSLOUT.turnAllOff();
	MSLOUT.turnOn("testPhiPsi");


	string phiPsiCounts = "./tables/phiPsiCounts.txt";
	MSLOUT.stream() << "Read "<<phiPsiCounts<<endl;
	PhiPsiReader ppr(phiPsiCounts);
	ppr.open();
	ppr.read();
	ppr.close();

	MSLOUT.stream() << "Read a string pdb 'fourHelixBundle'"<<endl;
	writePdbFile();

	PDBReader pdbin;
	pdbin.open("/tmp/pdbDimer.pdb");
	pdbin.read();
	pdbin.close();


	MSLOUT.stream() << "Get Statistics for each residue (# atoms "<<pdbin.getAtomPointers().size()<<")"<<endl;
	PhiPsiStatistics &pps = ppr.getPhiPsiStatistics();

	System sys(pdbin.getAtomPointers());
	MSLOUT.stream() << "Number of residues: "<<sys.getChain("A").positionSize()<<endl;
	for (uint i = 0; i < sys.getChain("A").positionSize();i++){

		if (!(i > 0 && i < sys.getChain("A").positionSize()-1)) continue;

		Residue & nm1 = sys.getChain("A").getResidue(i-1);
		Residue & n   = sys.getChain("A").getResidue(i);
		Residue & np1 = sys.getChain("A").getResidue(i+1);
		//double phi = PhiPsiStatistics::getPhi(nm1,n);
		//double psi = PhiPsiStatistics::getPsi(n,np1);
		int counts = pps.getCounts(nm1,n,np1);
		double prob = pps.getProbability(nm1,n,np1);
		double probAll = pps.getProbabilityAll(nm1,n,np1);
		double prop = pps.getPropensity(nm1,n,np1);

		MSLOUT.fprintf(stdout, "%-15s %1s %3d %5d  %5.2f  %5.2f  %5.2f\n", "/tmp/pdbDimer.pdb" , "A", n.getResidueNumber(), counts, prob, probAll, prop);
		
	}



	// Do get random phi/psi
	PhiPsiStatistics pps2;
	pps2.addStatisitics("foo", "-10","-50",100); 
	pps2.addStatisitics("foo", "-20","-50",10);
	pps2.addStatisitics("foo", "0","-40",1);
	pps2.addStatisitics("foo", "10","-30",5);

	// Counter is just for summarizing how many of each phi psi pairs were chosen
	std::map<std::string,int> counter;
	for (uint i = 0; i < 1000;i++){

	    // Get a random pair of Phi/Psi values
	    std::pair<double,double> angles = pps2.getRandomPhiPsi("foo");


	    // Store result in counter
	    stringstream ss;
	    ss << angles.first<<":"<<angles.second;

	    counter[ss.str()]++;
	}	


	// Output the results , how many of each phi/psi pairs got picked.
	MSLOUT.stream() << "\n *** SUMMARY OF GETTING RANDOM PHI/PSI RESULTS *** \n";
	
	std::map<std::string,int> refMap;
	refMap["-5:-45"]  = 100;
	refMap["-15:-45"] = 10;
	refMap["5:-35"]   = 1;
	refMap["15:-25"]  = 5;

        MSLOUT.fprintf(stdout,"%-20s  %-8s %-8s\n","AngleBinMid Phi:Psi","REF","COUNTED");
	std::map<std::string,int>::iterator it;
	for (it = counter.begin();it != counter.end();it++){
	    MSLOUT.fprintf(stdout,"%-20s  %-8d %-4d\n",it->first.c_str(),refMap[it->first],it->second);	    
	}


}
