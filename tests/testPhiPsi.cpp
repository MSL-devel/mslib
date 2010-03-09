/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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



int main(){

	cout << "Read /export/home/dwkulp/projects/PhiPsiPotential/v2/phiPsiCounts.txt"<<endl;
	PhiPsiReader ppr("/export/home/dwkulp/projects/PhiPsiPotential/v2/phiPsiCounts.txt");
	ppr.open();
	ppr.read();
	ppr.close();

	cout << "Read a string pdb 'fourHelixBundle'"<<endl;
	PDBReader pdbin(fourHelixBundle);
	pdbin.read();
	pdbin.close();

	cout << "Get Statistics for each residue"<<endl;
	PhiPsiStatistics &pps = ppr.getPhiPsiStatistics();

	System sys(pdbin.getAtoms());
	cout << "Number of residues: "<<sys.getChain("A").size()<<endl;
	for (uint i = 0; i < sys.getChain("A").size();i++){
		if (!(i > 0 && i < sys.getChain("A").size()-1)) continue;

		Residue & nm1 = sys.getChain("A").getResidue(i-1);
		Residue & n   = sys.getChain("A").getResidue(i);
		Residue & np1 = sys.getChain("A").getResidue(i+1);

		double phi = PhiPsiStatistics::getPhi(nm1,n);
		double psi = PhiPsiStatistics::getPsi(n,np1);
		int counts = pps.getCounts(nm1,n,np1);
		double prob = pps.getProbability(nm1,n,np1);
		double probAll = pps.getProbabilityAll(nm1,n,np1);
		double prop = pps.getPropensity(nm1,n,np1);

		fprintf(stdout, "%-15s %1s %3d %5d  %5.2f  %5.2f  %5.2f\n", "4HelixBundle" , "A", n.getResidueNumber(), counts, prob, probAll, prop);

		
		
	}



}
