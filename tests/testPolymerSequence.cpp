/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#include <iostream>
#include <string>

#include "PolymerSequence.h"

using namespace std;

int main() {

	PolymerSequence seq;

	string input = "\
A: {7}ALA ALA ALA ALA\\\n\
LEU LEU {35}[ LEU ILE VAL] ALA ALA [ LEU ILE ]\n\
\n\
B:ILE LEU ALA VAL\n\
C: {12}  ILE VAL GLY LEU\n";

	seq.setSequence(input);
	cout << seq << endl;
	


	PolymerSequence s1;
	s1.setSequence("A 73: GLY ALA VAL ILE SER THR");

	string refSeq = "VGSTALDKJASDFAAASDFASDF";
	//s1.setReferenceSequence(refSeq,"testRefSeq", 522,526);
	//s1.setReferenceSequence(refSeq,"testRefSeq2", 522, 523, 75);
	s1.setReferenceSequence(refSeq,"testRefSeq2", 522, 526, 75);
	cout << "Aligned PolymerSequence with a reference sequence.\n";
	cout << s1.getReferenceHeader();
	cout << s1<<endl;
	


	PolymerSequence p2;
	p2.setSequence("B: {34}ALA ILE TRP {37A}LYS {37B}ARG ALA [{39C}VAL ILE] GLY {99}SER");
	cout << p2<<endl;

	for (uint i = 0; i < p2.size();i++){
		for (uint j = 0; j < p2.chainSize(i);j++){
			string resNum = p2.getResidueNumber(i,j);
			int resn;
			string icode;
			MslTools::splitIntAndString(resNum,resn,icode);
			cout << "CHAIN "<<p2.getChainId(i)<<" RESI "<<resn<<" ICODE "<<icode<<endl;
		}
	}

	// Just testing the 'smart round' feature
	cout << "Smart round: "<<MslTools::smartRound(526,10)<<endl;
	return 0;
}
