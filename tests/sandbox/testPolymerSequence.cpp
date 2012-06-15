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

#include <iostream>
#include <string>
#include "AtomPointerVector.h"
#include "CharmmSystemBuilder.h"
#include "System.h"
#include "PDBReader.h"
#include "PolymerSequence.h"

using namespace std;

using namespace MSL;


int main(int argc, char* argv[]) {

	PolymerSequence seq;

	string input = "A: ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " A: ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "      A: ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = ": ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " : ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "     : ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "    ALA ALA ALA";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "A: ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " A: ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "      A: ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = ": ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " : ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "     : ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = " ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "    ALA ALA ALA B: LEU LEU LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "\
A: {7}ALA ALA ALA ALA\
LEU LEU [LEU ILE VAL] ALA ALA [LEU ILE ]\
\
B:ILE LEU ALA VAL C: ALA VAL LEU TYR \
D: {12}ILE VAL GLY LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "\
A: {7}ALA ALA ALA ALA\\\n\
LEU LEU [LEU ILE VAL] ALA ALA [LEU ILE ]\n\
\n\
B:ILE LEU ALA VAL C: ALA VAL LEU TYR \n\
D: {12}ILE VAL GLY LEU\n";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	input = "A:VAL GLY LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	
	seq.addPositionIdentity(0,1,"GLU");
	cout << "Adding GLU to chain A, 1" << endl;
	cout << seq << endl;
	cout << "=====" << endl;
	input = "{7}ALA ALA ALA ALA LEU LEU [LEU ILE VAL] ALA ALA [LEU ILE ] B:ILE LEU ALA VAL C: ALA VAL LEU TYR D: {12}ILE VAL GLY LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;
	input = "{7}ALA ALA ALA ALA LEU LEU [LEU ILE VAL] ALA ALA [LEU ILE ] :ILE LEU ALA VAL : ALA VAL LEU TYR : {12}ILE VAL GLY LEU";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;
	input = "\
A: {7}ALA ALA ALA ALA\\\n\
LEU LEU [LEU ILE VAL] ALA ALA [LEU ILE ]\n\
\n\
B:ILE LEU ALA VAL : ALA VAL LEU TYR \n\
D: {12}ILE VAL GLY LEU\n";
	cout << input << endl;
	seq.setSequence(input);
	cout << seq << endl;
	cout << "=====" << endl;	

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
	p2.setSequence("B: {34}[ALA] [ILE] [TRP] {37A}[LYS] {37B}[ARG] [ALA] {39C}[VAL ILE] [GLY] {99}[SER]");
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


	cout << "******* Test the setSequence(AtomPointerVector _atoms) Method *********" << endl;


	AtomPointerVector atoms;

	Atom* a1 = new Atom("CA");
	a1->setResidueName("LEU");
	a1->setResidueNumber(1);
	a1->setResidueIcode("");
	a1->setChainId("A");
	atoms.push_back(a1);

	Atom* a2 = new Atom("CA");
	a2->setResidueName("ILE");
	a2->setResidueNumber(2);
	a2->setResidueIcode("");
	a2->setChainId("A");
	atoms.push_back(a2);

	Atom* a3 = new Atom("CA");
	a3->setResidueName("TRP");
	a3->setResidueNumber(2);
	a3->setResidueIcode("");
	a3->setChainId("A");
	atoms.push_back(a3);

	Atom* a4 = new Atom("CA");
	a4->setResidueName("VAL");
	a4->setResidueNumber(3);
	a4->setChainId("A");
	atoms.push_back(a4);

	Atom* a5 = new Atom("CA");
	a5->setResidueName("PRO");
	a5->setResidueNumber(3);
	a5->setResidueIcode("A");
	a5->setChainId("A");
	atoms.push_back(a5);

	Atom* a6 = new Atom("CA");
	a6->setResidueName("GLN");
	a6->setResidueNumber(4);
	a6->setChainId("A");
	atoms.push_back(a6);

	Atom* a7 = new Atom("CA");
	a7->setResidueName("HIS");
	a7->setResidueNumber(5);
	a7->setChainId("A");
	atoms.push_back(a7);

	Atom* a8 = new Atom("CA");
	a8->setResidueName("GLY");
	a8->setResidueNumber(7);
	a8->setResidueIcode("");
	a8->setChainId("A");
	atoms.push_back(a8);

	Atom* a9 = new Atom("CA");
	a9->setResidueName("GLN");
	a9->setResidueNumber(8);
	a9->setResidueIcode("");
	a9->setChainId("A");
	atoms.push_back(a9);

	Atom* a10 = new Atom("CA");
	a10->setResidueName("TYR");
	a10->setResidueNumber(1);
	a10->setResidueIcode("");
	a10->setChainId("B");
	atoms.push_back(a10);

	Atom* a11 = new Atom("CA");
	a11->setResidueName("ALA");
	a11->setResidueNumber(2);
	a11->setResidueIcode("");
	a11->setChainId("B");
	atoms.push_back(a11);

	Atom* a12 = new Atom("CA");
	a12->setResidueName("SER");
	a12->setResidueNumber(3);
	a12->setResidueIcode("");
	a12->setChainId("B");
	atoms.push_back(a12);

	Atom* a13 = new Atom("CA");
	a13->setResidueName("LYS");
	a13->setResidueNumber(4);
	a13->setResidueIcode("");
	a13->setChainId("B");
	atoms.push_back(a13);

	cout << atoms << endl;

	PolymerSequence p3(atoms);
	cout << p3 << endl;


	if(argc == 2) {

		cout << "*********** Test the setSequence(AtomPointerVector _atoms) Method with a pdb file **********" << endl;
		string pdbPath = string(argv[1]); 


		PDBReader pRead;
		if(!pRead.open(pdbPath)) {
			cerr << "Unable to open " << pdbPath << endl;
			exit(0);
		}

		if(!pRead.read()) {
			cerr << "Unable to read " << pdbPath << endl;
			exit(0);
		}
		
		PolymerSequence p3;
		p3.setSequence(pRead.getAtomPointers());

//		csb.buildSystem(sys,p3);
		cout << p3 << endl;
	}

	return 0;
}
