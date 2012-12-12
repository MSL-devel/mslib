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

#include "AtomSelection.h"
#include "PDBWriter.h"
#include "MslOut.h"

using namespace MSL;
using namespace std;

// MslOut 
static MslOut MSLOUT("testAtomSelection");


int main(int argc, char *argv[]) {
  
        if (argc == 2){
	  MSLOUT.turnAllOn();
	}

	
	Atom a("A,1,ILE,CA", 0.3, 9.7, 5.3);
	Atom b("A,2,LEU,CB", 0.4, 9.7, 4.3);
	Atom c("A,3,LEU,CG", 0.5, 9.8, 5.3);
	Atom d("A,4,TYR,OH", 0.5, 9.8, 5.3);
	Atom e("A,5,LYS,NZ", 0.5, 9.8, 5.3);
	Atom f("A,6,LEU,CD1", 7.5, 9.8, 5.3);
	Atom g("A,7,THR,CA", 5.5, 9.8, 5.3);
	Atom h("A,7,THR,OG1", 0.5, 9.8, 5.3);
	Atom i("A,8,GLY,C", 0.5, 9.8, 5.3);
	Atom j("A,9,TRP,O", 0.5, 9.8, 5.3);
	Atom k("A,10,LEU,N", 0.5, 9.8, 5.3);
	Atom l("B,7,GLY,C", 0.5, 9.8, 5.3);
	Atom m("B,8,TRP,O", 0.5, 9.8, 5.3);
	Atom n("B,9,LEU,N", 0.5, 9.8, 5.3);

	AtomPointerVector av;	
	av.push_back(&a);
	av.push_back(&b);
	av.push_back(&c);
	av.push_back(&d);
	av.push_back(&e);
	av.push_back(&f);
	av.push_back(&g);
	av.push_back(&h);
	av.push_back(&i);
	av.push_back(&j);
	av.push_back(&k);
	av.push_back(&l);
	av.push_back(&m);
	
	PDBWriter pout;
	pout.open("/tmp/atomSel.pdb");
	pout.write(av);
	pout.close();
	
	// Atom Selection Object
	AtomSelection sel(av);
	sel.setDebugFlag(true);

	// By default atoms of 'all' built in, so we should be able to select it.
	AtomPointerVector allAts = sel.select("all");

	// Print out atoms that were selected
	cout << endl<<"Create selection named 'all' using the command: " << endl;
	cout << "   sel.select(\"all\")" << endl;
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("all")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("all")<<endl;
	}


	// Do a selection, name it 'foo'
	AtomPointerVector subset = sel.select("foo, name CA AND y 9.7");
	

	cout << endl;
	cout << "=====================================" << endl;
	// Print out atoms that were selected
	cout << endl<<"Create selection named 'foo' using the command: " << endl;
	cout << "   AtomPointerVector subset = sel.select(\"foo, name CA AND y 9.7\")" << endl;
	cout << "Selected " << subset.size() << "("<<sel.size("foo")<<") atoms out of "<< av.size() << endl;
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("foo")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("foo")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;

	if (subset.size() != 1){
	  cerr << "ERROR with previous selection, should have selected 1 atoms"<<endl;
	  exit(1);
	}
	// Retreive 'foo' named selection
	if (sel.selectionExists("foo")) {
		AtomPointerVector subset_again = sel.getSelection("foo");
		cout << endl<<"Retreive again the selected atoms by calling the selection with: " << endl;
		cout << "   AtomPointerVector subset_again = sel.select(\"foo\");" << endl;
		cout << "Retreived " << subset.size() << " atoms out of " <<av.size() << endl;
		for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
			//CartesianPoint coor = (*it)->getCoor();
			//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("foo")<<endl;
			cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("foo")<<endl;
		}
	}




	AtomPointerVector backbone = sel.select("bb, name CA OR name N OR name C OR name O");

	cout << endl;
	cout << endl<<"Create a backbone selection named 'bb' using the command: " << endl;
	cout << "   AtomPointerVector backbone = sel.select(\"bb, name CA OR name N OR name C OR name O\");" << endl;
	cout << "Selected " << backbone.size() << " atoms out of "<< av.size() << endl;
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("bb")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("bb")<<endl;
	}
	if (backbone.size() != 7){
	  cerr << "ERROR with previous selection, should have selected 7 atoms"<<endl;
	  exit(1);
	}
	

	cout << endl;
	cout << "=====================================" << endl;
	// Selections using distances...
	cout << "Create a CA selection named 'distSel', based on a distance from another selection 'name CB'\n";
	cout << "    AtomPointerVector sphereOfCas = sel.select(\"distSel1, name CA WITHIN 2 OF name CB\");"<<endl;

	AtomPointerVector sphereOfCas1 = sel.select("distSel1, name CA WITHIN 2 OF name CB");
	cout << "Selected " << sphereOfCas1.size() << " atoms out of "<< av.size() << endl;
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("distSel")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("distSel1")<<endl;
	}
	if (sphereOfCas1.size() != 1) {
	  cerr << "ERROR with previous selection, should have selected 1 CA atom"<<endl;
	  exit(1);
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Selections using distances...
	cout << "Create a CA selection named 'distSel', based on a distance from another selection 'name CB'\n";
	cout << "    AtomPointerVector sphereOfCas = sel.select(\"distSel2, (name CA OR resn LEU) WITHIN 4 OF name CB\");"<<endl;
	AtomPointerVector sphereOfCas2 = sel.select("distSel2, (name CA OR resn LEU) WITHIN 4 OF name CB");
	cout << "Selected " << sphereOfCas2.size() << " atoms out of "<< av.size() << endl;
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("distSel")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("distSel2")<<endl;
	}
	if (sphereOfCas2.size() != 4) {
	  cerr << "ERROR with previous selection, should have selected 4 CA atoms"<<endl;
	  exit(1);
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Selections using distances...
	cout << "Create a CA selection named 'distSel', based on a distance from another selection 'name CB'\n";
	cout << "    AtomPointerVector sphereOfCas = sel.select(\"distSel3, name CA WITHIN 4 OF (name CB or resn LEU)\");"<<endl;
	AtomPointerVector sphereOfCas3 = sel.select("distSel3, name CA WITHIN 4 OF (name CB or resn LEU)");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("distSel")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("distSel3")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Selections using distances...
	cout << "Create a CA selection named 'distSel', based on a distance from another selection 'name CB'\n";
	cout << "    AtomPointerVector sphereOfCas = sel.select(\"distSel4, (name CA OR resn LEU) WITHIN 4 OF (name CB or resn LEU)\");"<<endl;
	AtomPointerVector sphereOfCas4 = sel.select("distSel4, (name CA OR resn LEU) WITHIN 4 OF (name CB or resn LEU)");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("distSel")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("distSel4")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;

	// How many should be selected?


	// Selections using names....
	cout << "Create a selection using the name of a previous selection.."<<endl;
	cout << "    AtomPointerVector test = sel.select(\"new, bb AND name O\");"<<endl;
	AtomPointerVector test = sel.select("new, bb AND name O");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("new")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("new")<<endl;
	}


	cout << endl;
	cout << "=====================================" << endl;
	// Selections using '+' operator
	cout << "Create a selection using the '+' operator.."<<endl;
	cout << "    AtomPointerVector bb2 = sel.select(\"bb2, name CA+C+O+N\");"<<endl;
	AtomPointerVector bb2 = sel.select("bb2,name CA+C+O+N");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("bb2")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("bb2")<<endl;
	}


	cout << endl;
	cout << "=====================================" << endl;
	// Selections using NOT operator
	cout << "Create a selection using the NOT operator.."<<endl;
	cout <<"      AtomPointerVector notBB = sel.select(\"notBB, not name CA+C+O+N\")"<<endl;
	AtomPointerVector notBB = sel.select("notBB, not name CA+C+O+N");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("notBB")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("notBB")<<endl;
	}

	// Selections using range "-" operator.
	//int n = 1;
	//for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
	//	(*it)->setResidueNumber(n++);
	//}

	cout << endl;
	cout << "=====================================" << endl;
	cout << "Create a selection using the '-' operator.."<<endl;
	cout <<"      AtomPointerVector res26 = sel.select(\"res26, resi 2-6\")"<<endl;
	AtomPointerVector res26 = sel.select("res26, resi 2-6");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", ResNum "<<(*it)->getResidueNumber()<<", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("res26")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("res26")<<endl;

		// For the next test "HASCRD"
		if ((*it)->getName() == "C" || (*it)->getName() == "N" || (*it)->getName() == "O"){
			(*it)->wipeCoordinates();
		}
		
	}
	
	cout << endl;
	cout << "=====================================" << endl;
	cout << "Create a selection using the 'HASCRD' operator.."<<endl;
	cout <<"      AtomPointerVector hasCoors = sel.select(\"coor, HASCRD 0 and name CA+C+O+N+CB\");\n";
	AtomPointerVector hasCoors = sel.select("coor, HASCRD 0 and name CA+C+O+N+CB");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", ResNum "<<(*it)->getResidueNumber()<<", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("coor")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("coor")<<endl;
		
		// For the next test "HASCRD"
		if ((*it)->getName() == "CA"){
			(*it)->wipeCoordinates();
			
		}
	}


	cout << endl;
	cout << "=====================================" << endl;
	cout << "Repeat the selection using the 'HASCRD' operator after wiping the CAs coordinates..."<<endl;
	cout <<"      AtomPointerVector hasCoors = sel.select(\"coor, HASCRD 0 and name CA+C+O\");\n";
	hasCoors = sel.select("coor, HASCRD 0 and name CA+C+O");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		//CartesianPoint coor = (*it)->getCoor();
		//cout << "\tAtom name " << (*it)->getName() << ", ResNum "<<(*it)->getResidueNumber()<<", coor " << coor[0] << " " << coor[1] << " " << coor[2] << " ; flag: "<<(*it)->getSelectionFlag("coor")<<endl;
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("coor")<<endl;
		
	}
	if (hasCoors.size() != 6){
	  cerr << "ERROR with previous selection, should have selected 6 atoms"<<endl;
	  exit(1);
	}
	cout << endl;
	cout << "=====================================" << endl;



	// Complex selection 1....
	cout << "Create a complex selection"<<endl;
	cout << "    AtomPointerVector cplx1 = sel.select(\"complex1, ((resi 7 and (name CA OR chain B)) AND RESN THR\");"<<endl;
	AtomPointerVector cplx1 = sel.select("complex1, ((resi 7 and (name CA OR chain B)) and resn THR)");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex1")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Complex selection 2....
	cout << "Another complex selection"<<endl;
	cout << "    AtomPointerVector cplx2 = sel.select(\"complex2, (resn LEU and resi 3-7) OR chain B\");"<<endl;
	AtomPointerVector cplx2 = sel.select("complex2, (resn LEU and resi 3-7) OR chain B");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex2")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;

	exit(0);

	// Complex selection 3....
	cout << "A complex selection with two braket levels"<<endl;
	cout << "    AtomPointerVector cplx3 = sel.select(\"complex3, HASCRD AND (resn LEU and (resi 3-7 OR chain B)) OR NAME CD1\");"<<endl;
	AtomPointerVector cplx3 = sel.select("complex3, HASCRD AND (resn LEU and (resi 3-7 OR chain B)) OR NAME CD1");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex3")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Complex selection 4....
	cout << "A complex selection with two braket levels"<<endl;
	cout << "    AtomPointerVector cplx4 = sel.select(\"complex4, (HASCRD OR resn GLY) and NOT (resi 3-7 OR (chain B AND NAME C))\");"<<endl;
	AtomPointerVector cplx4 = sel.select("complex4, (HASCRD OR resn GLY) and NOT (resi 3-7 OR (chain B AND NAME C))");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex4")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Check correct AND/OR/XOR priority, it should follow this order AND->XOR->OR....
	cout << "A long selection without brakets, the order of priority is AND -> XOR -> OR"<<endl;
	cout << "    AtomPointerVector cplx5 = sel.select(\"complex5, HASCRD OR resn GLY and resi 3-8 AND NOT NAME CB XOR chain B AND NAME C\");"<<endl;
	AtomPointerVector cplx5 = sel.select("complex5, HASCRD OR resn GLY and resi 3-8 AND NOT NAME CB XOR chain B AND NAME C");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex5")<<endl;
	}


	cout << endl;
	cout << "=====================================" << endl;
	// Check correct AND/OR/XOR priority, it should follow this order AND->XOR->OR)
	cout << "The same long selection with braket levels to enforce the AND -> XOR -> OR priority, the results should be the same of above"<<endl;
	cout << "    AtomPointerVector cplx6 = sel.select(\"complex6, HASCRD OR (((resn GLY and resi 3-8) AND NOT NAME CB) XOR (chain B AND NAME C))\");"<<endl;
	AtomPointerVector cplx6 = sel.select("complex6, HASCRD OR (((resn GLY and resi 3-8) AND NOT NAME CB) XOR (chain B AND NAME C))");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex6")<<endl;
	}

	cout << endl;
	cout << "=====================================" << endl;
	// Check correct AND/OR/XOR priority, it should follow this order AND->XOR->OR....
	cout << "Another long selection without brakets to test priority"<<endl;
	cout << "    AtomPointerVector cplx7 = sel.select(\"complex7, HASCRD AND resn GLY XOR resi 3-8 XOR NOT NAME CB OR chain B OR NAME C\");"<<endl;
	AtomPointerVector cplx7 = sel.select("complex7, HASCRD AND resn GLY XOR resi 3-8 XOR NOT NAME CB OR chain B OR NAME C");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex7")<<endl;
	}


	cout << endl;
	cout << "=====================================" << endl;
	// Check correct AND/OR/XOR priority, it should follow this order AND->XOR->OR)
	cout << "The same long selection with braket, the results should be the same of above"<<endl;
	cout << "    AtomPointerVector cplx8 = sel.select(\"complex8, ((((HASCRD AND resn GLY) XOR resi 3-8) XOR NOT NAME CB) OR chain B) OR NAME C\");"<<endl;
	AtomPointerVector cplx8 = sel.select("complex8, ((((HASCRD AND resn GLY) XOR resi 3-8) XOR NOT NAME CB) OR chain B) OR NAME C");
	for (AtomPointerVector::iterator it = av.begin();it != av.end(); it++){
		cout << "\t" << *(*it) << " ; flag: "<<(*it)->getSelectionFlag("complex8")<<endl;
	}

}


