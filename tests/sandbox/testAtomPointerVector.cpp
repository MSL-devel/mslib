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

using namespace MSL;
using namespace std;


int main(){


	
	Atom a("C1", 0.3, 9.7, 5.3);
	Atom b("C2", 0.4, 9.7, 4.3);
	Atom c("C3", 0.5, 9.8, 5.3);
	Atom d("C4", 0.5, 9.8, 5.3);
	Atom e("C5", 0.5, 9.8, 5.3);
	Atom f("C6", 0.5, 9.8, 5.3);
	Atom g("C7", 0.5, 9.8, 5.3);
	Atom h("C8", 0.5, 9.8, 5.3);
	Atom i("C9", 0.5, 9.8, 5.3);
	Atom j("C10", 0.5, 9.8, 5.3);

	AtomPointerVector full;	
	full.push_back(&a);
	full.push_back(&b);
	full.push_back(&c);
	full.push_back(&d);
	full.push_back(&e);
	full.push_back(&f);
	full.push_back(&g);
	full.push_back(&h);
	full.push_back(&i);
	full.push_back(&j);
	
	AtomPointerVector partA;	
	partA.push_back(&a);
	partA.push_back(&b);
	partA.push_back(&c);
	partA.push_back(&d);
	partA.push_back(&e);
	AtomPointerVector partB;	
	partB.push_back(&f);
	partB.push_back(&g);
	partB.push_back(&h);
	partB.push_back(&i);
	partB.push_back(&j);
	
	AtomPointerVector overlapA;	
	overlapA.push_back(&a);
	overlapA.push_back(&b);
	overlapA.push_back(&c);
	overlapA.push_back(&d);
	overlapA.push_back(&e);
	overlapA.push_back(&f);
	overlapA.push_back(&g);
	AtomPointerVector overlapB;	
	overlapB.push_back(&c);
	overlapB.push_back(&d);
	overlapB.push_back(&e);
	overlapB.push_back(&f);
	overlapB.push_back(&g);
	overlapB.push_back(&h);
	overlapB.push_back(&i);
	overlapB.push_back(&j);
	
	cout << "AtomPointerVector with all atoms:" << endl;
	cout << "  full (" << full.size() << "):" << endl;
	for (AtomPointerVector::iterator k=full.begin(); k<full.end(); k++) {
		cout << **k << endl;
	}
	cout << "====" << endl;
	cout << "List and combine partA and partB: " << endl;
	cout << "  partA (" << partA.size() << "):" << endl;
	for (AtomPointerVector::iterator k=partA.begin(); k<partA.end(); k++) {
		cout << **k << endl;
	}
	cout << "  partB (" << partB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=partB.begin(); k<partB.end(); k++) {
		cout << **k << endl;
	}
	AtomPointerVector partAB = partA;
	partAB.insert(partAB.end(), partB.begin(), partB.end());
	cout << "Combine using stl vector insert operator \"partAB.insert(partAB.end(), partB.begin(), partB.end())\":" << endl;
	cout << "  partAB (" << partAB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=partAB.begin(); k<partAB.end(); k++) {
		cout << **k << endl;
	}
	partAB.clear();
	partAB = partA + partB;
	cout << "Combine using sum operator (partA + partB):" << endl;
	cout << "  partAB (" << partAB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=partAB.begin(); k<partAB.end(); k++) {
		cout << **k << endl;
	}

	cout << "====" << endl;
	cout << "List and combine overlapA and overlapB: " << endl;
	cout << "  overlapA (" << overlapA.size() << "):" << endl;
	for (AtomPointerVector::iterator k=overlapA.begin(); k<overlapA.end(); k++) {
		cout << **k << endl;
	}
	cout << "  overlapB (" << overlapB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=overlapB.begin(); k<overlapB.end(); k++) {
		cout << **k << endl;
	}
	AtomPointerVector overlapAB = overlapA;
	overlapAB.insert(overlapAB.end(), overlapB.begin(), overlapB.end());
	cout << "Combine using stl vector insert operator \"overlapAB.insert(overlapAB.end(), overlapB.begin(), overlapB.end())\":" << endl;
	cout << "  overlapAB (" << overlapAB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=overlapAB.begin(); k<overlapAB.end(); k++) {
		cout << **k << endl;
	}
	overlapAB.clear();
	overlapAB = overlapA + overlapB;
	cout << "Combine using sum operator (overlapA + overlapB):" << endl;
	cout << "  overlapAB (" << overlapAB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=overlapAB.begin(); k<overlapAB.end(); k++) {
		cout << **k << endl;
	}

	cout << "====" << endl;
	cout << "Remove partB from full using the - operator (full -= partB): " << endl;
	cout << "  full (" << full.size() << "):" << endl;
	for (AtomPointerVector::iterator k=full.begin(); k<full.end(); k++) {
		cout << **k << endl;
	}
	cout << "  partB (" << partB.size() << "):" << endl;
	for (AtomPointerVector::iterator k=partB.begin(); k<partB.end(); k++) {
		cout << **k << endl;
	}
	AtomPointerVector subtract = full - partB;
	cout << "  subtract (" << subtract.size() << "):" << endl;
	for (AtomPointerVector::iterator k=subtract.begin(); k<subtract.end(); k++) {
		cout << **k << endl;
	}
}
