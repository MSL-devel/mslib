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

#include <tr1/memory>
#include <iostream>
#include "Atom.h" 
using namespace std;
using namespace MSL;

typedef std::tr1::shared_ptr<Atom> AtomSP;

int main(){
	Atom *a = new Atom("A,37,ILE,CA");
	AtomSP sp1(a);

	// Easily access the address and the data of Atom 'a'
	cout << "sp1 address: " << sp1.get() << " toString(): "<<sp1->toString()<<endl;

	{
		AtomSP sp2(sp1);
		cout << "sp2 address: " << sp2.get() << " toString(): "<<sp2->toString()<<endl;
		cout << "Number of smart pointers to a, using sp1: " << sp1.use_count() << " or using sp2: "<<sp2.use_count()<<endl;

		// Change some value in the Atom 'a'
		sp2->setResidueName("PIZZA");
	}

	// Did the value change?
	cout << "a: " << a<< " "<<a->toString()<<endl;

	// We now have lost 'a' ... this would be a memory leak...
	a = NULL;
	cout << "a2: " << a<< " "<<endl;


	cout << "sp1 address: "<<sp1.get()<<" toString(): " <<sp1->toString()<<endl;
	cout << "Number of smart pointers: " << sp1.use_count() << endl;

	sp1.reset();
	cout << "Number of smart pointers at end: " << sp1.use_count() << endl;


	// Now create an AtomPointerVector-like construct
	cout << endl<<endl<< " *** TESTING AtomPointerVector-like construct with smart pointers ***"<<endl;

	vector<AtomSP> atomPointerVector;
	AtomSP aSmartPointer(new Atom("A,37,ILE,CA"));

	atomPointerVector.push_back(aSmartPointer);
	
	aSmartPointer.reset(new Atom("A,38,LEU,CA"));
	atomPointerVector.push_back(aSmartPointer);


	aSmartPointer.reset(new Atom("A,39,GLY,CA"));
	atomPointerVector.push_back(aSmartPointer);

	aSmartPointer.reset(new Atom("A,40,ASP,CA"));
	atomPointerVector.push_back(aSmartPointer);

	cout << "Use count of aSmartPointer: "<<aSmartPointer.use_count()<<endl;
	for (uint i =0; i < atomPointerVector.size();i++){
		cout << "Use count: "<< atomPointerVector[i].use_count()<<" "<<atomPointerVector[i].get()<<" toString(): "<<atomPointerVector[i]->toString()<<endl;
	}

	cout << "NOW I WILL CLEAR THE atomPointerVector, which will call the destructors for all the atoms except the last one!"<<endl;
	atomPointerVector.clear();

	cout << "aSmartPointer: "<<aSmartPointer.use_count()<<" "<<aSmartPointer.get()<<" "<<aSmartPointer->toString()<<endl;
	return 0;
}
