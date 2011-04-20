
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
