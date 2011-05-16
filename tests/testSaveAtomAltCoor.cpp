#include <iostream>

#include "Atom.h"

using namespace MSL;
using namespace std;

void printAtom(Atom & _atom) {
	unsigned int curr = _atom.getActiveConformation();
	for (unsigned int i=0; i<_atom.getNumberOfAltConformations(); i++) { 
		if (i == curr) {
			cout << "* ";
		} else {
			cout << "  ";
		}
		_atom.setActiveConformation(i);
		cout << _atom << endl;
	}
	_atom.setActiveConformation(curr);
	cout << endl;
}

int main(int argc, char *argv[]) {

	// set the cycle to a high number to check for memory leaks
	for (unsigned int i=0; i<1000000; i++) {
		cout << "Create an atom with one conf at 0 0 0 and save current coor as \"trial1\"" << endl;
		Atom atomA("atomA", 0.0, 0.0, 0.0);
		printAtom(atomA);
		atomA.saveCoor("trial1");

		cout << "Change the coor to 1 0 0 and save current coor as \"trial2\"" << endl;
		atomA.setCoor(1.0, 0.0, 0.0);
		atomA.saveCoor("trial2");
		printAtom(atomA);
		
		cout << "Restore \"trial1\" conf (should be 0 0 0)" << endl;
		atomA.applySavedCoor("trial1");
		printAtom(atomA);


		// Test with alt coors
		cout << "Add alternate coor: 2 0 0, 3 0 0, 4 0 0 and save all alt coors as \"trial3\"" << endl;
		atomA.addAltConformation(2.0, 0.0, 0.0);
		atomA.addAltConformation(3.0, 0.0, 0.0);
		atomA.addAltConformation(4.0, 0.0, 0.0);
		printAtom(atomA);
		atomA.saveAltCoor("trial3");

		cout << "remove all alt coor" << endl;
		atomA.removeAllAltConformations();

		cout << "print remaining coor, should be 0 0 0" << endl;
		printAtom(atomA);

		cout << "restore \"trial3\", should be *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		printAtom(atomA);

		cout << "change all coordinates, and print should be: 0 0 0, 1 1 1, 2 2 2, *3 3 3" << endl;
		for (unsigned int i=0; i < atomA.getNumberOfAltConformations(); i++) { 
			atomA.setActiveConformation(i);
			atomA.setCoor(i,i,i);
		}
		printAtom(atomA);

		cout << "save all alt coors as \"trial4\"" << endl;
		atomA.saveAltCoor("trial4");

		cout << "restore \"trial3\", should be: *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		printAtom(atomA);
		
		cout << "set the 3rd coor as the active one, should be: 0 0 0, 2 0 0, *3 0 0, 4 0 0" << endl;
		atomA.setActiveConformation(2);
		printAtom(atomA);

		cout << "restore \"trial2\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial2");
		printAtom(atomA);

		cout << "save all alt coors as \"trial5\"" << endl;
		atomA.saveAltCoor("trial5");

		cout << "add another alt conf 5 0 0, should be: 0 0 0, 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.addAltConformation(5.0, 0.0, 0.0);
		printAtom(atomA);
		cout << "save all alt coors as \"trial6\"" << endl;
		atomA.saveAltCoor("trial6");

		cout << "set the 2nd as the active conformation: should be 0 0 0, *2 0 0, 1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.setActiveConformation(1);
		printAtom(atomA);

		cout << "restore \"trial5\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial5");
		printAtom(atomA);

		cout << "Remove the 2nd alt conf, should be: 0 0 0, *1 0 0, 4 0 0" << endl;
		atomA.removeAltConformation(1);
		printAtom(atomA);
		cout << "save all alt coors as \"trial7\"" << endl;
		atomA.saveAltCoor("trial7");

		cout << "restore \"trial3\": should be *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		printAtom(atomA);

		cout << "restore \"trial7\": should be 0 0 0, *1 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial7");
		printAtom(atomA);

		cout << "remove all alt conf, should be: 1 0 0" << endl;
		atomA.removeAllAltConformations();
		printAtom(atomA);
	}

	return 0;


}

