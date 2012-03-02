#include <iostream>

#include "Atom.h"

using namespace MSL;
using namespace std;

void printAtom(Atom & _atom) {
	/*
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
	*/
	_atom.setToStringFormat(2);
	cout << _atom << endl;
	cout << endl;
}

int main(int argc, char *argv[]) {

	// set the cycle to a high number (1000000) to check for memory leaks
	for (unsigned int i=0; i<1; i++) {
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

		cout << "Hide alt coor with relative index 1: should be 0 0 0, [1 0 0], *4 0 0" << endl;
		atomA.hideAltCoorRelIndex(1);
		printAtom(atomA);
		cout << "save all alt coors as \"trial8\"" << endl;
		atomA.saveAltCoor("trial8");

		cout << "restore \"trial3\": should be *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		printAtom(atomA);

		cout << "Hide alt coor with absolute index 1: should be *0 0 0, [2 0 0 ], 3 0 0, 4 0 0" << endl;
		atomA.hideAltCoorAbsIndex(1);
		printAtom(atomA);

		cout << "Hide alt coor with absolute index 3: should be *0 0 0, [2 0 0 ], 3 0 0, [4 0 0]" << endl;
		atomA.hideAltCoorAbsIndex(3);
		printAtom(atomA);
		cout << "save all alt coors as \"trial9\"" << endl;
		atomA.saveAltCoor("trial9");

		cout << "restore \"trial6\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial6");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide alt coor with relative index 0: should be [0 0 0], 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.hideAltCoorRelIndex(0);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide alt coor with relative index 3: should be [0 0 0], 2 0 0, *1 0 0, 4 0 0, [5 0 0]" << endl;
		atomA.hideAltCoorRelIndex(3);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the relative index 0: should be [0 0 0], *2 0 0, [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoorButOneRelIndex(0);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the absolute index 3: should be [0 0 0], [2 0 0], [1 0 0], *4 0 0, [5 0 0]" << endl;
		atomA.hideAllAltCoorButOneAbsIndex(3);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Unhide alt coor with absolute index 1: should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, [5 0 0]" << endl;
		atomA.unhideAltCoorAbsIndex(1);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Unhide alt coor with absolute index 4: should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.unhideAltCoorAbsIndex(4);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << "save all alt coors as \"trial10\"" << endl;
		atomA.saveAltCoor("trial10");

		cout << "Hide all alt coor except the first 1: should be *0 0 0, [2 0 0], [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoodButFirstNAbsIndex(1);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the first 2: should be *0 0 0, 2 0 0, [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoodButFirstNAbsIndex(2);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the first 3: should be *0 0 0, 2 0 0, 1 0 0, [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoodButFirstNAbsIndex(3);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the first 4: should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, [5 0 0]" << endl;
		atomA.hideAllAltCoodButFirstNAbsIndex(4);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Hide all alt coor except the first 5 (which happens to be all): should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, 5 0 0" << endl;
		atomA.hideAllAltCoodButFirstNAbsIndex(5);
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "Unhide all alt coor: should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, 5 0 0" << endl;
		atomA.unhideAllAltCoor();
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "remove all alt conf, should be: 1 0 0" << endl;
		atomA.removeAllAltConformations();
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

		cout << "clear all saved alt conf, should still be: 1 0 0" << endl;
		atomA.clearSavedCoor();
		printAtom(atomA);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;

	}

	return 0;


}

