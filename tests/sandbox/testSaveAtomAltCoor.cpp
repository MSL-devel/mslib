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

#include "Atom.h"

using namespace MSL;
using namespace std;

int main(int argc, char *argv[]) {

	// set the cycle to a high number (1000000) to check for memory leaks
	for (unsigned int i=0; i<1; i++) {
		cout << "Create an atom with one conf at 0 0 0 and save current coor as \"trial1\"" << endl;
		Atom atomA("atomA", 0.0, 0.0, 0.0);
		atomA.setToStringFormat(2); // print all alternative coordinates including those hidden
		cout << atomA << endl;
		atomA.saveCoor("trial1");

		cout << "Change the coor to 1 0 0 and save current coor as \"trial2\"" << endl;
		atomA.setCoor(1.0, 0.0, 0.0);
		atomA.saveCoor("trial2");
		cout << atomA << endl;
		
		cout << "Restore \"trial1\" conf (should be 0 0 0)" << endl;
		atomA.applySavedCoor("trial1");
		cout << atomA << endl;


		// Test with alt coors
		cout << "Add alternate coor: 2 0 0, 3 0 0, 4 0 0 and save all alt coors as \"trial3\"" << endl;
		atomA.addAltConformation(2.0, 0.0, 0.0);
		atomA.addAltConformation(3.0, 0.0, 0.0);
		atomA.addAltConformation(4.0, 0.0, 0.0);
		cout << atomA << endl;
		atomA.saveAltCoor("trial3");

		cout << "remove all alt coor" << endl;
		atomA.removeAllAltConformations();

		cout << "print remaining coor, should be 0 0 0" << endl;
		cout << atomA << endl;

		cout << "restore \"trial3\", should be *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		cout << atomA << endl;

		cout << "change all coordinates, and print should be: 0 0 0, 1 1 1, 2 2 2, *3 3 3" << endl;
		for (unsigned int i=0; i < atomA.getNumberOfAltConformations(); i++) { 
			atomA.setActiveConformation(i);
			atomA.setCoor(i,i,i);
		}
		cout << atomA << endl;

		cout << "save all alt coors as \"trial4\"" << endl;
		atomA.saveAltCoor("trial4");

		cout << "restore \"trial3\", should be: *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		cout << atomA << endl;
		
		cout << "set the 3rd coor as the active one, should be: 0 0 0, 2 0 0, *3 0 0, 4 0 0" << endl;
		atomA.setActiveConformation(2);
		cout << atomA << endl;

		cout << "restore \"trial2\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial2");
		cout << atomA << endl;

		cout << "save all alt coors as \"trial5\"" << endl;
		atomA.saveAltCoor("trial5");

		cout << "add another alt conf 5 0 0, should be: 0 0 0, 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.addAltConformation(5.0, 0.0, 0.0);
		cout << atomA << endl;
		cout << "save all alt coors as \"trial6\"" << endl;
		atomA.saveAltCoor("trial6");

		cout << "set the 2nd as the active conformation: should be 0 0 0, *2 0 0, 1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.setActiveConformation(1);
		cout << atomA << endl;

		cout << "restore \"trial5\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial5");
		cout << atomA << endl;

		cout << "Remove the 2nd alt conf, should be: 0 0 0, *1 0 0, 4 0 0" << endl;
		atomA.removeAltConformation(1);
		cout << atomA << endl;
		cout << "save all alt coors as \"trial7\"" << endl;
		atomA.saveAltCoor("trial7");

		cout << "Hide alt coor with relative index 1: should be 0 0 0, [1 0 0], *4 0 0" << endl;
		atomA.hideAltCoorRelIndex(1);
		cout << atomA << endl;
		cout << "save all alt coors as \"trial8\"" << endl;
		atomA.saveAltCoor("trial8");

		cout << "restore \"trial3\": should be *0 0 0, 2 0 0, 3 0 0, 4 0 0" << endl;
		atomA.applySavedCoor("trial3");
		cout << atomA << endl;

		cout << "Hide alt coor with absolute index 1: should be *0 0 0, [2 0 0 ], 3 0 0, 4 0 0" << endl;
		atomA.hideAltCoorAbsIndex(1);
		cout << atomA << endl;

		cout << "Hide alt coor with absolute index 3: should be *0 0 0, [2 0 0 ], 3 0 0, [4 0 0]" << endl;
		atomA.hideAltCoorAbsIndex(3);
		cout << atomA << endl;
		cout << "save all alt coors as \"trial9\"" << endl;
		atomA.saveAltCoor("trial9");

		cout << "restore \"trial6\": should be 0 0 0, 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial6");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide alt coor with relative index 0: should be [0 0 0], 2 0 0, *1 0 0, 4 0 0, 5 0 0" << endl;
		atomA.hideAltCoorRelIndex(0);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide alt coor with relative index 3: should be [0 0 0], 2 0 0, *1 0 0, 4 0 0, [5 0 0]" << endl;
		atomA.hideAltCoorRelIndex(3);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the relative index 0: should be [0 0 0], *2 0 0, [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoorsButOneRelIndex(0);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the absolute index 3: should be [0 0 0], [2 0 0], [1 0 0], *4 0 0, [5 0 0]" << endl;
		atomA.hideAllAltCoorsButOneAbsIndex(3);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Unhide alt coor with absolute index 1: should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, [5 0 0]" << endl;
		atomA.unhideAltCoorAbsIndex(1);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Unhide alt coor with absolute index 4: should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.unhideAltCoorAbsIndex(4);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;
		cout << "save all alt coors as \"trial10\"" << endl;
		atomA.saveAltCoor("trial10");

		cout << "Hide all alt coor except the first 1: should be *0 0 0, [2 0 0], [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoorsButFirstN(1);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the first 2: should be *0 0 0, 2 0 0, [1 0 0], [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoorsButFirstN(2);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the first 3: should be *0 0 0, 2 0 0, 1 0 0, [4 0 0], [5 0 0]" << endl;
		atomA.hideAllAltCoorsButFirstN(3);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the first 4: should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, [5 0 0]" << endl;
		atomA.hideAllAltCoorsButFirstN(4);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Hide all alt coor except the first 5 (which happens to be all): should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, 5 0 0" << endl;
		atomA.hideAllAltCoorsButFirstN(5);
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "restore \"trial10\": should be [0 0 0], 2 0 0, [1 0 0], *4 0 0, 5 0 0" << endl;
		atomA.applySavedCoor("trial10");
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "Unhide all alt coor: should be 0 0 0, 2 0 0, 1 0 0, *4 0 0, 5 0 0" << endl;
		atomA.unhideAllAltCoors();
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "remove all alt conf, should be: 4 0 0" << endl;
		atomA.removeAllAltConformations();
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

		cout << "clear all saved alt conf, should still be: 4 0 0" << endl;
		atomA.clearSavedCoor();
		cout << "Not-hidden/total = " << atomA.getNumberOfAltConformations() << "/" << atomA.getNumberOfAltConformations(true) << endl;
		cout << atomA << endl;

	}

	return 0;


}

