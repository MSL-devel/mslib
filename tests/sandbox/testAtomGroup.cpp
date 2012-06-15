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
#include "AtomGroup.h"

using namespace std;

using namespace MSL;


int main() {

	AtomGroup AG;

	AG.push_back(new Atom);
	AG.push_back(new Atom);
	AG.push_back(new Atom);

	cout << "Create an AtomGroup, push_back 3 atoms and set the coordinates to:" << endl;
	cout << "Atom 0          : 2.3,  5.4,  2.1" << endl;
	cout << "Atom 0          : 5.2, -2.1,  7.3" << endl;
	cout << "Atom 0          : 9.3,  6.0, -4.3" << endl;
	cout << "                  ===============" << endl;
	cout << "Geometric center: 5.6,  3.1,  1.7" << endl;
	cout << endl;

	AG[0]->setCoor(2.3,  5.4,  2.1);
	AG[1]->setCoor(5.2, -2.1,  7.3);
	AG[2]->setCoor(9.3,  6.0, -4.3);

	cout << "The atom group with address " << &AG << " has size " << AG.size() << endl;
	cout << endl;

	cout << "Print the coordinates and the address of the parent group" << endl;
	for (AtomGroup::iterator k=AG.begin(); k!=AG.end(); k++) {
		cout << **k << " " << (*k)->getParentGroup() << endl;
	}
	cout << endl;
	cout << "Geometric center: " << AG.getGeometricCenter() << endl;

}
