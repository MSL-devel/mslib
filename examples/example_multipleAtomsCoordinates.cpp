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
#include <cstdlib>

#include "Atom.h"
#include "MslTools.h"

using namespace std;
using namespace MSL;

/*******************************************************************
 *  This program illustrates how to add alternative coordinates to
 *  an Atom and how to switch the active one
 * 
 *******************************************************************/

int main() {

	cout << "  ***************************************************************************************" << endl;
	cout << "" << endl;
	cout << "     Example on setting multiple atom coordinates (" << MslTools::getMSLversion() << ")   " << endl;
	cout << "" << endl;
	cout << "  ***************************************************************************************" << endl;
	cout << endl;
	cout << endl;

	// Create an atom named CA and put it at the origin
	Atom A("CA", 0.0, 0.0, 0.0);

	cout << "Atom A has " << A.getNumberOfAltConformations() << " conformations and the active conformation has index " << A.getActiveConformation() << endl;
	
	// atoms can be printed!
	cout << A << endl;
	cout << endl;
	cout << "=============================" << endl;

	// add four alternative positions (the first one will remain the active one)
	A.addAltConformation(1.0, 0.0, 0.0);
	A.addAltConformation(2.0, 0.0, 0.0);
	A.addAltConformation(3.0, 0.0, 0.0);
	A.addAltConformation(4.0, 0.0, 0.0);

	cout << "Added 4 more conformations to the atoms" << endl;
	cout << "Atom A has " << A.getNumberOfAltConformations() << " conformations and the active conformation has index " << A.getActiveConformation() << endl;
	cout << A << endl;
	cout << endl;
	cout << "=============================" << endl;

	// change the active conformation and print the atoms
	for (unsigned int i=0; i<A.getNumberOfAltConformations(); i++) {
		cout << "Set atom in the " << i << "-th conformation and print the atom" << endl;
		A.setActiveConformation(i);
		cout << A << endl;
		cout << endl;
	}


	return 0;
}
