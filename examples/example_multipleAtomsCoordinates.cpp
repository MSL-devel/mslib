#include <iostream>
#include <cstdlib>

#include "Atom.h"
#include "MslTools.h"

using namespace std;

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
