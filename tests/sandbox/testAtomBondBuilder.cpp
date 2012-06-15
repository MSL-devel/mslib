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

#include "AtomPointerVector.h"
#include "AtomBondBuilder.h"
#include "PDBReader.h"
#include "System.h"

using namespace std;

using namespace MSL;


int main(int argc,char *argv[]) {

	/*
	PDBResysr rAv;
	rAv.open(argv[1]);
	rAv.read();
	*/
//	if (argc == 1) {
//		cout << "Please specify a PDB file" << endl;
//		exit(0);
//	}
	string file = "exampleFiles/example0001.pdb";
	System sys;
	if (!sys.readPdb(file)) {
		cerr << "Cannot read pdb file " << file << endl;
		exit(1);
	}
	AtomPointerVector av = sys.getAtomPointers();
	//cout << av;

	AtomBondBuilder abb;
	abb.buildConnections(av);
	cout << "==========================" << endl;
	for (unsigned int i=0; i<av.size(); i++) {
		char c [1000];
		sprintf(c, "%04d", i);
		cout << c << " " << *av[i] << "  : ";
		for (unsigned int j=0; j<av.size(); j++) {
			if (i == j) {
				cout << "1";
			} else if (av[i]->isBoundTo(av[j])) {
				cout << "2";
			} else if (av[i]->isOneThree(av[j])) {
				cout << "3";
			} else if (av[i]->isOneFour(av[j])) {
				cout << "4";
			} else {
				cout << "-";
			}
		}
		cout << endl;
	}

	cout << "==========================" << endl;

	for (unsigned int i=0; i<av.size(); i++) {
		cout << av[i]->getAtomId() << " is bound to" << endl;
		vector<Atom*> bonded = av[i]->getBonds();
		for (unsigned int j=0; j<bonded.size(); j++) {
			cout << "   " << bonded[j]->getAtomId() << endl;
		}
	}
	cout << "==========================" << endl;

	av[6]->setUnboundFrom(av[5]);	// 1N 2CA

	cout << "==========================" << endl;
	cout << "Let's now manually remove some bonds" << endl;
	cout << "Removed bond between atoms 1N-2CA" << endl;
	for (unsigned int i=0; i<av.size(); i++) {
		char c [1000];
		sprintf(c, "%04d", i);
		cout << c << " " << *av[i] << "  : ";
		for (unsigned int j=0; j<av.size(); j++) {
			if (i == j) {
				cout << "1";
			} else if (av[i]->isBoundTo(av[j])) {
				cout << "2";
			} else if (av[i]->isOneThree(av[j])) {
				cout << "3";
			} else if (av[i]->isOneFour(av[j])) {
				cout << "4";
			} else {
				cout << "-";
			}
		}
		cout << endl;
	}

	av[7]->setUnboundFrom(av[6]);	// 2CA 2CB

	cout << "==========================" << endl;
	cout << "Removed bond between atoms 2CA-2CB" << endl;
	for (unsigned int i=0; i<av.size(); i++) {
		char c [1000];
		sprintf(c, "%04d", i);
		cout << c << " " << *av[i] << "  : ";
		for (unsigned int j=0; j<av.size(); j++) {
			if (i == j) {
				cout << "1";
			} else if (av[i]->isBoundTo(av[j])) {
				cout << "2";
			} else if (av[i]->isOneThree(av[j])) {
				cout << "3";
			} else if (av[i]->isOneFour(av[j])) {
				cout << "4";
			} else {
				cout << "-";
			}
		}
		cout << endl;
	}

	av[11]->setUnboundFrom(av[6]);	// 2C 2CA

	cout << "==========================" << endl;
	cout << "Removed bond between atoms 2CA-2C" << endl;
	for (unsigned int i=0; i<av.size(); i++) {
		char c [1000];
		sprintf(c, "%04d", i);
		cout << c << " " << *av[i] << "  : ";
		for (unsigned int j=0; j<av.size(); j++) {
			if (i == j) {
				cout << "1";
			} else if (av[i]->isBoundTo(av[j])) {
				cout << "2";
			} else if (av[i]->isOneThree(av[j])) {
				cout << "3";
			} else if (av[i]->isOneFour(av[j])) {
				cout << "4";
			} else {
				cout << "-";
			}
		}
		cout << endl;
	}

	av[11]->setUnboundFrom(av[13]);	// 2C 3N

	cout << "==========================" << endl;
	cout << "Removed bond between atoms 2C-3N" << endl;
	for (unsigned int i=0; i<av.size(); i++) {
		char c [1000];
		sprintf(c, "%04d", i);
		cout << c << " " << *av[i] << "  : ";
		for (unsigned int j=0; j<av.size(); j++) {
			if (i == j) {
				cout << "1";
			} else if (av[i]->isBoundTo(av[j])) {
				cout << "2";
			} else if (av[i]->isOneThree(av[j])) {
				cout << "3";
			} else if (av[i]->isOneFour(av[j])) {
				cout << "4";
			} else {
				cout << "-";
			}
		}
		cout << endl;
	}


}

