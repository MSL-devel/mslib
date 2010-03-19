/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
	if (argc == 1) {
		cout << "Please specify a PDB file" << endl;
		exit(0);
	}
	System sys;
	if (!sys.readPdb(argv[1])) {
		cerr << "Cannot read pdb file " << argv[1] << endl;
		exit(1);
	}
	AtomPointerVector av = sys.getAtomPointers();
	//cout << av;

	AtomBondBuilder abb;
	abb.buildConnections(av);

	for (unsigned int i=0; i<av.size(); i++) {
		cout << *av[i] << endl;
		vector<Atom*> bonded = av[i]->getBonds();
		for (unsigned int j=0; j<bonded.size(); j++) {
			cout << "   " << *bonded[j] << endl;
		}
	}


}

