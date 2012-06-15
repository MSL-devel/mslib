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

#include "EnvironmentDatabase.h"
#include "AtomPointerVector.h"
#include "AtomSelection.h"
#include "PDBReader.h"
#include "PDBWriter.h"
#include "testData.h"

#include <string>

using namespace std;

using namespace MSL;



int main() {
	       
		
	// Create a four Helix Bundle 
	AtomPointerVector av;
	stringstream ss;
	ss.str(fourHelixBundle);

	PDBReader rAv(ss);
	rAv.read();
	av = rAv.getAtomPointers();
	rAv.close();


	// Setup an Environment Database
	System sys(av);
	EnvironmentDatabase ed("LEU");
	ed.createDatabase(sys, "testDatabase");


	cout << "Number of descriptors is: "<<ed.getNumberDescriptors()<<endl;

	// Save as a ckpt file
    //	ed.save_checkpoint("envData.ckpt");
	
	EnvironmentDatabase ed2;

	// Load as a ckpt file
    //	ed2.load_checkpoint("envData.ckpt");
	cout << "Number of descriptors is: "<<ed2.getNumberDescriptors()<<endl;



	// Use loaded ckpt file to search
	vector<EnvironmentDescriptor *> allED = ed2.getAllDescriptors();

	if (ed2.searchForEnvironment(*allED[0],"CA")){
		vector<EnvironmentDescriptor *> foundED = ed2.getSearchResults();

		cout << "Num ED Results: "<<foundED.size()<<endl;

		// Print out envFrame, refFrame?
		
	} else {
		cout << "No search results were found."<<endl;
	}
	

}
