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


#include "PDBWriter.h"
#include "PDBReader.h"
#include "EnergySet.h"
#include "testData.h"
#include "AtomSelection.h"
#include "CharmmVdwInteraction.h"
#include "CharmmBondInteraction.h"

using namespace std;

using namespace MSL;


int main(){

	// Coiling an ideal, Z-aligned helix
	PDBReader pin(idealHelix);
	AtomPointerVector ideal;
	ideal = pin.getAtomPointers();
	pin.close();

	string filename = "/tmp/out.pdb";
	PDBWriter writer(filename);
    writer.open();
	writer.write(ideal);
	writer.close();
	cout << endl;
	cout << "=========================" << endl;
	cout << "Written output pdb " << filename << endl;

	EnergySet ESet;	
	CharmmBondInteraction * pCBI= new CharmmBondInteraction(*ideal[0], *ideal[1], 1, 2);
	cout << "Created bond with parameters (1, 2) between atoms " << *ideal[0] << *ideal[1] << endl;
	cout << endl;
	cout << "Distance = " << ideal[0]->distance(*ideal[1]) << endl;
	cout << "Energy = " << pCBI->getEnergy() << endl;
	ESet.addInteraction(pCBI);
	//ESet.calcEnergy();
	//ESet.printSummary();

	AtomSelection selector(ideal);
	AtomPointerVector res7 = selector.select("res7, resi 7 and chain A");
	cout << "Selected res7 with " << res7.size() << " atoms" << endl;
	AtomPointerVector res8 = selector.select("res8, resi 8 and chain A");
	cout << "Selected res8 with " << res8.size() << " atoms" << endl;
	CharmmVdwInteraction * pCVI= new CharmmVdwInteraction(*res7[0], *res7[4], 4, -2);
	ESet.addInteraction(pCVI);
	ESet.calcEnergy();
	ESet.printSummary();
	cout << "Total number of interactions = " << ESet.getTotalNumberOfInteractionsCalculated() << endl;

	EnergySet ESet2(ESet); //Copy doesnt work as of now
	ESet2.calcEnergy();
	ESet2.printSummary();
}
