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


#include "CRDReader.h"

using namespace MSL;
using namespace std;



/**
 * This method will delete all of the Atom pointers
 * held in this CRDReader.
 */
void CRDReader::deletePointers() {
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}


	atoms.clear();
	if(pTopReader) {
		delete pTopReader;
		pTopReader = NULL;
	}
}


bool CRDReader::read() {
	if (!is_open()) {
		cerr << "WARNING: 14567 CRDReader::read () File not open" << endl;
		return false;
	}

	try { 
		string currentResidue = "";

		bool doneWithTitle = false;
		while (!endOfFileTest()){
			string line = Reader::getLine();
			if(!doneWithTitle) {
				if(line.substr(0,1)=="*") {
					// Skip comments
					continue;
				} else {
					// the first line after the title is the number of atoms
					doneWithTitle = true;	
					continue;
				}
			}
			if(line.size() < CRDFormat::E_CHARGE-1) {
				continue;
			}
						
			CRDFormat::AtomData atom = CRDFormat::parseAtomLine(line);

			atoms.push_back(new Atom);

			// atom name, residue name, residue icode, chain id, coor, element
			atoms.back()->setName(atom.D_ATOM_NAME);
			atoms.back()->setResidueName(atom.D_RES_NAME);
			atoms.back()->setResidueNumber(atom.D_RES_NUM);
			atoms.back()->setChainId(atom.D_CHAIN_ID);
			atoms.back()->setCoor(atom.D_X,atom.D_Y, atom.D_Z);
//			Dont know what the element type is
			if(pTopReader) {
				// use the information in the Topology File
				atoms.back()->setElement(pTopReader->getElement(std::string(atom.D_ATOM_NAME)));
			} else {
				// default to first letter of the atomType
				atoms.back()->setElement(std::string(atom.D_ATOM_NAME,1));
			}
			atoms.back()->setCharge(atom.D_CHARGE);
		
		}
		return true;
		

	} catch(...){
		cerr << "ERROR 5623 in CRDReader::read()\n";
		return false;
	}

	return true;
}
