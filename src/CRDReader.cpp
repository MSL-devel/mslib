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


#include "CRDReader.h"

using namespace MSL;
using namespace std;



/**
 * This method will delete all of the Atom pointers
 * held in this CRDReader.
 */
void CRDReader::deletePointers() {
	if(pTopReader != NULL) {
		delete pTopReader;
	}
	pTopReader = NULL;
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}


	atoms.clear();
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
			string chainIdwaste = "";
			int resNum = 0;
			string icode = "";
			MslTools::parsePositionId(atom.D_RES_NUM, chainIdwaste, resNum, icode, 1);
			atoms.back()->setResidueNumber(resNum);
			atoms.back()->setResidueIcode(icode);
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
