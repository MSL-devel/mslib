/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#include "CRDWriter.h"

using namespace MSL;
using namespace std;


bool CRDWriter::write(AtomPointerVector &_av) {

	if( is_open() == false )
	return false;

	string prevPosId = "";
	unsigned int absres = 0;
	unsigned int atomNum = 0;
	string crdline = "*                                                                               ";
	if (!writeln(crdline)) {
		cerr << "WARNING 12487: cannot write atom line in bool CRDWriter::write(AtomPointerVector &_av)" << endl;
		return false;
	}
	for (AtomPointerVector::iterator it = _av.begin(); it != _av.end(); it++){
		if (atomNum < 99999) {
			atomNum++;
		}
		if (it == _av.begin() || (*it)->getPositionId() != prevPosId) {
			// absolute residue number, from 1 to N to the end of the file, independently from residue numbers and chains
			absres++;
		}

		CRDFormat::AtomData atom;

		atom = CRDFormat::createAtomData(**it);
		
		crdline = CRDFormat::createAtomLine(atom, atomNum, absres);

		if (!writeln(crdline)) {
			cerr << "WARNING 12491: cannot write atom line in bool CRDWriter::write(AtomPointerVector &_av)" << endl;
			return false;
		}
		
	}

	return true;
}
