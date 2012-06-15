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

#include "CRDWriter.h"

using namespace MSL;
using namespace std;


bool CRDWriter::write(AtomPointerVector &_av, bool _writeRemarks) {

	if( is_open() == false )
	return false;

	if(_writeRemarks) {
		writeREMARKS();
	}

	//The  first line is the number of atoms
	char tmp[10];
	if(_av.size() <= 99999) {
		sprintf(tmp,"%5lu",_av.size());
	} else {
		sprintf(tmp,"99999");
	}

	string numAtoms = string(tmp);
	writeln(numAtoms);

	string prevPosId = "";
	string crdline = "";

	unsigned int absres = 0;
	unsigned int atomNum = 0;

	for (AtomPointerVector::iterator it = _av.begin(); it != _av.end(); it++){
		if (atomNum < 99999) {
			atomNum++;
		}
		string thisPosId = (*it)->getPositionId();
		if (it == _av.begin() || thisPosId != prevPosId) {
			// absolute residue number, from 1 to N to the end of the file, independently from residue numbers and chains
			absres++;
			prevPosId = thisPosId;
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


void CRDWriter::writeREMARKS(){

	if (fileHandler != stringstyle && !fileStream.is_open()){
		cout << "File : "<<fileName<< " is not open for writing\n";
		return;
	}

	vector<string>::iterator it;
	for (it = remarks.begin(); it != remarks.end(); it++){
		vector<string> singleRemarks = MslTools::tokenize(*it, "\n");
		for (uint i = 0;i < singleRemarks.size();i++){
				string line = "*\t"+singleRemarks[i];
				writeln(line);
		}
	}
}
