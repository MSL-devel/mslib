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

#include "OptionParser.h"
#include "Timer.h"
#include "MIDReader.h"
#include "MoleculeInterfaceDatabase.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>


using namespace std;

using namespace MSL;


int main(int argc, char *argv[]){
	vector<string> required;
	required.push_back("flatfile");
	
	vector<string> optional;
	optional.push_back("reread");
	optional.push_back("outfile");
	OptionParser OP;
	OP.readArgv(argc,argv);
	OP.setRequired(required);
	OP.setAllowed(optional);

	string flatfile = OP.getString("flatfile");
	if (OP.fail()){
		cout << "Usage: ./testMolecularInterfaceDatabase --flatfile FLAT.mid [ --outfile mid.ckpt --reread --archive text ]\n";
		exit(0);
	}

	string outfile = OP.getString("outfile");
	if (OP.fail()){
		outfile = "mid.ckpt";
	}

	string archive = OP.getString("archive");
	if (OP.fail()){
		archive = "text";
	}
	
	// Create from flat file
	MIDReader mread;
	mread.open(flatfile);
	mread.read(archive);
	
	MoleculeInterfaceDatabase &mid = mread.getMoleculeInterfaceDatabase();
	cout << "MID size: "<<mid.size()<<endl;
	
	// Write out to binary 
	mid.setArchiveType(archive);
	mid.save_checkpoint(outfile);

	// Re-read if option set..
	bool reread = OP.getBool("reread");
	if (!OP.fail() && reread){
		cout << "Re-reading binary MolecularInterfaceDatabase.\n";
		Timer t;
		double start = t.getWallTime();
		MoleculeInterfaceDatabase mid2;
		mid2.setArchiveType(archive);
		mid2.load_checkpoint(outfile);
		cout << "MID2 size: "<<mid2.size()<<" read in "<<(t.getWallTime() - start)<<" seconds."<<endl;
	}
}
