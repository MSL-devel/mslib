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


#include "ResidueSelection.h"

using namespace MSL;
using namespace std;


ResidueSelection::ResidueSelection(){
	sys = NULL;
	debug = false;
	storedSelections.clear();
}

ResidueSelection::ResidueSelection(System &_sys){
	sys = &_sys;
	debug = false;
	storedSelections.clear();
}

ResidueSelection::~ResidueSelection(){
}


bool ResidueSelection::selectionExists(string _selectName){
	Hash<string,vector<Residue *> >::Table::iterator it = storedSelections.find(_selectName);
	return (it != storedSelections.end());
}
vector<Residue *> & ResidueSelection::getSelection(string _selectName){
	vector<Residue *>  *a = NULL;
	Hash<string,vector<Residue *> >::Table::iterator it = storedSelections.find(_selectName);
	if (it != storedSelections.end()){
		a = &(it->second);
	}

	return (*a);
}
vector<Residue *> & ResidueSelection::select(string _selectString){
	
	// Tokenize select string based on ',' character 
	vector<string> selectToks = MslTools::tokenize(_selectString, ",");

	string name         = MslTools::toUpper(selectToks[0]);
	string selectString = selectToks[0];
	if (selectToks.size() > 1){
		selectString = selectToks[1];
	}


	// Check for no data
	if (sys == NULL) {
		cerr << "ERROR no data to select from, most likely you forgot to give this ResidueSelection object a starting vector<Residue *>  in its constructor."<<endl;
		exit(-1);
	}


	vector<Residue *>  tmp;


	// Parse the Logic
	if (debug){
		cout << "Selection Statement: "<<selectString<<endl;
	}
	LogicalParser lp;
	lp.setDebugFlag(debug);
	lp.setLogicStatementInFix(selectString);
	lp.parse();

	if (debug){
		lp.printLogicTree();
	}

	for (uint i = 0;i < sys->positionSize();i++){

		Residue &r = sys->getResidue(i);
		if (lp.eval(r)){
			r.setSelectionFlag(name, true);
			tmp.push_back(&r);
		} else {
			r.setSelectionFlag(name, false);
		}
	}


	storedSelections[name] = tmp; // this does a copy
	vector<Residue *>  *a =  &storedSelections[name];


	return (*a);
}


