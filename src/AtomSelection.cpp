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

#include "AtomSelection.h"

AtomSelection::AtomSelection(){
	data = NULL;
	debug = false;
	storedSelections.clear();
}

AtomSelection::AtomSelection(AtomPointerVector &_data){
	data = &_data;
	debug = false;
	storedSelections.clear();
}

AtomSelection::~AtomSelection(){
}


bool AtomSelection::selectionExists(string _selectName){
	Hash<string,AtomPointerVector>::Table::iterator it = storedSelections.find(_selectName);
	return (it != storedSelections.end());
}
AtomPointerVector& AtomSelection::getSelection(string _selectName){
	AtomPointerVector *a = NULL;
	Hash<string,AtomPointerVector>::Table::iterator it = storedSelections.find(_selectName);
	if (it != storedSelections.end()){
		a = &(it->second);
	}

	return (*a);
}
AtomPointerVector& AtomSelection::select(string _selectString, bool _selectAllAtoms){
	
	// Tokenize select string based on ',' character 
	vector<string> selectToks = MslTools::tokenize(_selectString, ",");

	string name         = MslTools::toUpper(selectToks[0]);
	string selectString = selectToks[0];
	if (selectToks.size() > 1){
		selectString = selectToks[1];
	}

	AtomPointerVector *a = NULL;

	// Check for no data
	if (data == NULL) {
		cerr << "ERROR no data to select from, most likely you forgot to give this AtomSelection object a starting AtomPointerVector in its constructor."<<endl;
		exit(-1);
	}


	// Create AtomPointerVector for data
	AtomPointerVector tmp;

	/*
	  Check for complex selections:  
	  1.   SEL1 within 5 of SEL2 
	*/


	int wpos = selectString.find("WITHIN",0);
	int opos = selectString.find("OF",0);
	if (wpos != std::string::npos && opos != std::string::npos) {

		// Complex Selection form 'SEL1 WITHIN X of SEL2'
		//   Create 2 logical parsers.  Create two tmp AtomPointerVectors.  Do combination.
			
		string sel1 = MslTools::trim(selectString.substr(0,wpos));
		double dist = MslTools::toDouble(MslTools::trim(selectString.substr(wpos+6,wpos-opos)),"AtomSelection::select is trying to find a distance in selection statement");
		string sel2 = MslTools::trim(selectString.substr(opos+2));

		//cout << "sel1: '"<<sel1<<"'"<<endl;
		//cout << "dist: '"<<dist<<"'"<<endl;
		//cout << "sel2: '"<<sel2<<"'"<<endl;

		// Parse the Logic
		LogicalParser lp1;
		lp1.setDebugFlag(debug);
		lp1.setLogicStatementInFix(sel1);
		lp1.parse();

		LogicalParser lp2;
		lp2.setDebugFlag(debug);
		lp2.setLogicStatementInFix(sel2);
		lp2.parse();


		// Locally store sel1,sel2
		AtomPointerVector tmp1;
		AtomPointerVector tmp2;
		AtomPointerVector::iterator avIt;
		for (avIt = data->begin();avIt != data->end();avIt++){

			if (lp1.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())){
				tmp1.push_back(*avIt);
			}

			if (lp2.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())) {
				tmp2.push_back(*avIt);
			}

		}

		AtomPointerVector::iterator avIti;
		AtomPointerVector::iterator avItj;

		for (avIti = tmp1.begin(); avIti != tmp1.end(); avIti++){

			bool withinDistance = false;
			for (avItj = tmp2.begin(); avItj != tmp2.end(); avItj++){

				if ((*avIti)->distance((**avItj)) <= dist) {
					withinDistance = true;
					break;
				}
			}

			if (withinDistance){
				(*avIti)->setSelectionFlag(name, true);
				tmp.push_back(*avIti);
			} else {
				(*avIti)->setSelectionFlag(name, false);
			}
				
		}

	} else {


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


		AtomPointerVector::iterator avIt;
		for (avIt = data->begin();avIt != data->end();avIt++){
			//cout  << "Data: "<<**avIt<<endl;

			if (lp.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())) {
				(*avIt)->setSelectionFlag(name, true);
				tmp.push_back(*avIt);
			} else {
				(*avIt)->setSelectionFlag(name, false);
			}

		}
	}


	storedSelections[name] = tmp; // this does a copy
	a =  &storedSelections[name];
		



	return (*a);
}


