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

#ifndef CONFORMATIONEDITOR_H
#define CONFORMATIONEDITOR_H
/*

 */

//MSL Includes
#include "System.h"
#include "MslTools.h"
#include "DegreeOfFreedomReader.h"

/**
 * This class is an object for altering the conformation of a protein.
 * 
 */
namespace MSL { 
class ConformationEditor {

	public:
		// Constructors/Destructors
		ConformationEditor();
		ConformationEditor(System & _sys);
		~ConformationEditor();

		bool readDefinitionFile(std::string _defiFile); 

		void applyPDBdefinitions();

		// the following take a positionId ("A,37") and a definition
		// such as "N,CA" (bond), "N,CA,CB" (angle), "N,CA,CB,CG" (dihedral),
		// or with a label "chi1", or "phi".  The labels need to be defined
		// normally by reading a definition file
		bool editIC(std::string _positionId, std::string _deegreOfFreedom, double _value);
		// this one takes 2-4 atomIds ("A,27,N A,27,CA A,27,CB A,27,CG")
		bool editIC(std::string _deegreOfFreedom, double _value);
		// this one takes the atom pointers
		bool editIC(vector<Atom*> _pAtoms, double _value);

		bool applyConformation();

	private:

		System * pSys;
		IcTable * pIcTable;

		set<Atom*> buildAtoms;

		std::map<std::string, std::map<std::string, vector<std::string> > > degOfFreedomlabels;

		void defineDegreesOfFreedom();
};

inline bool ConformationEditor::applyConformation() {
	bool out = true;
	for (set<Atom*>::iterator k=buildAtoms.begin(); k!=buildAtoms.end(); k++) {
		if (!(*k)->buildFromIc()) {
			cerr << "WARNING 34823: failed building " << **k << endl;
			out = false;
		}
	}
	return out;
}

//Inlines go HERE





}
#endif
