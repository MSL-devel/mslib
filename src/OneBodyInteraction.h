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


#ifndef ONEBODYINTERACTION_H
#define ONEBODYINTERACTION_H

#include <iostream>
#include <vector>
#include <string>

#include "Interaction.h"


namespace MSL { 
	class OneBodyInteraction: public Interaction {

		/*******************************************************
		 *   Inherits from Interaction (a prototype object
		 *      AtomPointerVector pAtoms  (the atoms) 
		 *      std::vector<double> param  (the parameters
		 *******************************************************/

		public:

			virtual ~OneBodyInteraction();

			/* setting and getting the atoms */
			void setAtoms(std::vector<Atom*> _atoms);
			void setAtoms(Atom & _a1);

			bool isSelected(std::string _selection1, std::string _selection2) const;
			bool isActive() const;

			virtual double getEnergy()=0;
			virtual double getEnergy(double _param, std::vector<double> *paramDerivatives=NULL)=0;
			virtual std::vector<double> getEnergyGrad()=0;

			friend std::ostream & operator<<(std::ostream &_os, OneBodyInteraction & _term) {_os << _term.toString(); return _os;};
			virtual std::string toString() =0;

		protected:
			OneBodyInteraction();

	};


	inline void OneBodyInteraction::setAtoms(std::vector<Atom*> _atoms) { if (_atoms.size() != 1) {std::cerr << "ERROR 37967: invalid number of atoms in inline void OneBodyInteraction::setAtoms(std::vector<Atom*> _atoms)" << std::endl; exit(37967);} pAtoms = _atoms;}
	inline void OneBodyInteraction::setAtoms(Atom & _a1) { pAtoms[0] = &_a1;}
	inline bool OneBodyInteraction::isSelected(std::string _selection1, std::string _selection2) const {
		if (pAtoms[0]->getSelectionFlag(_selection1) && pAtoms[0]->getSelectionFlag(_selection2)) {
			return true;
		} else {
			return false;
		}
	}

	inline bool OneBodyInteraction::isActive() const {return pAtoms[0]->getActive();}


}

#endif

