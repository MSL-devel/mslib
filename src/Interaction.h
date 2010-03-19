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

#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>

#include "Atom.h"

/* ERROR CODE 34xxx */


namespace MSL { 
class Interaction {
	public:
		virtual ~Interaction();

		virtual Atom * operator[](size_t _n);
		virtual double operator()(size_t _n);
		
		std::vector<Atom*> & getAtomPointers();
		std::vector<double> & getParams();
		Atom * getAtom(size_t _n) const;
		double getParam(size_t _n) const;
		unsigned int getAtomSize() const;
		unsigned int getParamSize() const;
		void setAtoms(std::vector<Atom*> _atoms);
		void setParams(std::vector<double> _params);

		virtual bool isSelected(std::string _sele1, std::string _sele2) const=0;
		virtual bool isActive() const=0;
		virtual double getEnergy()=0;
		virtual double getEnergy(double _param)=0;
		
		// print atom information
		virtual std::string toString() const=0;
	//	virtual unsigned int getType() const=0;
		virtual std::string getName() const=0;

		bool atomsHaveCoordinates() const;

		void update();
	
	protected:
		/*************************************************
		 *   The following enum is for the selection type
		 *   - none: no atoms, always selected
		 *   - single: one atom, if it is in both selections
		 *   - distance: two atoms, one in one selection, the other in the other selection
		 *   - angle: 3 atoms, center atoms in both selections
		 *   - dihedral: 4 atoms, the second and third atoms in different selections
		 *   - improper: 4 atoms, the first (improper center) atom in both selections
		 *************************************************/
		enum SelectionTypes { none=0, single=1, distance=2, angle=3, dihedral=4, improper=5 };
		Interaction();
		std::vector<Atom*> pAtoms;
		std::vector<double> params;
		double energy;

};

inline std::vector<Atom*> & Interaction::getAtomPointers() {return pAtoms;}
inline std::vector<double> & Interaction::getParams() {return params;}
inline Atom * Interaction::getAtom(size_t _n) const {return pAtoms[_n];}
inline double Interaction::getParam(size_t _n) const {return params[_n];}
inline unsigned int Interaction::getAtomSize() const {return pAtoms.size();};
inline unsigned int Interaction::getParamSize() const {return params.size();}
inline Atom * Interaction::operator[](size_t _n) {return pAtoms[_n];}
inline double Interaction::operator()(size_t _n) {return params[_n];}
inline void Interaction::setAtoms(std::vector<Atom*> _atoms) {pAtoms = _atoms;}
inline void Interaction::setParams(std::vector<double> _params) {params = _params;}
inline void Interaction::update() {} // emtpy function, some terms, like charmm elec might need to update
inline bool Interaction::atomsHaveCoordinates() const {
	for (std::vector<Atom*>::const_iterator k=pAtoms.begin(); k!=pAtoms.end(); k++) {
		if (*k == NULL || !(*k)->hasCoor()) {
			return false;
		}
	}
	return true;
}

}

#endif

