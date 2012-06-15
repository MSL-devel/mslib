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
		bool hasAtom(Atom * _pAtom) const;

		virtual bool isSelected(std::string _sele1, std::string _sele2) const=0;
		virtual bool isActive() const=0;
		virtual double getEnergy()=0;

		virtual std::pair<double,std::vector<double> > partialDerivative() = 0; // computes d(_param)/dx1,d(_param)/dy1,.. etc..,
		virtual double getEnergy(double _param, std::vector<double> *paramDerivatives=NULL)=0; // computes dE/d(_param)

		// this function is used by the minimizer. It computes energy without the switching function even if cutoffs are in place
		virtual double getEnergy(std::vector<double> *_paramDerivatives)=0; // computes dE/dx1,dE/dy1,dE/dz1.....and stores in _paramDerivatives and returns energy
		virtual std::vector<double> getEnergyGrad()=0; // computes  and returns dE/dx1,dE/dy1,dE/dz1....


		virtual bool reset();
		
		// print atom information
		virtual std::string toString() =0;
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
inline bool Interaction::reset() {
	std::cerr << "Reset not Implemented for this interaction type" << std::endl;
	return false;
}
inline bool Interaction::hasAtom(Atom * _pAtom) const {
	for (std::vector<Atom*>::const_iterator k=pAtoms.begin(); k!=pAtoms.end(); k++) {
		if (_pAtom == *k) {
			return true;
		}
	}
	return false;
}

}

#endif

